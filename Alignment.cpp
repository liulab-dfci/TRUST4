#include <map>

#include "Alignment.h"

extern int radius ;
extern char numToNuc[26] ;

int CompPairBInc( const void *p1, const void *p2 )
{
	return ( ( struct _pair *)p1 )->b - ( ( struct _pair *)p2 )->b ;
}

void ReverseComplement( char *s, int len )
{
	int i, j ;
	for ( i = 0, j = len - 1 ; i < j ; ++i, --j )
	{
		char tmp ;
		tmp = s[j] ;
		s[j] = s[i] ;
		s[i] = tmp ;
	}
	for ( i = 0 ; i < len ; ++i )
	{
		s[i] = numToNuc[ 3 - (int)nucToNum[ s[i] - 'A' ] ] ;
	}
}

bool IsSameStrand( char *sa, index_t pa, char *sb, index_t pb, int kmerLength )
{
	int i, j ;
	// Since we make sure kmer length is odd, it can't be palidrome
	for ( i = pa, j = pb ; i < pa + kmerLength  ; ++i, ++j )	
	{
		if ( sa[i] != sb[j] )
			return false ;
	}
	return true ;	
}

// Return the first index whose hits.b is smaller or equal to valB
int BinarySearch_LIS( int top[], int size, int valB, SimpleVector<struct _pair> &hits )
{
	int l = 0, r = size - 1 ;
	int m ;
	while ( l <= r )
	{
		m = ( l + r ) / 2 ;
		if ( valB == hits[ top[m] ].a )
		{
			return m ;
		}
		else if ( valB < hits[ top[m] ].a )
		{
			r = m - 1 ;
		}
		else
		{
			l = m + 1 ;
		}
	}
	return l - 1 ;
}

// The O(nlogn) method for solving LIS problem, suppose there are n hits.
// Return the LIS, the LIS's length is returned by the function
int LongestIncreasingSubsequence( SimpleVector<struct _pair> &hits, struct _pair LIS[] ) 
{
	// Only use the first hit of each qhit
	// Bias towards left

	int i, j, k ;
	int ret = 0 ;
	int size = hits.Size() ;
	
	int *record = new int[size] ; // The index of the selected hits
	int *top = new int[size] ; // record the index of the hits with smallest valB of the corresponding LIS length. kind of the top element.
	int *link = new int[size] ; // used to retrieve the LIS
	
	int rcnt = 1 ;
	record[0] = 0 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( hits[i].b == hits[i - 1].b )
			continue ;
		record[rcnt] = i ;
		++rcnt ;
	}
	top[0] = 0 ;
	link[0] = -1 ;
	ret = 1 ;
	for ( i = 1 ; i < rcnt ; ++i )
	{
		int tag = BinarySearch_LIS( top, ret, hits[ record[i] ].a, hits ) ;				
		if ( tag == -1 )
		{
			top[0] = i ;
			link[i] = -1 ;
		}
		else if ( hits[ record[i] ].a > hits[ top[tag] ].a )
		{
			if ( tag == ret - 1 )
			{
				top[ret] = i ;
				++ret ;
				link[i] = top[tag] ;
			}
			else if ( hits[ record[i] ].a < hits[ top[tag + 1] ].a )
			{
				top[ tag + 1 ] = i ;
				link[i] = top[tag] ;
			}
		}
	}


	k = top[ret - 1] ;
	for ( i = ret - 1 ; i >= 0 ; --i )
	{
		//printf( "%d %d\n", i, k ) ;
		LIS[i] = hits[k] ;
		k = link[k] ;	
	}

	delete []top ;
	delete []record ;
	delete []link ;

	return ret ;
}


// Is two sequence a and b overlap?
// In general, a is for base read, b is for query, and in hits, hits.b is already sorted.
// b is always in the forward strand.
// Typically, identicalLength is the kmerLength+minimizer_windowSize-1
// rangeA holds the overlap for a; and rangeB for b.
bool IsOverlap( SimpleVector<struct _pair> &hits, bool sameStrand, int minHitRequired, struct _pair &rangeA, struct _pair &rangeB, SimpleVector<struct _pair> *retLIS ) 
{
	int len ;

	int i ;
	int hitCnt = hits.Size() ;
	if ( hitCnt < minHitRequired )
		return false ;

	// Pick a set of kmer hits such that they are kind of colinear.
	SimpleVector< struct _pair > hitCoordDiff ; // a-index of the hit, b-the difference.
	hitCoordDiff.Reserve( hitCnt ) ;

	for ( i = 0 ; i < hitCnt ; ++i )
	{
		struct _pair np ;
		np.a = i ;
		if ( sameStrand )
			np.b = hits[i].a - hits[i].b ;
		else
			np.b = hits[i].a + hits[i].b ;
		hitCoordDiff.PushBack( np ) ;
	}
	hitCoordDiff.QSort( CompPairBInc ) ;

	struct _pair bestWindow, currWindow ;
	bestWindow.a = bestWindow.b = -1 ;
	currWindow.a = currWindow.b = 0 ;
	int bestDiffSum = 0 ;
	int currDiffSum = 0 ;
	for ( i = 1 ; i < hitCnt ; ++i )
	{
		int diff = hitCoordDiff[i].b - hitCoordDiff[i - 1].b ;
		if ( diff < 0 )
			diff = -diff ;
		if ( diff > radius )
		{
			if ( currWindow.b - currWindow.a > bestWindow.a - bestWindow.b ||
					( currWindow.b - currWindow.a == bestWindow.a - bestWindow.b && 
					  currDiffSum < bestDiffSum ) )
			{
				bestWindow = currWindow ;
				bestDiffSum = currDiffSum ;
			}
			currDiffSum = 0 ;
			currWindow.a = currWindow.b = i ;
		}
		else
		{
			currDiffSum += diff ;
			currWindow.b = i ;
		}

	}
	// notice that the currWindow contains j-i-1, if it creates a new window but 
	// replaced the best window, then we have very small window and can filter this bundle.
	if ( currWindow.b - currWindow.a > bestWindow.a - bestWindow.b ||
			( currWindow.b - currWindow.a == bestWindow.a - bestWindow.b && 
			  currDiffSum < bestDiffSum ) )
	{
		bestWindow = currWindow ;
		bestDiffSum = currDiffSum ;
	}
	if ( bestWindow.b - bestWindow.a + 1 < minHitRequired )
		return false ;

	// Now, bestWindow gives the index of hit that are chosen, which LIS will be computed on.	
	int *diffRank = new int[ hitCnt ] ;
	SimpleVector<struct _pair> selectedHits ;
	selectedHits.Reserve( hitCnt ) ;

	for ( i = 0 ; i < hitCnt ; ++i )
		diffRank[ hitCoordDiff[i].a ] = i ;
	for ( i = 0 ; i < hitCnt ; ++i )
	{
		if ( diffRank[i] >= bestWindow.a && diffRank[i] <= bestWindow.b )
			selectedHits.PushBack( hits[i] ) ;
	}
	
	
	if ( !sameStrand )
	{
		int size = selectedHits.Size() ;
		for ( i = 0 ; i < size ; ++i )
			selectedHits[i].a = -selectedHits[i].a ;	
	}
	delete[] diffRank ;

	//for ( i = 0 ; i < hitCnt ; ++i )
	//	printf( "%d %d\n", selectedHits[i].a, selectedHits[i].b ) ;
	struct _pair *LIS = new struct _pair[ selectedHits.Size() ] ;

	int lengthLIS = LongestIncreasingSubsequence( selectedHits, LIS ) ;
	
	if ( !sameStrand )
		for ( i = 0 ; i < lengthLIS ; ++i )
		{
			LIS[i].a = -LIS[i].a ;
		}
	int LISstart, LISend ;
	LISstart = 0 ;
	LISend = lengthLIS - 1 ;
	//printf( "lengthLIS=%d\n", lengthLIS ) ;
	/*while ( LISstart + 2 * ( kmerLength - q + 1 ) < LISend )
	{
		len = ( LIS[ LISend ].a - LIS[ LISstart ].a + 1 + LIS[ LISend ].b - LIS[ LISstart ].b + 1 ) / 2 + q - 1 ;
		if ( LISend - LISstart + 1 < len * 0.02 ) // TODO: compute 0.02 in future.
		{
			// The LIS is too sparse.
			int len1 = ( LIS[ LISend - 1].a - LIS[ LISstart ].a + 1 
					+ LIS[ LISend - 1].b - LIS[ LISstart ].b + 1 ) / 2 + q - 1 ;
			int len2 = ( LIS[ LISend ].a - LIS[ LISstart + 1 ].a + 1 
					+ LIS[ LISend ].b - LIS[ LISstart + 1 ].b + 1 ) / 2 + q - 1 ;
			if ( len1 < len2 )
				--LISend ;
			else
				++LISstart ;
		}
		else
			break ;
		
	}*/
	bool ret = false ;
	if ( LISend - LISstart + 1 >= minHitRequired )
	{
		if ( sameStrand )
		{
			rangeA.a = LIS[ LISstart ].a ;
			rangeA.b = LIS[ LISend ].a ;
		}
		else
		{
			rangeA.a = LIS[ LISend ].a ;
			rangeA.b = LIS[ LISstart ].a ;
		}
	

		rangeB.a = LIS[ LISstart ].b ;
		rangeB.b = LIS[ LISend ].b ;

		ret = true ;
	}

	if ( retLIS != NULL )
	{
		retLIS->Reserve( LISend - LISstart + 1 ) ;
		for ( i = LISstart ; i <= LISend ; ++i )
		{
			retLIS->PushBack( LIS[i] ) ;
		}
	}

	delete []LIS ;
	return ret ;
}


// The implementation of Gene Myers bit-vector alignment method.
// Return: the number of different, and align holds the information on how the two strings aligned
// 0-match, 1-mismatch, 2-p has extra char(insert), 3-p misses a char(delete)
int GlobalAlignment( char *t, int lent, char *p, int lenp, char *align )
{
	int i, j ;
	int w = sizeof( uint64_t ) ;

	int bmax = lenp / w ;
	if ( bmax * w < lenp )
		++bmax ;

	int bid, boffset, blockSize ;
	
	blockSize = w * 8 ;
	
	uint64_t *Pv, *Mv, *Ph, *Mh ;
	struct _pair *computedBlockRange ;
	
	Pv = new uint64_t[ lent * bmax ] ;
	Mv = new uint64_t[ lent * bmax ] ;
	Ph = new uint64_t[ lent * bmax ] ;
	Mh = new uint64_t[ lent * bmax ] ;
	computedBlockRange = new struct _pair[ lent ] ;

	// Initialize Peq.
	uint64_t *Peq[4] ;
	for ( i = 0 ; i < 4 ; ++i )
	{
		Peq[i] = new uint64_t[ bmax] ;
		memset( Peq[i], 0, sizeof( uint64_t ) * bmax ) ;

		bid = 0 ;
		boffset = 0 ;
		for ( j = 0 ; j < lenp ; ++j )
		{
			if ( p[j] == numToNuc[i] )
				Peq[i][ bid ] |= ( (uint64_t)1 << boffset ) ;
			++boffset ;
			if ( boffset >= blockSize )
			{
				++bid ;
				boffset = 0 ;
			}
		}
		
		// padding
		while ( boffset < blockSize )
		{
			Peq[i][ bid ] |= ( (uint64_t)1 << boffset ) ;
			++boffset ;
		}
	}

	
	for ( j = 1 ; j <= lent ; ++j )
	{
			
	}
	
	delete[] Pv ;
	delete[] Mv ;
	delete[] Ph ;
	delete[] Mh ;
	delete[] computedBlockRange ;	
	for ( i = 0 ; i < 4 ; ++i )
		delete[] Peq[i] ;
}

// The text-book classical implementation, slow, just for debug purpose.
int GlobalAlignment_classic( char *t, int lent, char *p, int lenp, char *align )
{
	int i, j ;
	int *score = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
	int bmax = ( lent + 1 ) ;
	for ( i = 0 ; i <= lenp ; ++i )
		score[i * bmax + 0] = i ;
	for ( j = 0 ; j <= lent ; ++j )
		score[0 * bmax + j] = j ;

	for ( i = 1 ; i <= lenp ; ++i )
	{	
		for ( j = 1 ; j <= lent ; ++j )
		{
			int min = score[ ( i - 1 ) * bmax + j] + 1 ;
			if ( score[i * bmax + j - 1] + 1 < min )
				min = score[i * bmax + j - 1] + 1 ;
			int diag = score[ ( i - 1 ) * bmax + j - 1] + ( t[j - 1] == p[i - 1] ? 0 : 1 ) ;
			if ( diag < min )
				min = diag ;
			score[i * bmax + j] = min ;
		}
	}

	int tagi = lenp, tagj = lent ;
	int tag = 0 ;
	while ( tagi > 0 || tagj > 0 )
	{
		int min = score[tagi * bmax + tagj] ;
		int a = 0 ;
		if ( tagj > 0 && score[tagi * bmax + tagj - 1] + 1 == min )
			a = EDIT_DELETE ;	
		if ( tagi > 0 && score[ ( tagi - 1 ) * bmax + tagj] + 1 == min )
			a = EDIT_INSERT ; 
		if ( tagj > 0 && tagi > 0 )
		{
			int diff = ( t[tagj - 1] == p[tagi - 1] ? 0 : 1 ) ;
			if ( score[ ( tagi - 1) * bmax + tagj - 1 ] + diff == min )
			{
				if ( diff == 0 )
					a = EDIT_MATCH ;
				else
					a = EDIT_MISMATCH ;
			}
		}
		
		align[ tag ] = a ;
		++tag ;
		if ( a == EDIT_DELETE )
			--tagj ;
		else if ( a == EDIT_INSERT )
			--tagi ;
		else
		{
			--tagi ;
			--tagj ;
		}
	}
	align[tag] = -1 ;
	for ( i = 0, j = tag - 1 ; i < j ; ++i, --j )
	{	
		char tmp = align[i] ;
		align[i] = align[j] ;
		align[j] = tmp ;
	}
	//printf( "%d\n", score[lenp * bmax + lent] ) ;
	return score[lenp * bmax + lent] ;
}

void VisualizeAlignment( char *t, int lent, char *p, int lenp, char *align )
{
	int i, k, j;
	int tagt, tagp, taga ;
	taga = 0 ;
	int width = 100 ;
	k = j = 0 ;
	/*for ( i = 0 ; align[i] != -1 ; ++i )
		printf( "%d", align[i] ) ;
	printf( "\n" ) ;*/
	while ( align[taga] != -1 )
	{
		for ( i = taga ; i < taga + width && align[i] != -1 ; ++i  )
		{
			if ( align[i] == EDIT_INSERT )
				printf( "-" ) ;
			else
			{
				printf( "%c", t[k] ) ;
				++k ;
			}
		}
		printf( "\n" ) ;
		for ( i = taga ; i < taga + width && align[i] != -1 ; ++i )
		{
			if ( align[i] == EDIT_MATCH )
				printf( "|" ) ;
			else
				printf( " " ) ;
		}
		printf( "\n" ) ;
		for ( i = taga ; i < taga + width && align[i] != -1 ; ++i  )
		{
			if ( align[i] == EDIT_DELETE )
				printf( "-" ) ;
			else
			{
				printf( "%c", p[j] ) ;
				++j ;
			}
		}
		printf( "\n\n" ) ;
		taga = i ;
	}
}
