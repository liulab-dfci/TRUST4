#ifndef _MOURISL_POASET_HEADER
#define _MOURISL_POASET_HEADER

#include <string.h>
#include <algorithm>

#include "poa.hpp"
#include "SimpleVector.hpp"
#include "KmerIndex.hpp"
#include "ReadFiles.hpp"

enum POAType 
{
	NOVEL, VGENE, JGENE, CGENE
} ;

struct _poaWrapper
{
	POA poa ;
	char *name ;
} ;

struct _hit
{
	struct _indexInfo indexHit ;
	int readOffset ;
	int readStrand ;
	
	bool operator<( const struct _hit &b ) const
	{
		if ( indexHit.strand != b.indexHit.strand )
			return indexHit.strand < b.indexHit.strand ;
		else if ( indexHit.idx != b.indexHit.idx )
			return indexHit.idx < b.indexHit.idx ;
		else if ( readOffset != b.readOffset )
			return  readOffset < b.readOffset ;
		else if ( indexHit.offset != b.indexHit.offset )
			return indexHit.offset < b.indexHit.offset ;
		
		return false ;
	}
} ;

struct _overlap
{
	int poaIdx ;
	int readStart, readEnd ;
	int poaStart, poaEnd ;
	int poaStrand ;
} ;

class POASet
{
private:
	SimpleVector<struct _poaWrapper> poas ;
	KmerIndex poaIndex ;
	int kmerLength ;
	int minHitRequired ;

	static bool CompSortPairBInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.b < p2.b ;
	}
	

	// Return the first index whose hits.a is smaller or equal to valA
	int BinarySearch_LIS( int top[], int size, int valA, SimpleVector<struct _pair> &hits )
	{
		int l = 0, r = size - 1 ;
		int m ;
		while ( l <= r )
		{
			m = ( l + r ) / 2 ;
			if ( valA == hits[ top[m] ].a )
			{
				return m ;
			}
			else if ( valA < hits[ top[m] ].a )
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
	int LongestIncreasingSubsequence( SimpleVector<struct _pair> &hits, SimpleVector<struct _pair> &LIS ) 
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
					top[ret] = record[i] ;
					++ret ;
					link[ record[i] ] = top[tag] ;
				}
				else if ( hits[ record[i] ].a < hits[ top[tag + 1] ].a )
				{
					top[ tag + 1 ] = record[i] ;
					link[ record[i] ] = top[tag] ;
				}
			}
		}


		k = top[ret - 1] ;
		for ( i = ret - 1 ; i >= 0 ; --i )
		{
			LIS.PushBack( hits[k] ) ;
			k = link[k] ;	
		}
		LIS.Reverse() ;

		delete []top ;
		delete []record ;
		delete []link ;

		return ret ;
	}

public:
	POASet( int kl ) 
	{
		kmerLength = kl ;
		minHitRequired = 3 ;
	}
	~POASet() 
	{
		int size ;
		int i ;
		size = poas.Size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			//poas[i].poa.Release() ;
			free( poas[i].name ) ;	
		}
	}

	int Size()
	{
		return poas.Size() ;
	}
	
	// Input some baseline sequence to match against.
	void InputRefFa( char *filename ) 
	{
		ReadFiles fa ;
		fa.AddReadFile( filename ) ;
		
		KmerCode kmerCode( kmerLength ) ;
		while ( fa.Next() )
		{
			// Insert the kmers 
			struct _poaWrapper np ;
			np.name = strdup( fa.id ) ;

			int id = poas.Size() ;
			poas.PushBack( np ) ;

			struct _poaWrapper &pw = poas[id] ;
			int seqLen = strlen( fa.seq ) ;
			poas[id].poa.Initialize( fa.seq, seqLen, true  ) ;
			poaIndex.BuildIndexFromRead( kmerCode, fa.seq, seqLen, id ) ;
		}
	}

	// Compute the length of hit from the read, take the overlaps of kmer into account 
	int GetTotalHitLength( SimpleVector<struct _hit> &hits )
	{
		int hitSize = hits.Size() ;
		int i, j ;
		int ret = 0 ;
		//for ( i = 0 ; i < hitSize ; ++i )
		//	printf( "%d %d\n", i, hits[i].readOffset) ;
		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
			{
				if ( hits[j].readOffset > hits[j - 1].readOffset + kmerLength - 1 )	
					break ;
			}

			ret += hits[j - 1].readOffset - hits[i].readOffset + kmerLength ;

			i = j ;
		}
		return ret ;
	}

	// Use the hits to extract overlaps from PoaSet 
	int GetPoaOverlaps( SimpleVector<struct _hit> &hits, SimpleVector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int hitSize = hits.Size() ;
		
		SimpleVector<struct _pair> hitCoordDiff ;
		hitCoordDiff.Reserve( hitSize ) ;
		SimpleVector<struct _pair> concordantHitCoord ;
		SimpleVector<struct _pair> hitCoordLIS ;
		SimpleVector<struct _hit> finalHits ;

		int radius = 10 ;
		int hitLenRequired = 31  ;

		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].indexHit.strand != hits[i].indexHit.strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
					break ;
			//[i,j) holds the hits onto the same poa on the same strand.	
			if ( j - i < minHitRequired )
			{
				i = j ;
				continue ;
			}

			hitCoordDiff.Clear() ;
			for ( k = i ; k < j ; ++k )
			{
				struct _pair nh ;
				nh.a = k ;
				if ( hits[k].indexHit.strand == hits[k].readStrand )
					nh.b = hits[k].readOffset - hits[k].indexHit.offset ;
				else
					nh.b = hits[k].readOffset + hits[k].indexHit.offset ;
				hitCoordDiff.PushBack( nh ) ;
			}
			std::sort( hitCoordDiff.BeginAddress(), hitCoordDiff.EndAddress(), CompSortPairBInc ) ;

			// Pick the best concordant hits.
			int s, e ;
			struct _pair bestWindow ;
			int bestDiffSum = 10000 ;
			bestWindow.a = -1 ;
			bestWindow.b = -1 ;
			for ( s = i ; s < j ; )
			{
				int diffSum = 0 ;
				for ( e = s + 1 ; e < j ; ++e )
				{
					int diff = hitCoordDiff[e - i].b - hitCoordDiff[e - 1 - i].b ;
					if ( diff < 0 )
						diff = -diff ;

					if ( diff > radius ) 
						break ;
					diffSum += diff ; 
				}
			
				if ( e - s > bestWindow.b - bestWindow.a 
					|| ( e - s == bestWindow.b - bestWindow.a && diffSum < bestDiffSum ) )
				{
					bestWindow.a = s ;
					bestWindow.b = e ;
					bestDiffSum = diffSum ;
				}

				s = e ;
			}

			if ( bestWindow.b - bestWindow.a < minHitRequired
				|| ( bestWindow.b - bestWindow.a ) * kmerLength < hitLenRequired )
			{
				i = j ;
				continue ;
			}

			// [s, e) holds the candidate in the array of hitCoordDiff 
			// Note that hitCoordDiff is a subset from hits that lies on the same poa.
			s = bestWindow.a - i ;
			e = bestWindow.b - i ;

			concordantHitCoord.Clear() ;
			for ( k = s ; k < e ; ++k )
			{
				struct _pair nh ;
				int hitId = hitCoordDiff[k].a ; 
				nh.a = hits[ hitId ].readOffset ;
				nh.b = hits[ hitId ].indexHit.offset ;
				concordantHitCoord.PushBack( nh ) ;
			}

			if ( hits[i].indexHit.strand != hits[i].readStrand )
			{
				for ( k = 0 ; k < e - s ; ++k )
					concordantHitCoord[k].b = -concordantHitCoord[k].b ;
			}

			std::sort( concordantHitCoord.BeginAddress(), concordantHitCoord.EndAddress(), CompSortPairBInc ) ;
			//for ( k = 0 ; k < e - s ; ++k )	
			//	printf( "%d: %d %d %d\n", i, hits[i].indexHit.idx, concordantHitCoord[k].a, concordantHitCoord[k].b ) ;
			
			
			// Compute the longest increasing subsequence.
			hitCoordLIS.Clear() ;
			int lisSize = LongestIncreasingSubsequence( concordantHitCoord, hitCoordLIS ) ; 
			if ( lisSize * kmerLength < hitLenRequired )
			{
				i = j ;
				continue ;
			}

			// Rebuild the hits.
			finalHits.Clear() ;
			for ( k = 0 ; k < lisSize ; ++k )
			{
				struct _hit nh = hits[i];
				nh.readOffset = hitCoordLIS[k].a ;
				nh.indexHit.offset = hitCoordLIS[k].b ;
				if ( nh.readStrand != nh.indexHit.strand )
					nh.indexHit.offset = -nh.indexHit.offset ;
				//printf( "%d: %d %d %d %d\n", i, nh.readOffset, nh.indexHit.idx, nh.indexHit.offset, nh.indexHit.strand ) ;
				finalHits.PushBack( nh ) ;
			}

			int hitLen = GetTotalHitLength( finalHits ) ;
			if ( hitLen < hitLenRequired )
			{
				i = j ;
				continue ;
			}

			struct _overlap no ;
			no.poaIdx = hits[i].indexHit.idx ;
			no.readStart = finalHits[0].readOffset ;
			no.readEnd = finalHits[ lisSize - 1 ].readOffset + kmerLength - 1 ;
			no.poaStrand = finalHits[0].indexHit.strand ;
			if ( no.poaStrand == 1 )
			{
				no.poaStart = finalHits[0].indexHit.offset ;
				no.poaEnd = finalHits[ lisSize - 1 ].indexHit.offset + kmerLength - 1 ;
			}
			else
			{
				// if on other strand, the index starting the kmer offset from the other end.
				no.poaStart = finalHits[0].indexHit.offset - kmerLength + 1 ;
				no.poaEnd = finalHits[ lisSize - 1].indexHit.offset ;
			}
			overlaps.PushBack( no ) ;
			i = j ;
		}
		return overlaps.Size() ;
	}

	
	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set.
	int AddRead( char *seq )
	{
		//printf( "%s\n", seq ) ;
		int i, j, k ;
		int len = strlen( seq ) ;
		if ( len < kmerLength )
			return -1 ;

		SimpleVector<struct _hit> hits ;		    	

		KmerCode kmerCode( kmerLength ) ;

		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( seq[i] ) ;

		for ( ; i < len ; ++i )
		{
			kmerCode.Append( seq[i] ) ;
			SimpleVector<struct _indexInfo> &indexHit = *poaIndex.Search( kmerCode ) ; 
			
			int size = indexHit.Size() ;

			//if ( size >= 40 )
			//	continue ;
			
			for ( j = 0 ; j < size ; ++j )
			{
				struct _hit nh ;
				nh.indexHit = indexHit[j] ;
				nh.readOffset = i - kmerLength + 1 ;
				nh.readStrand = 1 ;
				hits.PushBack( nh ) ;
			}
		}
		
		std::sort( hits.BeginAddress(), hits.EndAddress() ) ;
		//for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
		//	printf( "- %d %d %d %d\n", it->readOffset, it->indexHit.idx, it->indexHit.offset, it->indexHit.strand ) ;
		
		SimpleVector<struct _overlap> overlaps ;
		int overlapCnt = GetPoaOverlaps( hits,  overlaps ) ;
		
		// Determine whether we want to add this reads.
		for ( i = 0 ; i < overlapCnt ; ++i )
			printf( "%d: %d %s. %d %d %d %d\n", i, overlaps[i].poaIdx, poas[ overlaps[i].poaIdx ].name, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].poaStart, overlaps[i].poaEnd ) ;

		return 0 ;
	}
} ;


#endif
