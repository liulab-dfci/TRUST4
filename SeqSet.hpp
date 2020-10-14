// The data structure holds the set of sequences (can be "assembled" from several reads)
#ifndef _MOURISL_SEQSET_HEADER
#define _MOURISL_SEQSET_HEADER

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <string>

#include "SimpleVector.hpp"
#include "KmerIndex.hpp"
#include "ReadFiles.hpp"
#include "AlignAlgo.hpp"


struct _seqWrapper
{
	char *name ;
	char *consensus ; // This should be handled by malloc/free.
	int consensusLen ;
	SimpleVector<struct _posWeight> posWeight ;
	bool isRef ; // this is from reference.

	int minLeftExtAnchor, minRightExtAnchor ; // only overlap with size larger than this can be counted as valid extension.

	struct _triple info[3] ; // For storing extra information. for ref, info[0,1] contains the coordinate for CDR1,2 and info[2].a for CDR3
					// In extending seqs with mate pair information, these are used to store rough V, J, C read coordinate.	
	int barcode ; // transformed barcode. -1: no barcode 
} ;

struct _hit
{
	struct _indexInfo indexHit ;
	int readOffset ;
	int strand ; // -1: different strand, 1: same strand. When strand==-1, the readOffset is the offset in the rcSeq.
	
	int repeats ; // record how many times this hit with other index part.

	bool operator<( const struct _hit &b ) const
	{
		if ( strand != b.strand )
			return strand < b.strand ;
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
	int seqIdx ;
	int readStart, readEnd ; // When strand ==-1, the start,end is the offset in the rcSeq.
	int seqStart, seqEnd ;
	int strand ;
	
	int matchCnt ; // The number of matched bases, count TWICE.
	double similarity ;

	SimpleVector<struct _pair> *hitCoords ;
	SimpleVector<int> *info ; // store extra informations 

	bool operator<( const struct _overlap &b ) const
	{
		// The overlap with more matched bases should come first.
		if ( matchCnt > b.matchCnt + 2 || matchCnt < b.matchCnt - 2 )
			return matchCnt > b.matchCnt ;
		else if ( similarity != b.similarity )
			return similarity > b.similarity ; 
		else if ( readEnd - readStart != b.readEnd - b.readStart )
			return readEnd - readStart > b.readEnd - b.readStart ;
		else if ( strand !=  b.strand )
			return strand < b.strand ;
		else if ( seqStart != b.seqStart )
			return seqStart < b.seqStart ;
		else if ( seqEnd != b.seqEnd )
			return seqEnd < b.seqEnd ; 
		else if ( readStart != b.readStart )
			return readStart < b.readStart ;
		else if ( readEnd != b.readEnd )
			return readEnd < b.readEnd ;
		else 
			return seqIdx < b.seqIdx ;

		return false ;
	}
} ;

struct _assignRead
{
	char *id ;
	char *read ;
	int barcode ;
	int umi ;

	int info ; 
	struct _overlap overlap ;
} ;

struct _assignHistory
{
	int curIdx ;
	int curMateIdx ;

	SimpleVector<int> readsAssignedTo ;
	std::vector< std::vector<int> >  readsInSeq ;
} ;


class SeqSet
{
private:
	std::vector<struct _seqWrapper> seqs ;
	KmerIndex seqIndex ;
	int kmerLength ;
	int radius ;
	int hitLenRequired ;
	int gapN ;
	bool isLongSeqSet ; // Whether this seq set is built from long reads. Long reads may require more drastic filtration.

	// Some threshold
	double novelSeqSimilarity ;
	double refSeqSimilarity ;
	double repeatSimilarity ; // e.g., the repeat when building the branch graph.

	struct _overlap prevAddInfo ; 

	static bool CompSortPairBInc( const struct _pair &p1, const struct _pair &p2 )
	{
		if ( p1.b != p2.b )
			return p1.b < p2.b ;
		else
			return p1.a < p2.a ;
	}
	
	static bool CompSortPairAInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.a < p2.a ;
	}

	static bool CompSortOverlapsOnReadCoord( const struct _overlap &a, const struct _overlap &b )
	{
		return a.readStart < b.readStart ; 
	}

	static bool CompSortAssignedReadById( const struct _assignRead &a, const struct _assignRead &b )
	{
		return strcmp( a.id, b.id ) < 0 ;
	}
		
	static bool CompSortOverlapByCoord( const struct _overlap &a, const struct _overlap &b )	
	{
		if ( a.seqIdx != b.seqIdx )
			return a.seqIdx < b.seqIdx ;
		else if ( a.readStart != b.readStart )
			return a.readStart < b.readStart ;
		else
			return a.readEnd < b.readEnd ;
	}

	static bool CompSortHitCoordDiff( const struct _triple &a, const struct _triple &b )
	{
		if ( a.c != b.c )
			return a.c < b.c ;
		else if ( a.b != b.b )
			return a.b < b.b ;
		else
			return a.a < b.a ;
	}

	bool IsReverseComplement( char *a, char *b )
	{
		int i, j ;
		int len = strlen( a ) ;
		if ( len != strlen( b) ) 
			return false ;
		for ( i = 0, j = len - 1 ; i < len ; ++i, --j )
			if ( a[i] == 'N' && b[j] == 'N' )
				continue ;
			else if ( a[i] != 'N' && b[j] != 'N' )
			{
				if ( 3 - nucToNum[ a[i] - 'A' ] != nucToNum[ b[j] - 'A' ] )
					return false ;
			}
			else
				return false ;
		return true ;
	}

	void Reverse( char *r, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
			r[i] = seq[len - 1 - i] ;  
		r[i] = '\0' ;
	}
	
	bool IsPosWeightCompatible( const struct _posWeight &a, const struct _posWeight &b )
	{
		int sumA = a.Sum() ;
		int sumB = b.Sum() ;
	
		if ( sumA == 0 || sumB == 0 
			|| ( sumA < 3 * a.count[0] && sumB < 3 * b.count[0] ) 
			|| ( sumA < 3 * a.count[1] && sumB < 3 * b.count[1] )
			|| ( sumA < 3 * a.count[2] && sumB < 3 * b.count[2] )
			|| ( sumA < 3 * a.count[3] && sumB < 3 * b.count[3] ) )
			return true ;
		return false ;	
	}

	bool IsOverlapIntersect( const struct _overlap &a, const struct _overlap &b )
	{
		if ( a.seqIdx == b.seqIdx && 
			( ( a.seqStart <= b.seqStart && a.seqEnd >= b.seqStart ) 
			|| ( b.seqStart <= a.seqStart && b.seqEnd >= a.seqStart ) ) )
			return true ;
		return false ;
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
			//if ( hits[i].b == hits[i - 1].b )
			//	continue ;
			record[rcnt] = i ;
			++rcnt ;
		}
		top[0] = 0 ;
		link[0] = -1 ;
		ret = 1 ;
		for ( i = 1 ; i < rcnt ; ++i )
		{
			int tag = 0 ;
			if ( hits[ top[ ret - 1 ] ].a <= hits[ record[i] ].a )
				tag = ret - 1 ;
			else
				tag = BinarySearch_LIS( top, ret, hits[ record[i] ].a, hits ) ;			
			
			if ( tag == -1 )
			{
				top[0] = record[i] ;
				link[ record[i] ] = -1 ;
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
		//for ( i = 0 ; i < ret ; ++i )
		//	LIS.PushBack( hits[ top[i] ] ) ;
		
		// Remove elements with same b.
		if ( ret > 0 )
		{
			k = 1 ;
			for ( i = 1 ; i < ret ; ++i )
			{
				if ( LIS[i].b == LIS[k - 1].b )
					continue ;
				LIS[k] = LIS[i] ;
				++k ;
			}
			ret = k ;
		}

		delete []top ;
		delete []record ;
		delete []link ;

		return ret ;
	}

	void GetAlignStats( char *align, bool update, int &matchCnt, int &mismatchCnt, int &indelCnt)
	{
		int k ;
		if ( !update )
		{
			matchCnt = mismatchCnt = indelCnt = 0 ;
		}

		for ( k = 0 ; align[k] != -1 ; ++k )
		{
			if ( align[k] == EDIT_MATCH )
				++matchCnt ;
			else if ( align[k] == EDIT_MISMATCH )
				++mismatchCnt ;
			else 
				++indelCnt ;
		}
	}


	bool IsOverlapLowComplex( char *r, struct _overlap &o )
	{
		int cnt[4] = {0, 0, 0, 0} ;
		int i ;
		for ( i = o.readStart ; i <= o.readEnd ; ++i )
		{
			if ( r[i] == 'N' )
				continue ;
			++cnt[ nucToNum[ r[i] - 'A' ] ] ;
		}
		int len = o.readEnd - o.readStart + 1 ;
		int lowCnt = 0 ; 
		int lowTotalCnt = 0 ;
		for ( i = 0 ; i < 4 ; ++i )
		{
			if ( cnt[i] <= 2 )
			{
				++lowCnt ;
				lowTotalCnt += cnt[i] ;
			}
		}
		if ( lowTotalCnt * 7 >= o.readEnd - o.readStart + 1 )
			return false ;

		if ( lowCnt >= 2 )
			return true ;
		return false ;
	}

	bool IsEquivalentConstantGene( char *a, char *b )
	{
		int i ;
		for ( i = 0 ; i < 4 ; ++i )
			if ( a[i] != b[i] )
				return false ;
		return true ;
	}

	void SetPrevAddInfo( int seqIdx, int readStart, int readEnd, int seqStart, int seqEnd, int strand )
	{
		prevAddInfo.seqIdx = seqIdx ;
		prevAddInfo.readStart = readStart ;
		prevAddInfo.readEnd = readEnd ;
		prevAddInfo.seqStart = seqStart ;
		prevAddInfo.seqEnd = seqEnd ;
		prevAddInfo.strand = strand ;
	}
	
	char DnaToAa( char a, char b, char c )
	{
		if ( a == 'N' || b == 'N' || c == 'N' )
			return '-' ;
		if ( a == 'M' || b == 'M' || c == 'M' )
			return '-' ;

		if ( a == 'A' )
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return 'K' ;
				else
					return 'N' ;
			}
			else if ( b == 'C' )
			{
				return 'T' ;
			}
			else if ( b == 'G' )
			{
				if ( c == 'A' || c == 'G' )
					return 'R' ;
				else
					return 'S' ;
			}
			else
			{
				if ( c == 'G' )
					return 'M' ;
				else 
					return 'I' ;
			}
		}
		else if ( a == 'C' )
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return 'Q' ;
				else
					return 'H' ;
				
			}
			else if ( b == 'C' )
			{
				return 'P' ;
			}
			else if ( b == 'G' )
			{
				return 'R' ;
			}
			else
			{
				return 'L' ;
			}
		}
		else if ( a == 'G' )
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return 'E' ;
				else
					return 'D' ;
			}
			else if ( b == 'C' )
			{
				return 'A' ;
			}
			else if ( b == 'G' )
			{
				return 'G' ;
			}
			else
			{
				return 'V' ;
			}
		}
		else
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return '*' ;
				else
					return 'Y' ;
			}
			else if ( b == 'C' )
			{
				return 'S' ;
			}
			else if ( b == 'G' )
			{
				if ( c == 'A' )
					return '*' ;
				else if ( c == 'G' )
					return 'W' ;
				else
					return 'C' ;
				
			}
			else
			{
				if ( c == 'A' || c == 'G' )
					return 'L' ;
				else
					return 'F' ;
			}
		}
	}

	void ReleaseSeq( int idx )
	{
		if ( seqs[idx].consensus == NULL )
			return ;
		free( seqs[idx].name ) ;
		free( seqs[idx].consensus ) ;
		seqs[idx].posWeight.Release() ;

		seqs[idx].name = seqs[idx].consensus = NULL ;
	}
	
	// Use the hits to extract overlaps from SeqSet
	int GetOverlapsFromHits( SimpleVector<struct _hit> &hits, int hitLenRequired, int filter, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int hitSize = hits.Size() ;
		
		SimpleVector<struct _triple> hitCoordDiff ;
		hitCoordDiff.Reserve( hitSize ) ;
		SimpleVector<struct _pair> concordantHitCoord ;
		SimpleVector<struct _pair> hitCoordLIS ;
		SimpleVector<struct _hit> finalHits ;
		
		// Compute the minHitRequired. 
		// NOTE: each strand should have its own minHitRequired, it could be that on one strand,
		//    each hit is matched to too many places and the skip hits mechanism is triggered.
		int novelMinHitRequired[2] = {3, 3} ;
		int refMinHitRequired[2] = {3, 3} ;
		bool removeOnlyRepeats[2] = {false, false} ; // Remove the hits on a seq that are all repeats hit.
		int possibleOverlapCnt[2] = {0, 0} ;
		if ( filter == 1 )
		{
			int longestHits[2] = {0, 0} ;
			for ( i = 0 ; i < hitSize ; ++i )
			{
				int isPlusStrand = ( 1 + hits[i].strand ) / 2 ; 
				for ( j = i + 1 ; j < hitSize ; ++j )
					if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
						break ;
				if ( !seqs[ hits[i].indexHit.idx].isRef )
				{
					if ( j - i > novelMinHitRequired[ isPlusStrand ] )
						++possibleOverlapCnt[ isPlusStrand ] ;
					if ( j - i > longestHits[ isPlusStrand ] )
						longestHits[ isPlusStrand] = j - i  ;
				}
				
				if ( !removeOnlyRepeats[isPlusStrand] )
				{
					int cnt = 0 ;
					for ( k = i ; k < j ; ++k )
						if ( hits[k].repeats <= 10000 )
							++cnt ;
					if ( cnt >= novelMinHitRequired[ isPlusStrand ] )
					{
						removeOnlyRepeats[ isPlusStrand ] = true ;
					}
				}

				i = j ;
			}
			// filter based on the repeatability of overlaps.
			for ( i = 0 ; i <= 1 ; ++i )
			{
				if ( possibleOverlapCnt[i] > 100000 )
					novelMinHitRequired[i] = longestHits[i] * 0.75 ;
				else if ( possibleOverlapCnt[i] > 10000 )
					novelMinHitRequired[i] = longestHits[i] / 2 ;
				else if ( possibleOverlapCnt[i] > 1000 )
					novelMinHitRequired[i] = longestHits[i] / 3 ;
				else if ( possibleOverlapCnt[i] > 100 )
					novelMinHitRequired[i] = longestHits[i] / 4 ;
			}
		}

		//if ( novelMinHitRequired > 3 )
		//	printf( "novelMinHitRequired=%d\n", novelMinHitRequired ) ;
		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
					break ;

			int minHitRequired = novelMinHitRequired[ ( 1 + hits[i].strand ) / 2 ] ;
			if ( seqs[ hits[i].indexHit.idx].isRef )
				minHitRequired = refMinHitRequired[ ( 1 + hits[i].strand ) / 2 ];
			

			/*if ( filter == 1 && readLen > 0 )
			{
				int readStart = hits[i].readOffset, readEnd = hits[j - 1].readOffset + kmerLength - 1 ;
				int seqStart = seqs[ hits[j - 1].indexHit.idx ].consensusLen, seqEnd = -1 ;

				for ( k = i ; k < j ; ++k )
				{
					if ( hits[k].indexHit.offset < seqStart )
						seqStart = hits[k].indexHit.offset ;
					if ( hits[k].indexHit.offset > seqEnd )
						seqEnd = hits[k].indexHit.offset ;
				}
				seqEnd += kmerLength - 1 ;

				int leftOverhangSize = MIN( readStart, seqStart ) ;
				int rightOverhangSize = MIN( readLen - 1 - readEnd, 
						seqs[ hits[i].indexHit.idx ].consensusLen - 1 - seqEnd ) ;
			
				if ( leftOverhangSize > 2 * hitLenRequired || rightOverhangSize > 2 * hitLenRequired )
				{
					i = j ;
					continue ;
				}
			}*/

			//[i,j) holds the hits onto the same seq on the same strand.	
			if ( j - i < minHitRequired )
			{
				i = j ;
				continue ;
			}
			
			if ( removeOnlyRepeats[( 1 + hits[i].strand ) / 2] )
			{
				bool hasUnique = false ;
				for ( k = i ; k < j ; ++k )
				{
					if ( hits[k].repeats <= 10000 )
					{
						hasUnique = true ;
						break ;
					}
				}
				if ( !hasUnique )
				{
					i = j ;
					continue ;
				}
			}

			hitCoordDiff.Clear() ;
			for ( k = i ; k < j ; ++k )
			{
				struct _triple nh ;
				nh.a = hits[k].readOffset ;
				nh.b = hits[k].indexHit.offset ;
				nh.c = hits[k].readOffset - hits[k].indexHit.offset ;
				hitCoordDiff.PushBack( nh ) ;
			}
			std::sort( hitCoordDiff.BeginAddress(), hitCoordDiff.EndAddress(), CompSortHitCoordDiff ) ;

			// Pick the best concordant hits.
			int s, e ;
			int adjustRadius = radius ;
			if ( !seqs[ hits[i].indexHit.idx ].isRef )
				adjustRadius = 0 ;

			for ( s = 0 ; s < j - i ; )
			{
				int diffSum = 0 ;
				for ( e = s + 1 ; e < j - i ; ++e )
				{
					int diff = hitCoordDiff[e].c - hitCoordDiff[e - 1].c ;
					if ( diff < 0 )
						diff = -diff ;
					
					if ( diff > adjustRadius ) 
						break ;

					diffSum += diff ; 
				}
				//printf( "%d %d: %d %d\n", i, j, s, e ) ;


				if ( e - s < minHitRequired 
					|| ( e - s ) * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}


				if ( removeOnlyRepeats[( 1 + hits[i].strand ) / 2] )
				{
					bool hasUnique = false ;
					for ( k = s ; k < e ; ++k )
					{
						if ( hits[k].repeats <= 10000 )
						{
							hasUnique = true ;
							break ;
						}
					}
					if ( !hasUnique )
					{
						s = e ;
						continue ;
					}
				}
				// [s, e) holds the candidate in the array of hitCoordDiff 
				concordantHitCoord.Clear() ;
				for ( k = s ; k < e ; ++k )
				{
					struct _pair nh ;
					nh.a = hitCoordDiff[k].a ;
					nh.b = hitCoordDiff[k].b ;
					concordantHitCoord.PushBack( nh ) ;
				}
				if ( adjustRadius > 0 )
					std::sort( concordantHitCoord.BeginAddress(), concordantHitCoord.EndAddress(), CompSortPairBInc ) ;
				//for ( k = 0 ; k < e - s ; ++k )	
				//	printf( "%d (%d-%d): %d %s %d %d\n", i, s, e, hits[i].indexHit.idx, seqs[ hits[i].indexHit.idx ].name, concordantHitCoord[k].a, concordantHitCoord[k].b ) ;

				// Compute the longest increasing subsequence.
				//printf( "lis for %d (%d %d; %d %d). strand=%d (%d)\n", e - s, i, j, s, e, hits[i].strand, seqs.size() ) ;
				hitCoordLIS.Clear() ;
				int lisSize = LongestIncreasingSubsequence( concordantHitCoord, hitCoordLIS ) ; 
				if ( lisSize * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}

				// Rebuild the hits.
				int lisStart = 0 ;
				int lisEnd = lisSize - 1 ;
				// Ignore long insert gaps.
				if ( isLongSeqSet )
				{
					int maxGap = 2 * hitLenRequired + 3 * kmerLength ;
					if ( filter == 0 )//&& possibleOverlapCnt[( 1 + hits[i].strand ) / 2] > 1000 )
						maxGap *= 4 ;
					if ( maxGap < 200 )
						maxGap = 200 ;
					int max = -1 ;
					for ( k = 0 ; k < lisSize ; )
					{
						int l ;
						for ( l = k + 1 ; l < lisSize ; ++l )
						{
							if ( hitCoordLIS[l].a - hitCoordLIS[l - 1].a > maxGap )
								break ;
						}
						if ( l - k > max )
						{
							max = l - k ;
							lisStart = k ;
							lisEnd = l - 1 ;
						}

						k = l ;	
					}
				}

				finalHits.Clear() ;
				for ( k = lisStart ; k <= lisEnd ; ++k )
				{
					struct _hit nh = hits[i];
					nh.readOffset = hitCoordLIS[k].a ;
					nh.indexHit.offset = hitCoordLIS[k].b ;
					//if (seqs.size() == 1 )
					//	printf( "%d: %d %d %d %d\n", i, nh.readOffset, nh.indexHit.idx, nh.indexHit.offset, nh.strand ) ;
					finalHits.PushBack( nh ) ;
				}
				lisSize = lisEnd - lisStart + 1 ;

				int hitLen = GetTotalHitLengthOnRead ( finalHits ) ;
				if ( hitLen < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				else if ( GetTotalHitLengthOnSeq( finalHits ) < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				
				struct _overlap no ;
				no.seqIdx = hits[i].indexHit.idx ;
				no.readStart = finalHits[0].readOffset ;
				no.readEnd = finalHits[ lisSize - 1 ].readOffset + kmerLength - 1 ;
				no.strand = finalHits[0].strand ;
				no.seqStart = finalHits[0].indexHit.offset ;
				no.seqEnd = finalHits[ lisSize - 1 ].indexHit.offset + kmerLength - 1 ;
				no.matchCnt = 2 * hitLen ;
				no.similarity = 0 ;

				if ( !seqs[ no.seqIdx ].isRef && hitLen * 2 < no.seqEnd - no.seqStart + 1 )
				{
					s = e ; 
					continue ;
				}
				no.hitCoords = new SimpleVector<struct _pair> ;
				no.hitCoords->Reserve( lisSize ) ;
				for ( k = 0 ; k < lisSize ; ++k )
				{
					struct _pair nh ;
					nh.a = finalHits[k].readOffset ;
					nh.b = finalHits[k].indexHit.offset ;
					no.hitCoords->PushBack( nh ) ;
				}
				overlaps.push_back( no ) ;

				s = e ;
			} // iterate through concordant hits.
			i = j ;
		}
		return overlaps.size() ;
	}
	
	// Find the overlaps from hits if it possibly span the CDR3 region and anchor paritally on V and J gene 
	int GetVJOverlapsFromHits( SimpleVector<struct _hit> &hits, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		SimpleVector<struct _hit> VJhits ; 		
		
		int hitSize = hits.Size() ;
		
		// Filter hits that are out of VJ junction region.
		VJhits.Reserve( hitSize ) ;
		for ( i = 0 ; i < hitSize ; ++i )
		{
			int seqIdx = hits[i].indexHit.idx ;
			if ( !seqs[ seqIdx ].isRef )
				continue ;

			if ( seqs[ seqIdx ].name[3] == 'V' && hits[i].indexHit.offset >= seqs[ seqIdx ].consensusLen - 31 )
			{
				VJhits.PushBack( hits[i] ) ;
			}
			else if ( seqs[ seqIdx ].name[3] == 'J' && hits[i].indexHit.offset < 31 )
			{
				VJhits.PushBack( hits[i] ) ;
			}
		}

		GetOverlapsFromHits( VJhits, 17, 0, overlaps ) ;
		
		// Extract the best VJ pair 
		int overlapCnt = overlaps.size() ;
		int maxMatchCnt = 0 ;
		int tagi = 0, tagj = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			for ( j = i + 1 ; j < overlapCnt ; ++j )
			{
				int seqIdxI = overlaps[i].seqIdx ;
				int seqIdxJ = overlaps[j].seqIdx ;

				if ( seqs[ seqIdxI ].name[0] != seqs[ seqIdxJ ].name[0] ||
					seqs[ seqIdxI ].name[1] != seqs[ seqIdxJ ].name[1] ||
					seqs[ seqIdxI ].name[2] != seqs[ seqIdxJ ].name[2] ||
					seqs[ seqIdxI ].name[3] == seqs[ seqIdxJ ].name[3] )
					continue ;			

				if ( seqs[ seqIdxI ].name[3] == 'V' )
				{
					if ( overlaps[i].readStart > overlaps[j].readStart )
						continue ;
				}
				else 
				{
					if ( overlaps[i].readStart < overlaps[j].readStart )
						continue ;
				}

				if ( overlaps[i].matchCnt + overlaps[j].matchCnt > maxMatchCnt )
				{
					maxMatchCnt = overlaps[i].matchCnt + overlaps[j].matchCnt ;
					tagi = i ;
					tagj = j ;
				}
			}
		}

		if ( maxMatchCnt == 0 )
		{
			int size = overlaps.size() ;
			for ( i = 0 ; i < size ; ++i )
			{
				overlaps[i].hitCoords->Release() ;
				delete overlaps[i].hitCoords ;
				overlaps[i].hitCoords = NULL ;
			}
			overlaps.clear() ;
			return 0 ;
		}
		int size = overlaps.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( i == tagi || i == tagj )
				continue ;

			overlaps[i].hitCoords->Release() ;
			delete overlaps[i].hitCoords ;
			overlaps[i].hitCoords = NULL ;
		}


		std::vector<struct _overlap> ret ;
		ret.push_back( overlaps[ tagi ] ) ;
		ret.push_back( overlaps[ tagj ] ) ;
		
		overlaps = ret ;

		return 2 ;
	}

	// Extend the overlap to include the overhang parts and filter the overlaps if the overhang does not match well.
	// return: whether this is a valid extension or not
	int ExtendOverlap( char *r, int len, struct _seqWrapper &seq, double mismatchThresholdFactor, 
		char *align, struct _overlap &overlap, struct _overlap &extendedOverlap )
	{
		// Check whether the overhang part is compatible with each other or not.
		// Extension to 5'-end ( left end )
		int matchCnt, mismatchCnt, indelCnt ;
		int leftOverhangSize = MIN( overlap.readStart, overlap.seqStart ) ;
		int ret = 1 ;
		int i, k ;
		int goodLeftOverhangSize = 0 ;
		//AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqStart - leftOverhangSize,
		AlignAlgo::GlobalAlignment_PosWeight( seq.posWeight.BeginAddress() + overlap.seqStart - leftOverhangSize, 
				leftOverhangSize, 
				r + overlap.readStart - leftOverhangSize, leftOverhangSize, align ) ;
		GetAlignStats( align, false, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > 0 )
		{
			leftOverhangSize = 0 ;
			ret = 0 ;
		}
		for ( i = 0 ; align[i] != -1 ; ++i )
			;
		int tmpMatchCnt = 0  ;
		for ( i = i - 1, k = 1 ; i >= 0 ; --i, ++k )
		{
			if ( align[i] == EDIT_MATCH )
			{
				++tmpMatchCnt ;
				if ( tmpMatchCnt > 0.75 * k )
					goodLeftOverhangSize = k ;
			}
			else if ( align[i] != EDIT_MISMATCH )
				break ;
		}
		// Extension to 3'-end ( right end )
		int rightOverhangSize = MIN( len - 1 - overlap.readEnd, seq.consensusLen - 1 - overlap.seqEnd ) ;
		int goodRightOverhangSize = 0 ;
		//AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqEnd + 1, 
		AlignAlgo::GlobalAlignment_PosWeight( seq.posWeight.BeginAddress() + overlap.seqEnd + 1, 
				rightOverhangSize,
				r + overlap.readEnd + 1, rightOverhangSize, align ) ;
		int oldIndelCnt = indelCnt ;
		GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > oldIndelCnt )
		{
			rightOverhangSize = 0 ;
			ret = 0 ;
		}
		tmpMatchCnt = 0 ;
		for ( i = 0 ; align[i] != -1 ; ++i )
		{
			if ( align[i] == EDIT_MATCH )
			{
				++tmpMatchCnt ;
				if ( tmpMatchCnt > 0.75 * ( i + 1 ) )
					goodRightOverhangSize = i + 1 ;
			}
			else if ( align[i] != EDIT_MISMATCH )
				break ;
		}

		int mismatchThreshold = 2 ;
		if ( leftOverhangSize >= 2 )
			++mismatchThreshold ;
		if ( rightOverhangSize >= 2 )
			++mismatchThreshold ;

		double densityThreshold = 1.5 / kmerLength ;
		//if ( leftOverhangSize + rightOverhangSize > 2 * kmerLength )
		//	densityThreshold = 1.0 / 6 ;
		mismatchThreshold *= mismatchThresholdFactor ;
		if ( mismatchCnt > mismatchThreshold && (double)mismatchCnt / ( leftOverhangSize + rightOverhangSize ) > densityThreshold ) 
			ret = 0 ;
		/*if ( !strcmp(r, "AACCTTCACCTACACGCCCTGCAGCCAGAAGACTCAGCCCTGTATCTCTGCGCCAGCAGCGGGACATTTGGTTCACCCCTCCACTCTGGGAACGGGACCCGGCTCTCAGTG"))
		{
			fprintf(stderr, "%s: %d %d %d %d. %d\n", seq.name, mismatchCnt, mismatchThreshold, leftOverhangSize, rightOverhangSize,
				kmerLength );
		}*/
		extendedOverlap.seqIdx = overlap.seqIdx ;
		extendedOverlap.readStart = overlap.readStart - leftOverhangSize ;
		extendedOverlap.readEnd = overlap.readEnd + rightOverhangSize ;
		extendedOverlap.seqStart = overlap.seqStart - leftOverhangSize ;
		extendedOverlap.seqEnd = overlap.seqEnd + rightOverhangSize ;
		extendedOverlap.strand = overlap.strand ;	
		extendedOverlap.matchCnt = 2 * matchCnt + overlap.matchCnt ;
		extendedOverlap.similarity = (double)( 2 * matchCnt + overlap.matchCnt ) / 
			( extendedOverlap.readEnd - extendedOverlap.readStart + 1 + extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ) ;	
			
		if ( ( seqs[ extendedOverlap.seqIdx ].isRef && extendedOverlap.similarity < refSeqSimilarity ) 
			|| ( !seqs[ extendedOverlap.seqIdx ].isRef && extendedOverlap.similarity < novelSeqSimilarity ) )
		{
			extendedOverlap = overlap ;
			ret = 0 ;
		}

		if ( ret == 0 )
		{
			extendedOverlap.readStart = overlap.readStart - goodLeftOverhangSize ;
			extendedOverlap.readEnd = overlap.readEnd + goodRightOverhangSize ;
			extendedOverlap.seqStart = overlap.seqStart - goodLeftOverhangSize ;
			extendedOverlap.seqEnd = overlap.seqEnd + goodRightOverhangSize ;
		}

		/*if ( !seqs[ extendedOverlap.seqIdx ].isRef && 
			extendedOverlap.readEnd - extendedOverlap.readStart + 1 - extendedOverlap.matchCnt / 2 >= 5 )
		{
			// Exceed the number of mismatch allowed
			extendedOverlap = overlap ;
			ret = 0 ;
		}*/
		//printf( "%d: %d %d\n", ret, extendedOverlap.readStart, extendedOverlap.readEnd ) ;
		return ret ;
	}
	
	// Test at sequences level, whether overlap a is a substring of b. strict controls whether allow them to equal 
	bool IsOverlapSubstringOf( struct _overlap a, struct _overlap b, bool strict, int maxMismatch )
	{
		int seqIdxA = a.seqIdx ;
		int seqIdxB = b.seqIdx ;
		int i, j ;
		if ( seqIdxA == -1 || seqIdxB == -1 )
			return false ;

		if ( a.readStart < b.readStart || a.readEnd > b.readEnd )
			return false ;

		if ( strict && a.readEnd - a.readStart == b.readEnd - b.readStart )
			return false ;

		int offset = a.readStart - b.readStart ;
		int mismatchCnt = 0 ;
		for ( i = a.seqStart, j = b.seqStart + offset ; i <= a.seqEnd ; ++i, ++j )
		{
			if ( seqs[seqIdxA].consensus[i] != seqs[seqIdxB].consensus[j] )
				++mismatchCnt ;
			if ( mismatchCnt > maxMismatch )
				return false ;
		}
		return true ;
	}

	void SortHits( SimpleVector<struct _hit> &hits, bool alreadyReadOrder )
	{
		int i, k ;
		if ( hits.Size() > 2 * seqs.size() && alreadyReadOrder ) 
		{
			// Bucket sort.
			int hitCnt = hits.Size() ;
			int seqCnt = seqs.size() ;
			SimpleVector<struct _hit> *buckets[2] ;
			buckets[0] = new SimpleVector<struct _hit>[seqCnt] ;
			buckets[1] = new SimpleVector<struct _hit>[seqCnt] ;

			for ( i = 0 ; i < hitCnt ; ++i )
			{
				int tag = hits[i].strand == 1 ? 1 : 0 ;
				buckets[tag][ hits[i].indexHit.idx ].PushBack( hits[i] ) ;
			}
			
			hits.Clear() ;
			for ( k = 0 ; k <= 1 ; ++k )
			{
				for ( i = 0 ; i < seqCnt ; ++i )
				{
					hits.PushBack( buckets[k][i] ) ;
				}
			}

			delete[] buckets[0] ;
			delete[] buckets[1] ;
		}
		else
			std::sort( hits.BeginAddress(), hits.EndAddress() ) ;

	}

	int GetHitsFromRead( char *read, char *rcRead, int strand, int barcode, bool allowTotalSkip, SimpleVector<struct _hit> &hits, SimpleVector<bool> *puse )
	{
		int i, j ;
		int len = strlen( read ) ;

		KmerCode kmerCode( kmerLength ) ;
		KmerCode prevKmerCode( kmerLength ) ;

		// Locate the hits from the same-strand case.
		//int skipLimit = 3 ;
		int skipLimit = kmerLength / 2 ; 

		int skipCnt = 0 ;
		int downSample = 1 ;
		if ( len > 200 && isLongSeqSet )
		{
			downSample = 1 + len / 200 ;
			//skipLimit /= downSample ;
			//if ( skipLimit < 2 )
			//	skipLimit = 2 ; 
		}

		if ( strand != -1 )
		{
			for ( i = 0 ; i < kmerLength - 1 ; ++i )
				kmerCode.Append( read[i] ) ;

			for ( ; i < len ; ++i )
			{
				kmerCode.Append( read[i] ) ;
				if ( downSample > 1 && ( i - kmerLength + 1 ) % downSample != 0 )
					continue ;

				if ( i == kmerLength - 1 || !prevKmerCode.IsEqual( kmerCode ) )
				{
					SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

					int size = indexHit.Size() ;
					if ( size >= 100 && puse == NULL /*&& barcode == -1*/ && i != kmerLength - 1 && i != len - 1 )
					{
						if ( skipCnt < skipLimit )
						{
							++skipCnt ;
							continue ;
						}
					}

					if ( size >= 100 && allowTotalSkip )
						continue ;

					skipCnt = 0 ;
					int repeats = size ;
					if ( puse != NULL )
					{
						repeats = 0 ;
						for ( j = 0 ; j < size ; ++j )
						{
							if ( !puse->Get( indexHit[j].idx ) ) 
								continue ;
							++repeats ;
						}
					}

					if ( barcode != -1 )
						repeats = 1 ;

					for ( j = 0 ; j < size ; ++j )
					{
						struct _hit nh ;
						nh.indexHit = indexHit[j] ;
						nh.readOffset = i - kmerLength + 1 ;
						nh.strand = 1 ;
						nh.repeats = repeats ;
						if ( puse != NULL && !puse->Get( indexHit[j].idx ) )
							continue ;
						if ( barcode != -1 && seqs[ indexHit[j].idx ].barcode != barcode )
							continue ;
						hits.PushBack( nh ) ;
					}
				}

				prevKmerCode = kmerCode ;
			}
		}
		// Locate the hits from the opposite-strand case.
		ReverseComplement( rcRead, read, len ) ;		

		if ( strand != 1 )
		{
			kmerCode.Restart() ;
			for ( i = 0 ; i < kmerLength - 1 ; ++i )
				kmerCode.Append( rcRead[i] ) ;

			skipCnt = 0 ; 
			for ( ; i < len ; ++i )
			{
				kmerCode.Append( rcRead[i] ) ;
				if ( downSample > 1 && ( i - kmerLength + 1 ) % downSample != 0 )
					continue ;
				if ( i == kmerLength - 1 || !prevKmerCode.IsEqual( kmerCode ) )
				{
					SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

					int size = indexHit.Size() ;

					if ( size >= 100 && puse == NULL /*&& barcode == -1*/ && i != kmerLength - 1 && i != len - 1 )
					{
						if ( skipCnt < skipLimit )
						{
							++skipCnt ;
							continue ;
						}
					}
					if ( size >= 100 && allowTotalSkip )
						continue ;

					skipCnt = 0 ;
					
					int repeats = size ;
					if ( puse != NULL )
					{
						repeats = 0 ;
						for ( j = 0 ; j < size ; ++j )
						{
							if ( !puse->Get( indexHit[j].idx ) ) 
								continue ;
							++repeats ;
						}
					}

					if ( barcode != -1 )
						repeats = 1 ;

					for ( j = 0 ; j < size ; ++j )
					{
						struct _hit nh ;
						nh.indexHit = indexHit[j] ;
						nh.readOffset = i - kmerLength + 1 ;
						nh.strand = -1 ;
						nh.repeats = repeats ;
						if ( puse != NULL && !puse->Get( indexHit[j].idx ) )
							continue ;
						if ( barcode != -1 && seqs[ indexHit[j].idx ].barcode != barcode )
							continue ;

						hits.PushBack( nh ) ;
						if ( seqs[indexHit[j].idx].name == NULL)
						{
							printf( "%d %d\n", indexHit[j].idx, indexHit[j].offset ) ;
						}
						assert( seqs[indexHit[j].idx].name != NULL ) ;
					}
				}

				prevKmerCode = kmerCode ;
			}
		}
		return hits.Size() ;
	}

	// Obtain the overlaps, each overlap further contains the hits induce the overlap. 
	// readType: 0(default): sqeuencing read. 1:seqs, no need filter.
	// ## forExtension (not used): further filter to ignore regions that could not result in sequence extension.
	// skipRepeats: skip the repeated hits.
	// Return: the number of overlaps.
	int GetOverlapsFromRead( char *read, int strand, int barcode, int readType, bool skipRepeats, std::vector<struct _overlap> &overlaps, 
		SimpleVector<bool> *puse = NULL )
	{
		int i, j, k ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return -1 ;
		
		int overlapCnt = 0 ;
		SimpleVector<struct _hit> hits ;
		char *rcRead =  new char[len + 1] ;
		
		if ( skipRepeats /*&& barcode == -1*/ && puse == NULL )
		{
			GetHitsFromRead( read, rcRead, strand, barcode, true, hits, puse ) ;
			SortHits( hits, true ) ;
			overlapCnt = GetOverlapsFromHits( hits, hitLenRequired, 0, overlaps ) ;

			if ( overlapCnt == 0 )
			{
				hits.Clear() ;
				overlaps.clear() ;
			}
		}
		
		if ( overlapCnt == 0 ) // readType == 1 or failed aggressive skipping for readType 0.
		{
			GetHitsFromRead( read, rcRead, strand, barcode, false, hits, puse ) ;
			SortHits( hits, true ) ;

			// Find the overlaps.

			//if ( seqs.size() == 1 )
			//	for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
			//		printf( "- %d %s %d %d\n", it->readOffset, seqs[ it->indexHit.idx ].name, it->indexHit.offset, it->strand ) ;
			//if ( seqs.size() == 1 )
			//	for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
			//		printf( "- %d %s %d %d\n", it->readOffset, seqs[ it->indexHit.idx ].name, it->indexHit.offset, it->strand ) ;

			//int hitLenRequired = 31 ;
			int filterHits = 0 ;
			if ( readType == 0 )
			{
				//hitLenRequired = ( len / 3 < hitLenRequired ? hitLenRequired : ( len / 3 ) ) ;
				filterHits = 1 ;
			}

			overlapCnt = GetOverlapsFromHits( hits, hitLenRequired, filterHits, overlaps ) ;
		}
		delete[] rcRead ;
		//if ( seqs.size() != 620 )
		//	printf( "overlapCnt = %d\n", overlapCnt ) ;
		/*if ( overlapCnt > 100 )
		{
			printf( "%s\n", read ) ;
			exit( 1 ) ;
		}*/
		//if ( !strcmp( read, "TAGTAATCACTACTGGGCCTGGATCCGCCAGCCCCCAGGGAAAGGGCTGGAGTGGATTGGGAGTATCCATTCTAGTGGGAGCACCTACTTCAACCCGTCCCTCAAGAGTCGAGTCTCCACATCCGTAGACACGTCCGACAATCAAGTCTCCCTGAAGCTGAGGTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTGCGAGACAGTTTCTCCATCTGGACCCCATGTCCAACTGGTTCGACCCCCGG") && filterHits == 0  )
		//for ( i = 0 ; i < overlapCnt ; ++i )
		//	fprintf( stderr, "small test %d: %d %s %d. %d %d %d %d\n", i, overlaps[i].seqIdx,seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ; 
		
		// Determine whether we want to add this reads by looking at the quality of overlap
		if ( overlapCnt == 0 )
		{
			overlapCnt = GetVJOverlapsFromHits( hits, overlaps ) ;
			if ( overlapCnt == 0 )
				return 0 ;
		}
		// Filter out overlaps that is not a real overlap. 
		/*k = 0 ; 
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			int seqLen = strlen( seqs[ overlaps[i].seqIdx ].consensus ) ; 
			// TODO: adjust the boundary effect for V, J gene.
			if ( ( overlaps[i].readStart - radius >= 0 && overlaps[i].readEnd + radius > len - 1 
						&& overlaps[i].seqStart - radius < 0 ) // the last part of the read overlaps with the index
					|| ( overlaps[i].readStart - radius < 0 && overlaps[i].readEnd + radius <= len - 1 
						&& overlaps[i].seqEnd + radius > seqLen - 1 ) // the first part of the read overlaps with the index
					|| ( overlaps[i].readStart - radius < 0 && overlaps[i].readEnd + radius > len - 1 ) // the read is contained.
			   )
			{
				overlaps[k] = overlaps[i] ;
				++k ;
			}
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;*/

		// Since the seqs are all from the same strand, the overlaps should be on the same strand.
		std::sort( overlaps.begin(), overlaps.end() ) ;
		
		//for ( i = 0 ; i < overlapCnt ; ++i )
		//	printf( "%d: %d %s %d. %d %d %d %d. %d\n", i, overlaps[i].seqIdx,seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd, overlaps[i].matchCnt ) ; 
		if ( readType == 0 )
		{
			k = 1 ;
			for ( i = 1 ; i < overlapCnt ; ++i )
			{
				if ( overlaps[i].strand != overlaps[0].strand )
				{
					delete overlaps[i].hitCoords ;
					overlaps[i].hitCoords = NULL ;
					continue ;		
				}
				if ( i != k )
					overlaps[k] = overlaps[i] ;
				++k ;
			}
		}
		else
		{
			k = 0 ;
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				if ( overlaps[i].strand != 1 )
				{
					delete overlaps[i].hitCoords ;
					overlaps[i].hitCoords = NULL ;
					continue ;		
				}
				if ( i != k )
					overlaps[k] = overlaps[i] ;
				++k ;
			}
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		/*for ( i = 1 ; i < overlapCnt ; ++i )
		{
			if ( overlaps[i].strand != overlaps[i - 1].strand )
			{
				overlaps.clear() ;
				for ( i = 0 ; i < overlapCnt ; ++i )
					delete overlaps[i].hitCoords ;
				return 0 ;
			}
		}*/

		rcRead = new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		

		// Compute similarity overlaps
		//if ( overlaps.size() > 1000 )
		/*int bestNovelOverlap = -1 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef )
				continue ;

			if ( bestNovelOverlap == -1 || overlaps[i] < overlaps[ bestNovelOverlap ] )
				bestNovelOverlap = i ;
		}
		if ( bestNovelOverlap != -1 )
		{
			struct _overlap tmp ;
			tmp = overlaps[0] ;
			overlaps[0] = overlaps[ bestNovelOverlap ] ;
			overlaps[ bestNovelOverlap ] = tmp ;
		}*/
		
		int firstRef = -1 ;
		int bestNovelOverlap = -1 ;
		SimpleVector<struct _pair> readOverlapRepresentatives ; // The non-subset best overlaps
		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *r ;
			if ( overlaps[i].strand == 1 )
				r = read ;
			else
				r = rcRead ;

			SimpleVector<struct _pair> &hitCoords = *overlaps[i].hitCoords ; 	
			int hitCnt = hitCoords.Size() ;
			int matchCnt = 0, mismatchCnt = 0, indelCnt = 0  ;
			double similarity = 1 ;
			
			//if ( overlaps[i].seqIdx == 4 )
			//	printf( "%d %d %d. %d %d\n", seqs[ overlaps[i].seqIdx ].isRef, overlapCnt, bestNovelOverlap,
			//			overlaps[i].matchCnt, overlaps[ bestNovelOverlap ].matchCnt ) ;
			// Some fast pre-filters
			if ( seqs[ overlaps[i].seqIdx ].isRef )
			{
				if ( firstRef == -1 )
					firstRef = i ;
				/*else
				{
					if ( overlaps[i].matchCnt < 0.9 * overlaps[ firstRef ].matchCnt )
					{
						overlaps[i].similarity = 0 ;  // No need to look into this.
						continue ;
					}
				}*/
			}
			else if ( bestNovelOverlap != -1 && readType == 0 && overlapCnt > 50 )
			{
				// If we already found a perfect match.
				if ( overlaps[ bestNovelOverlap ].readStart == 0 && overlaps[ bestNovelOverlap ].readEnd == len - 1 ) 
				{
					//printf( "fast perfect filter %d: %lf %d %d\n", overlaps[ bestNovelOverlap ].seqIdx,
					//	overlaps[ bestNovelOverlap ].similarity, overlaps[ bestNovelOverlap ].matchCnt, 
					//	overlaps[i].matchCnt ) ;
					if ( overlaps[ bestNovelOverlap ].similarity == 1 )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
					else if ( overlaps[ bestNovelOverlap ].similarity > repeatSimilarity &&
						overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
				}
				
				// Almost fully cover case
				if ( overlaps[ bestNovelOverlap ].readStart + len - 1 - overlaps[ bestNovelOverlap].readEnd < radius )
				{
					if ( overlaps[ bestNovelOverlap].similarity == 1 
						&& overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
					else if ( ( overlaps[ bestNovelOverlap ].similarity > repeatSimilarity 
							|| isLongSeqSet ) 
						&& overlaps[i].matchCnt < 0.8 * overlaps[ bestNovelOverlap ].matchCnt )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
				}
				
				// Directly filter the bad overlaps if it is inside of a novel seq and worse than the best one. 
				if ( overlaps[i].seqStart - overlaps[i].readStart >= radius 
						&& overlaps[i].seqEnd + ( len - 1 - overlaps[i].readEnd ) + radius < 
							seqs[ overlaps[i].seqIdx ].consensusLen 
						&& overlaps[ bestNovelOverlap ].matchCnt > 0.97 * ( 2 * len ) 
						&& overlaps[ bestNovelOverlap ].similarity > repeatSimilarity 
						&& overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
				{
					//printf( "fast filter\n" ) ;
					overlaps[i].similarity = 0 ;	
					continue ;
				}
				
				// Filter overlaps that is subset of some presentative overlaps. 
				if ( readOverlapRepresentatives.Size() > 0 && isLongSeqSet )
				{
					int size = readOverlapRepresentatives.Size() ;
					for ( j = 0 ; j < size ; ++j )
					{
						k = readOverlapRepresentatives[j].a ;
						if ( overlaps[i].readStart >= overlaps[k].readStart 
								&& overlaps[i].readEnd <= overlaps[k].readEnd 
								&& ( overlaps[i].matchCnt < 0.8 * overlaps[k].matchCnt || 
								IsOverlapSubstringOf( overlaps[i], overlaps[k], true, 1 ) ) )
						
						{
							break ;
						}
					}
					if ( j < size )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
				}

				if ( overlaps[i].matchCnt < 0.4 * overlaps[ bestNovelOverlap ].matchCnt )
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}

				if ( overlapCnt > 1000 && overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}

				//printf( "%d %d: %d %d\n", overlapCnt, overlaps[i].seqIdx, overlaps[i].matchCnt, overlaps[ bestNovelOverlap ].matchCnt ) ;
				
			}
			
			/*if ( bestNovelOverlap != -1 && forExtension )
			{

				int maxGap = 3 * hitLenRequired + 3 * kmerLength ;
				if ( maxGap < 100 )
					maxGap = 100 ;
				
				if ( overlaps[i].readStart + 1 > maxGap && overlaps[i].readEnd + maxGap < len )
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}
				
				int seqLen = seqs[ overlaps[i].seqIdx ].consensusLen ;
				
				// Possible extension to right
				if ( overlaps[i].readStart + 1 <= maxGap 
					&& overlaps[i].seqEnd + maxGap < seqLen )// But could not extend
					//&& len - overlaps[i].readStart > seqLen - overlaps[i].seqStart ) // Make sure not contained in
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}
				// Possible extension to left
				if ( overlaps[i].readEnd + maxGap >= len 
					&& overlaps[i].seqStart + 1 > maxGap )
					//&& overlaps[i].readStart < overlaps[i].seqStart )
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}
			}*/

			matchCnt += 2 * kmerLength ;
			char *align = new char[ overlaps[i].readEnd - overlaps[i].readStart + 1 + 
				overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 1] ;
			for ( j = 1 ; j < hitCnt ; ++j )
			{
				if ( hitCoords[j - 1].b - hitCoords[j - 1].a == hitCoords[j].b - hitCoords[j].a )
				{
					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a )
					{
						matchCnt += 2 * ( hitCoords[j].a - hitCoords[j - 1].a ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ; 
						
						if ( seqs[ overlaps[i].seqIdx ].isRef  )
						{
							//printf( "Use ref %d %d.\n", hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
							//	hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) ) ;
							AlignAlgo::GlobalAlignment( 
								seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;
						}
						else
						{
						//AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
							//printf( "Use novel %d %d.\n", hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
							//	hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) ) ;
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;

							/*if ( seqs.size() == 1 )
							{
								AlignAlgo::VisualizeAlignment( 
										seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
										hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
										r + hitCoords[j - 1].a + kmerLength, 
										hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
										align ) ;

							}*/
						}

						int count[3] ;
						GetAlignStats( align, false, count[0], count[1], count[2] ) ;
						matchCnt += 2 * count[0] ;
						mismatchCnt += count[1] ;
						indelCnt += count[2] ;
						

						if ( ( radius == 0 || !seqs[ overlaps[i].seqIdx ].isRef ) && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}

						// Check whether the mismatch happens within a short window (9). 
						// Note that in this case, there is no indel.
						/*int l ;
						int windowSize = 9 ;
						int windowMaxError = 3 ;
						k = 0 ;
						for ( l = 0 ; align[l] != -1 ; ++l )
						{
							if ( align[l] == EDIT_MISMATCH )
								++k ;
							if ( l >= windowSize && align[l - windowSize] == EDIT_MISMATCH )
								--k ;	
							
							if ( k > windowMaxError )
							{
								similarity = 0 ;
								break ;
							}
						}*/
					}
				}
				else
				{
					if ( radius == 0 || !seqs[ overlaps[i].seqIdx ].isRef )
					{
						similarity = 0 ;
						break ;
					}

					//printf( "%d %d=>%d %d\n", hitCoords[j - 1].a, hitCoords[j - 1].b, hitCoords[j].a, hitCoords[j].b ) ;
					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 < hitCoords[j].b )
					{
						matchCnt += 2 * ( hitCoords[j].a - hitCoords[j - 1].a ) ; //+ kmerLength ;
						// Make the two kmer hit match on coordinate.
						indelCnt += ( hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) + 
							( hitCoords[j].a + kmerLength - hitCoords[j - 1].a )  ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 < hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						//matchCnt += kmerLength + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						matchCnt += 2 * ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ( hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) +
							( hitCoords[j].b + kmerLength - hitCoords[j - 1].b ) ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a &&
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						//matchCnt += ( hitCoords[j].a - hitCoords[j - 1].a ) + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						matchCnt += 2 * MIN( hitCoords[j].a - hitCoords[j - 1].a, hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ABS( ( hitCoords[j].a - hitCoords[j].b ) - 
							( hitCoords[j - 1].a - hitCoords[j - 1].b ) ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ;
						 
						if ( seqs[ overlaps[i].seqIdx ].isRef )
						{
							//printf( "Use ref2 %d %d.\n", hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
							//	hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) ) ;
							AlignAlgo::GlobalAlignment( 
								seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) , 
								align ) ;	
						}
						else
						{
							//AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) , 
								align ) ;	
						}
						
						int count[3] ;
						GetAlignStats( align, false, count[0], count[1], count[2] ) ;
						matchCnt += 2 * count[0] ;
						mismatchCnt += count[1] ;
						indelCnt += count[2] ;
						if ( !seqs[ overlaps[i].seqIdx ].isRef && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}

					}
				}
			} // for j
			delete[] align ;
			
			//printf( "%d: %d %d %d %lf\n", overlaps[i].seqIdx, matchCnt, overlaps[i].seqEnd - overlaps[i].seqStart + 1, overlaps[i].readEnd - overlaps[i].readStart + 1, similarity ) ;
			overlaps[i].matchCnt = matchCnt ;
			if ( similarity == 1 )
				overlaps[i].similarity = (double)matchCnt / ( overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 
								overlaps[i].readEnd - overlaps[i].readStart + 1 ) ;
			else
				overlaps[i].similarity = 0 ;
			
			if ( IsOverlapLowComplex( r, overlaps[i]) )
				overlaps[i].similarity = 0 ;
			

			//printf( "%d: %d %d %d %lf\n", overlaps[i].seqIdx, matchCnt, overlaps[i].seqEnd - overlaps[i].seqStart + 1, overlaps[i].readEnd - overlaps[i].readStart + 1, similarity ) ;
			overlaps[i].matchCnt = matchCnt ;
			if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity > 0 )
			{
				if ( bestNovelOverlap == -1 || overlaps[i] < overlaps[ bestNovelOverlap ] ) // the less than means has higher priority
				{
					bestNovelOverlap = i ;
				}
			}
			/*if ( overlaps[i].similarity > 0 )
			{
				printf( "%d: %d %d %d %d %d %lf\n", matchCnt, overlaps[i].seqIdx, overlaps[i].readStart, overlaps[i].readEnd, 
							overlaps[i].seqStart, overlaps[i].seqEnd, similarity ) ;
			}
			assert( overlaps[i].similarity <= 1 ) ;*/

			if ( !seqs[ overlaps[i].seqIdx ].isRef && readType == 1 
				 && overlaps[i].similarity < novelSeqSimilarity )
			{
				// Look for the core part that with high identity.
				// since similarity is greater than 0, we already know there is no indel in the hits.
				int maxLen = 0 ;
				int maxS = 0, maxE = 0 ;
				for ( j = 0 ; j < hitCnt ;  )
				{
					for ( k = j + 1 ; k < hitCnt ; ++k )
					{
						if ( hitCoords[k].a > hitCoords[k - 1].a + kmerLength - 1 
							|| hitCoords[k].a - hitCoords[k].b != hitCoords[k-1].a - hitCoords[k-1].b )	
							break ;
					}

					int len = hitCoords[k - 1].a - hitCoords[j].a + kmerLength ;
					if ( len > maxLen )
					{
						maxLen = len ;
						maxS = j ;
						maxE = k - 1 ;
					}
					j = k ;
				}
				
				if ( maxLen >= hitLenRequired )
				{
					overlaps[i].readStart = hitCoords[maxS].a ;
					overlaps[i].readEnd = hitCoords[maxE].a + kmerLength - 1 ;
					overlaps[i].seqStart = hitCoords[ maxS ].b ;
					overlaps[i].seqEnd = hitCoords[ maxE ].b + kmerLength - 1 ;
					overlaps[i].similarity = 1.0 ;
					overlaps[i].matchCnt = 2 * maxLen ;
				}
			}

			if ( overlaps[i].similarity > 0 )
			{
				int size = readOverlapRepresentatives.Size() ;
				for ( j = 0 ; j < size ; ++j )
				{
					int k = readOverlapRepresentatives[j].a ;
					if ( overlaps[i].readStart >= overlaps[k].readStart 
						&& overlaps[i].readEnd <= overlaps[k].readEnd )
						break ;		
				}

				if ( j >= size )
				{
					struct _pair np ;
					np.a = i ;
					np.b = 0 ;
					readOverlapRepresentatives.PushBack( np ) ;
				}
			}
		} // for i
		delete[] rcRead ;

		// Release the memory for hitCoords.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			overlaps[i].hitCoords->Release() ;
			delete overlaps[i].hitCoords ;
			overlaps[i].hitCoords = NULL ;
		}

		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < refSeqSimilarity )
				continue ;
			else if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < novelSeqSimilarity )
				continue ;

			//printf( "%d %s: %d-%d %d-%d: %d %d %lf\n", overlaps[i].seqIdx, seqs[overlaps[i].seqIdx].name, 
			//	overlaps[i].readStart, overlaps[i].readEnd, 
			//	overlaps[i].seqStart, overlaps[i].seqEnd, overlaps[i].matchCnt, overlaps[i].strand, overlaps[i].similarity ) ;
			overlaps[k] = overlaps[i] ;
			++k ;
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		
		//printf( "return: %d\n", overlapCnt) ;
		return overlapCnt ;
	}
	
	// Figure out whether a seq is a (almost) substring of another seq.
	int BuildSeqSubstringRelation( std::vector<struct _overlap> &subsetOf )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		for ( k = 0 ; k < seqCnt ; ++k )
			subsetOf[k].seqIdx = -1 ;
		
		std::map<int, int> seqHitCnt ;
		SimpleVector<struct _pair> firstSeqHit ; 
		firstSeqHit.ExpandTo( seqCnt ) ;

		for ( k = 0 ; k < seqCnt ; ++k )
		{
			if ( seqs[k].consensus == NULL )
				continue ;

			char *consensus = seqs[k].consensus ;
			int len = seqs[k].consensusLen ;
			if ( len < kmerLength )
				return -1 ;

			//SimpleVector<struct _hit> hits ;		    	

			KmerCode kmerCode( kmerLength ) ;
			KmerCode prevKmerCode( kmerLength ) ;

			//int skipLimit = 3 ;
			int skipLimit = kmerLength / 2 ; 

			int skipCnt = 0 ;
			int hitCnt ;
			for ( i = 0 ; i < kmerLength - 1 ; ++i )
				kmerCode.Append( consensus[i] ) ;
			seqHitCnt.clear() ;
			hitCnt = 0 ;
			for ( ; i < len ; ++i )
			{
				kmerCode.Append( consensus[i] ) ;
				if ( i == kmerLength || !prevKmerCode.IsEqual( kmerCode ) )
				{
					SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

					int size = indexHit.Size() ;
					if ( size >= 100 )
					{
						if ( skipCnt < skipLimit )
						{
							++skipCnt ;
							continue ;
						}
					}

					skipCnt = 0 ;

					for ( j = 0 ; j < size ; ++j )
					{
						struct _hit nh ;
						if ( indexHit[j].idx == k )
							continue ;
												
						if ( seqHitCnt.find( indexHit[j].idx ) != seqHitCnt.end() )
						{
							if ( hitCnt >= 50 && seqHitCnt[ indexHit[j].idx ] < hitCnt * 0.5 )
							{
								seqHitCnt.erase( indexHit[j].idx ) ;
							}
							else
								++seqHitCnt[ indexHit[j].idx ] ;
						}
						else if ( hitCnt < 50 )
						{
							seqHitCnt[ indexHit[j].idx ] = 1 ;
							firstSeqHit[ indexHit[j].idx ].a = i - kmerLength + 1 ;
							firstSeqHit[ indexHit[j].idx ].b = indexHit[j].offset ;
						}
					}
					++hitCnt ;
				}

				prevKmerCode = kmerCode ;
			}
			
			for ( std::map<int, int>::iterator it = seqHitCnt.begin() ; it != seqHitCnt.end() ; ++it )
			{
				//printf( "%d %d\n", it->second, hitCnt ) ;
				if ( it->second < hitCnt * 0.6 )
					continue ;
				int seqIdx = it->first ; 
				// Test whether k is a substring of seqIdx
				if ( firstSeqHit[ seqIdx ].b - firstSeqHit[ seqIdx ].a < 0 )
					continue ;

				int start = firstSeqHit[ seqIdx ].b - firstSeqHit[ seqIdx ].a ;
				if ( start + seqs[k].consensusLen - 1 >= seqs[seqIdx].consensusLen )
					continue ;
				int matchCnt = 0 ;
				int mismatchCnt = 0 ;
				int l ;
				for ( j = 0, l = start ; j < seqs[k].consensusLen ; ++j, ++l )
				{
					if ( seqs[k].consensus[j] != seqs[seqIdx].consensus[l] )
						++mismatchCnt ;
					else
						++matchCnt ;

					if ( mismatchCnt >= 2 )
						break ;
				}
				//char *p = strstr( seqs[ it->first ].consensus, consensus ) ;
				//printf( "test %d\n%s\n%s\n", mismatchCnt, seqs[k].consensus, seqs[ seqIdx ].consensus  ) ;
				if ( mismatchCnt < 2 ) // some mismatch are allowed because we allow mismatch in the overlaps of branch graph.
				{
					subsetOf[k].seqIdx = it->first ;
					subsetOf[k].readStart = 0 ;
					subsetOf[k].readEnd = seqs[k].consensusLen - 1 ;
					subsetOf[k].seqStart = start ;
					subsetOf[k].seqEnd = subsetOf[k].seqStart + seqs[k].consensusLen - 1 ; 
					break ;
				}
			}
		}

		return 0 ;
	}


	// adj is the adjacent list for each seq. The size of the adj array is seqCnt.
	int BuildSeqOverlapGraph( int overlapLength, std::vector<struct _overlap> *adj )
	{
		// Build overlap graph.
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int maxLen = 0 ;
		char *align ;
		for ( i = 0 ; i < seqCnt ; ++i )
			if ( seqs[i].consensusLen > maxLen )
				maxLen = seqs[i].consensusLen ;
		align = new char[2 * maxLen + 2] ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL ) 
				continue ;

			std::vector<struct _overlap> overlaps ;
			int overlapCnt ;

			overlapCnt = GetOverlapsFromRead( seqs[i].consensus, 1, seqs[i].barcode, 1, false, overlaps ) ;
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( overlaps[j].strand == -1 ) // Note that all the seqs are from 5'->3' on its strand.
					continue ;

				if ( i == overlaps[j].seqIdx  )
					continue ;
				struct _overlap extendedOverlap ;

				if ( ExtendOverlap( seqs[i].consensus, seqs[i].consensusLen, seqs[ overlaps[j].seqIdx ], 1.0,  
					align, overlaps[j], extendedOverlap ) == 1 ) 
				{
					if ( extendedOverlap.readEnd - extendedOverlap.readStart + 1 >= overlapLength && 
						( extendedOverlap.readStart > 0 || // i before j
						( extendedOverlap.readStart == 0 && extendedOverlap.readEnd == seqs[i].consensusLen - 1 )// i contained in j 
						) ) 
					{
						if ( ( extendedOverlap.readStart == 0 && extendedOverlap.readEnd == seqs[i].consensusLen - 1 ) 
							&& ( extendedOverlap.seqStart == 0 && 
								extendedOverlap.seqEnd == seqs[ overlaps[j].seqIdx ].consensusLen - 1 ) )
						{
							// The case of two almost identical seqs.
							if ( i < overlaps[j].seqIdx )
								continue ;
						}
						adj[i].push_back( extendedOverlap ) ;	
					}
				}
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
			std::sort( adj[i].begin(), adj[i].end() ) ;
		delete[] align ;
		return 1 ;
	}
	
	// Only build the graph for the seqs with id marked true in "use" array.
	int BuildBranchGraph( std::vector<struct _overlap> *adj, int leastOverlapLen,
		std::vector<struct _overlap> *prevAdj = NULL, std::vector<struct _overlap> *nextAdj = NULL )
	{
		// Build overlap graph.
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int maxLen = 0 ;
		char *align ;
		for ( i = 0 ; i < seqCnt ; ++i )
			if ( seqs[i].consensusLen > maxLen )
				maxLen = seqs[i].consensusLen ;
		align = new char[2 * maxLen + 2] ;
		SimpleVector<bool> use ; // Buffer to represent whether to use the seq when computing overlaps
		use.ExpandTo( seqCnt ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
			use[i] = false ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL ) 
				continue ;
			//if ( use[i] == false )
			//	continue ;
			if ( prevAdj != NULL && nextAdj != NULL && prevAdj[i].size() + nextAdj[i].size() == 0 )
				continue ;

			// Update the use buffer 
			if ( prevAdj != NULL && nextAdj != NULL )
			{
				int prevCnt = prevAdj[i].size() ;
				int nextCnt = nextAdj[i].size() ;

				for ( k = 0 ; k < prevCnt ; ++k )
					use[ prevAdj[i][k].seqIdx ] = true ;

				for ( k = 0 ; k < nextCnt ; ++k )
					use[ nextAdj[i][k].seqIdx ] = true ;
			}


			std::vector<struct _overlap> overlaps ;
			int overlapCnt ;

			double backupSimilarity = novelSeqSimilarity ;
			novelSeqSimilarity = repeatSimilarity ;
			overlapCnt = GetOverlapsFromRead( seqs[i].consensus, 1, seqs[i].barcode, 1, false, overlaps, &use ) ;
			novelSeqSimilarity = backupSimilarity ;

			//printf( "%d %d\n", i, overlapCnt ) ;
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				//printf( "%d: %d %d %d %d\n", i, j, overlaps[j].seqIdx, overlaps[j].readStart, overlaps[j].readEnd ) ;
				if ( overlaps[j].strand == -1 ) // Note that all the seqs are from 5'->3' on its strand.
					continue ;

				if ( i == overlaps[j].seqIdx || use[ overlaps[j].seqIdx ] == false )
					continue ;

				struct _overlap extendedOverlap ;
				int addMatchCnt = 0 ;				
				// Locally extend the overlap. In this case, the overlap actually just means share.
				// Allow "radius" overhang.
				int seqIdx = overlaps[j].seqIdx ;

				int a, b ;
				int matchCnt = 0 ;
				int rightExtend = 0;
				int rightExtendMatchCnt = 0 ; 
				for ( k = 1, a = overlaps[j].readEnd + 1, b = overlaps[j].seqEnd + 1 ; 
						a < seqs[i].consensusLen && b < seqs[ seqIdx ].consensusLen ; ++a, ++b, ++k )
				{
					if ( IsPosWeightCompatible( seqs[i].posWeight[a], seqs[ seqIdx ].posWeight[b] ) )
					{
						++matchCnt ;

						if ( matchCnt > k * 0.75 )
						{
							rightExtendMatchCnt = 2 * matchCnt ;
							rightExtend = k ;
						}
					}
				}

				matchCnt = 0 ;
				int leftExtend = 0 ;
				int leftExtendMatchCnt = 0 ;
				for ( k = 1, a = overlaps[j].readStart - 1, b = overlaps[j].seqStart - 1 ; 
						a >= 0 && b >= 0 ; --a, --b, ++k )
				{
					/*if ( b >= seqs[ seqIdx ].consensusLen )
					  fprintf( stderr, "%d %d %d %d\n%s\n%s\n", overlaps[j].readStart, overlaps[j].readEnd,
					  overlaps[j].seqStart, overlaps[j].seqEnd,
					  seqs[i].consensus, seqs[ seqIdx ].consensus ) ;*/
					if ( IsPosWeightCompatible( seqs[i].posWeight[a], seqs[ seqIdx ].posWeight[b] ) )
					{
						++matchCnt ;

						if ( matchCnt > k * 0.75 )
						{
							leftExtendMatchCnt = 2 * matchCnt ;
							leftExtend = k ;
						}
					}
				}

				extendedOverlap = overlaps[j] ;
				extendedOverlap.readStart -= leftExtend ;
				extendedOverlap.seqStart -= leftExtend ;
				extendedOverlap.readEnd += rightExtend ;
				extendedOverlap.seqEnd += rightExtend ;
				extendedOverlap.matchCnt += rightExtendMatchCnt + leftExtendMatchCnt ;
				extendedOverlap.similarity = (double)extendedOverlap.matchCnt / 
					( extendedOverlap.readEnd - extendedOverlap.readStart + 1 +
					  extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ) ;
				//printf( "%d %d %d %lf\n", extendedOverlap.seqIdx, extendedOverlap.readStart,
				//	extendedOverlap.readEnd, extendedOverlap.similarity ) ;
				if ( extendedOverlap.readEnd - extendedOverlap.readStart + 1 < leastOverlapLen )
					continue ;

				//printf( "hahaha %d %d\n", extendedOverlap.readEnd - extendedOverlap.readStart + 1, extendedOverlap.matchCnt ) ;
				if ( extendedOverlap.similarity < repeatSimilarity )
					extendedOverlap = overlaps[j] ;
				if ( extendedOverlap.similarity >= repeatSimilarity ) 
				{
					adj[i].push_back( extendedOverlap ) ;
#ifdef DEBUG
					printf( "branch %d %d: %d %d %d %d %d %lf\n", i, j, extendedOverlap.seqIdx, 
							extendedOverlap.readStart, extendedOverlap.readEnd,
							extendedOverlap.seqStart, extendedOverlap.seqEnd,
							extendedOverlap.similarity ) ;
#endif				
				}	
			}
			
			// Reset the use buffer 
			if ( prevAdj != NULL && nextAdj != NULL )
			{
				int prevCnt = prevAdj[i].size() ;
				int nextCnt = nextAdj[i].size() ;

				for ( k = 0 ; k < prevCnt ; ++k )
					use[ prevAdj[i][k].seqIdx ] = false ;

				for ( k = 0 ; k < nextCnt ; ++k )
					use[ nextAdj[i][k].seqIdx ] = false ;
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
			std::sort( adj[i].begin(), adj[i].end() ) ;
		delete[] align ;
		return 1 ;
	}


	void UpdatePosWeightFromRead( SimpleVector<struct _posWeight> &posWeight, int offset, char *read )
	{
		int i ;
		for ( i = 0 ; read[i] ; ++i )
		{
			if ( read[i] != 'N' )
				++posWeight[i + offset].count[ nucToNum[ read[i] - 'A' ] ] ;
		}
	}
public:
	SeqSet( int kl ) 
	{
		kmerLength = kl ;
		radius = 10 ;
		hitLenRequired = 31 ;
		isLongSeqSet = false ;

		novelSeqSimilarity = 0.9 ;
		refSeqSimilarity = 0.75 ; 
		repeatSimilarity = 0.95 ; 

		gapN = 7 ;

		prevAddInfo.readStart = -1 ;
	}
	~SeqSet() 
	{
		int size ;
		int i ;
		size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].consensus != NULL )
				free( seqs[i].consensus ) ;
			if ( seqs[i].name != NULL )
				free( seqs[i].name ) ;	
		}
	}

	int Size()
	{
		return seqs.size() ;
	}

	int SetRadius( int r )  
	{
		return radius = r ;
	}
	
	int SetHitLenRequired( int l )
	{
		return hitLenRequired = l ;
	}

	double SetNovelSeqSimilarity( double s )
	{
		return novelSeqSimilarity = s ;
	}

	void ReverseComplement( char *rcSeq, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
		{
			if ( seq[len - 1 - i] != 'N' )
				rcSeq[i] = numToNuc[ 3 - nucToNum[seq[len - 1 - i] - 'A'] ];
			else
				rcSeq[i] = 'N' ;
		}
		rcSeq[i] = '\0' ;
	}
	
	// Input some baseline sequence to match against.
	void InputRefFa( char *filename, bool isIMGT = false ) 
	{
		int i, j, k ;
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		
		KmerCode kmerCode( kmerLength ) ;
		while ( fa.Next() )
		{
			// Insert the kmers 
			struct _seqWrapper ns ;
			// Filter the gene with "/" sign
			for ( i = 0 ; fa.id[i] ; ++i )
				if ( fa.id[i] == '/' )
					break ;
			if ( fa.id[i] == '/' )
				continue ;

			ns.name = strdup( fa.id ) ;
			ns.isRef = true ;

			int id = seqs.size() ;
			seqs.push_back( ns ) ;

			struct _seqWrapper &sw = seqs[id] ;
			int seqLen = strlen( fa.seq ) ;
			sw.consensus = strdup( fa.seq ) ;	
					
			
			// Remove "." from IMGT annotation.
			k = 0 ;
			for ( i = 0 ; i < seqLen ; ++i )
				if ( sw.consensus[i] != '.' )
				{
					sw.consensus[k] = sw.consensus[i] ;
					++k ;
				}
			sw.consensus[k] = '\0' ;

			// Use IMGT documented coordinate to infer CDR1,2,3 coordinate.
			if ( isIMGT && GetGeneType( fa.id ) == 0 && seqLen >= 66 * 3 )
			{	
				// Infer the coordinate for CDR1
				for ( i = 0, k = 0 ; i < 3 * ( 27 - 1 ) ; ++i )
					if ( fa.seq[i] != '.' )		
						++k ;
				sw.info[0].a = k ;
				for ( ; i < 3 * ( 38 ) ; ++i )
					if ( fa.seq[i] != '.' )	
						++k ;
				sw.info[0].b = k - 1 ;
				if (sw.info[0].a > sw.info[0].b) // in case the annotation does not have information on the cdr1
					sw.info[0].a = sw.info[0].b = -1 ;

				// Infer the coordinate for CDR2
				for ( ; i < 3 * ( 56 - 1 ) ; ++i )
					if ( fa.seq[i] != '.' )		
						++k ;
				sw.info[1].a = k ;
				for ( ; i < 3 * ( 65 ) ; ++i )
					if ( fa.seq[i] != '.' )	
						++k ;
				sw.info[1].b = k - 1 ;
				if (sw.info[1].a > sw.info[1].b) // in case the annotation does not have information on the cdr2
					sw.info[1].a = sw.info[1].b = -1 ;
				
				if ( seqLen >= 3 * ( 104 - 1 ) + 1 )
				{
					for ( ; i < 3 * ( 104 - 1 ) ; ++i )
						if ( fa.seq[i] != '.' )
							++k ;
					sw.info[2].a = sw.info[2].b = k ;
				}
				else
					sw.info[2].a = sw.info[2].b = -1 ;
			}
			else if ( isIMGT && GetGeneType( fa.id ) == 2 ) // Found the end position for CDR3
			{
				bool found = false ;
				for ( i = 0 ; i + 11 < seqLen ; ++i )
				{
					if ( ( DnaToAa( fa.seq[i], fa.seq[i + 1], fa.seq[i + 2] ) == 'W' || 
								DnaToAa( fa.seq[i], fa.seq[i + 1], fa.seq[i + 2] ) == 'F' )  
							&& DnaToAa( fa.seq[i + 3], fa.seq[i + 4], fa.seq[i + 5] ) == 'G' 
							&& DnaToAa( fa.seq[i + 9], fa.seq[i + 10], fa.seq[i + 11] ) == 'G' )
					{
						found = true ;
						break ;
					}
				}

				if ( found )
					sw.info[2].a = sw.info[2].b = i ;
				else
					sw.info[2].a = sw.info[2].b = -1 ;
			}
			else
			{
				sw.info[0].a = sw.info[0].b = -1 ;
				sw.info[1].a = sw.info[1].b = -1 ;
				sw.info[2].a = sw.info[2].b = -1 ;
			}


			sw.consensusLen = strlen( sw.consensus );
			sw.barcode = -1 ;
			seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, sw.consensusLen, id ) ;
		}
	}
	
	void InputNovelFa( char *filename ) 
	{
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		
		while ( fa.Next() )
			InputNovelRead( fa.id, fa.seq, 1, -1 ) ;
	}

	int InputRefSeq( char *id, char *read )
	{
		struct _seqWrapper ns ;
		ns.name = strdup( id ) ;
		ns.isRef = true ;

		int seqIdx = seqs.size() ;
		seqs.push_back( ns ) ;

		struct _seqWrapper &sw = seqs[ seqIdx ] ;
		int seqLen = strlen( read ) ;
		sw.consensus = strdup( read ) ;
		sw.consensusLen = seqLen ;	
		int i ;
		sw.posWeight.ExpandTo( sw.consensusLen ) ;
		sw.posWeight.SetZero( 0, sw.consensusLen ) ;
		for ( i = 0 ; i < sw.consensusLen ; ++i )
		{
			if ( sw.consensus[i] != 'N' )
				sw.posWeight[i].count[ nucToNum[ sw.consensus[i] - 'A' ] ] = 1 ;
		}
		KmerCode kmerCode( kmerLength ) ;
		seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, seqLen, seqIdx ) ;
		
		sw.barcode = -1 ;
		sw.minLeftExtAnchor = sw.minRightExtAnchor = 0 ;
		SetPrevAddInfo( seqIdx, 0, seqLen - 1, 0, seqLen - 1, 1 ) ;
		return seqIdx ;
	}	
	
	int InputNovelRead( const char *id, char *read, int strand, int barcode )
	{
		struct _seqWrapper ns ;
		ns.name = strdup( id ) ;
		ns.isRef = false ;
		
		int seqIdx = seqs.size() ;
		seqs.push_back( ns ) ;

		struct _seqWrapper &sw = seqs[ seqIdx ] ;
		int seqLen = strlen( read ) ;
		sw.consensus = strdup( read ) ;
		if ( strand == -1 )
			ReverseComplement( sw.consensus, read, seqLen ) ;
		sw.consensusLen = seqLen ;	
		int i ;
		sw.posWeight.ExpandTo( sw.consensusLen ) ;
		sw.posWeight.SetZero( 0, sw.consensusLen ) ;
		for ( i = 0 ; i < sw.consensusLen ; ++i )
		{
			if ( sw.consensus[i] != 'N' )
				sw.posWeight[i].count[ nucToNum[ sw.consensus[i] - 'A' ] ] = 1 ;
		}
		KmerCode kmerCode( kmerLength ) ;
		seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, seqLen, seqIdx ) ;
	
		sw.minLeftExtAnchor = sw.minRightExtAnchor = 0 ;
		sw.barcode = barcode ;
		SetPrevAddInfo( seqIdx, 0, seqLen - 1, 0, seqLen - 1, strand ) ;


#ifdef DEBUG
		printf( "add novel seq: %d\n", seqIdx ) ;
#endif
		return seqIdx ;
	}

	int InputNovelSeq( char *name, char *seq, SimpleVector<struct _posWeight> &posWeight )
	{
		struct _seqWrapper ns ;
		ns.name = strdup( name ) ;
		ns.isRef = false ;
		int seqIdx = seqs.size() ;
		seqs.push_back( ns ) ;

		struct _seqWrapper &sw = seqs[ seqIdx ] ;
		int seqLen = strlen( seq ) ;
		sw.consensus = strdup( seq ) ;
		sw.consensusLen = seqLen ;	
		sw.posWeight = posWeight ;
		KmerCode kmerCode( kmerLength ) ;
		seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, seqLen, seqIdx ) ;
		sw.barcode = -1 ;	
		sw.minLeftExtAnchor = sw.minRightExtAnchor = 0 ;
		SetPrevAddInfo( seqIdx, 0, seqLen - 1, 0, seqLen - 1, 1 ) ;

		return seqIdx ;
	}

	// Input sequence from another SeqSet
	void InputSeqSet( const SeqSet &in, bool inputRef )
	{
		int i, k ;
		int seqCnt = in.seqs.size() ;
		KmerCode kmerCode( kmerLength ) ;
		
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( !inputRef && in.seqs[i].isRef )
				continue ;
			if ( in.seqs[i].consensus == NULL )
				continue ;
			
			struct _seqWrapper ns ;
			ns.name = strdup( in.seqs[i].name ) ;
			ns.isRef = in.seqs[i].isRef ;
			ns.barcode = in.seqs[i].barcode ;
			
			int id = seqs.size() ;
			ns.consensus = strdup( in.seqs[i].consensus ) ;
			ns.consensusLen = in.seqs[i].consensusLen ;
			ns.posWeight = in.seqs[i].posWeight ;
			ns.minLeftExtAnchor = in.seqs[i].minLeftExtAnchor ;
			ns.minRightExtAnchor = in.seqs[i].minRightExtAnchor ;
			seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, ns.consensusLen, id ) ;
			seqs.push_back( ns ) ;
		}
	}
	
	// Test whether the read share a kmer hit on the seqs.
	bool HasHitInSet( char *read, char *rcRead )
	{
		int i, k ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return false ;

		SimpleVector<struct _hit> hits ;	
		GetHitsFromRead( read, rcRead, 0, -1, false, hits, NULL  ) ;

		int hitCnt = hits.Size() ;
		if ( hitCnt == 0 )
			return false ;

		// Bucket sort.
		int seqCnt = seqs.size() ;
		SimpleVector<struct _hit> *buckets[2] ;
		buckets[0] = new SimpleVector<struct _hit>[seqCnt] ;
		buckets[1] = new SimpleVector<struct _hit>[seqCnt] ;

		for ( i = 0 ; i < hitCnt ; ++i )
		{
			int tag = hits[i].strand == 1 ? 1 : 0 ;
			buckets[tag][ hits[i].indexHit.idx ].PushBack( hits[i] ) ;
		}
		
		// Find the best bucket.
		int max = -1 ;
		int maxTag = -1 ;
		int maxSeqIdx = -1 ;
		for ( k = 0 ; k <= 1 ; ++k )
		{
			for ( i = 0 ; i < seqCnt ; ++i )
			{
				int size = buckets[k][i].Size() ;
				if ( size > 0 && size > max )
				{
					maxTag = k ;
					maxSeqIdx = i ;
					max = size ;
				}
			}
		}
		
		std::vector<struct _overlap> overlaps ;
		GetOverlapsFromHits( buckets[maxTag][maxSeqIdx], hitLenRequired, 1, overlaps ) ;
		delete[] buckets[0] ;
		delete[] buckets[1] ;
		//printf( "%d %d\n", hitCnt, overlaps.size() ) ;	
		for ( i = 0 ; i < overlaps.size() ; ++i )
		{
			//printf( "%s\n", seqs[ overlaps[i].seqIdx ].name ) ;
			delete overlaps[i].hitCoords ;
		}
		if ( overlaps.size() == 0 )
			return false ;
		/*printf( "%s %d %d %lf\n", seqs[ overlaps[0].seqIdx ].name, overlaps[0].readStart, overlaps[0].readEnd, overlaps[0].similarity ) ;
		for ( i = 0 ; i < overlaps[0].hitCoords->Size() ; ++i )
			printf( "%d %d\n", overlaps[0].hitCoords->Get(i).a, overlaps[0].hitCoords->Get(i).b ) ;*/
		return true ;
	}

	// Compute the length of hit from the read, take the overlaps of kmer into account 
	int GetTotalHitLengthOnRead( SimpleVector<struct _hit> &hits )
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

	int GetTotalHitLengthOnSeq( SimpleVector<struct _hit> &hits )
	{
		int hitSize = hits.Size() ;
		int i, j ;
		int ret = 0 ;

		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].indexHit.offset > hits[j - 1].indexHit.offset + kmerLength - 1 )
					break ;
			ret += hits[j - 1].indexHit.offset - hits[i].indexHit.offset + kmerLength ;
			i = j ;
		}
		return ret ;
	}

	// b comes after a, test wheter their gene names is compatible
	bool IsNameCompatible( char *a, char *b )
	{
		int i, j ;
		int maxA = -1 ;
		int minB = 10 ;
		for ( i = 0 ; a[i] ; )	
		{
			if ( a[i] == '+' )
			{
				++i ;
				continue ;
			}

			for ( j = i ; a[j] && a[j] != '+' ; ++j )
				;
			char tmp = a[j] ;
			a[j] = '\0' ;
			int gt = GetGeneType( a + i ) ;
			if ( gt > maxA )
				maxA = gt ;
			a[j] = tmp ;

			i = j ;
		}
		
		for ( i = 0 ; b[i] ; )	
		{
			if ( b[i] == '+' )
			{
				++i ;
				continue ;
			}

			for ( j = i ; b[j] && b[j] != '+' ; ++j )
				;
			char tmp = b[j] ;
			b[j] = '\0' ;
			int gt = GetGeneType( b + i ) ;
			if ( gt < minB )
				minB = gt ;
			b[j] = tmp ;

			i = j ;
		}

		if ( maxA <= minB )
			return true ;
		else
			return false ;
	}
	

	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set. 
	//	   -1: not add. -2: only overlapped with novel seq and could not be extended.
	int AddRead( char *read, char *geneName, int &strand, int barcode, int minKmerCount, bool repetitiveData, double similarityThreshold )
	{
		//printf( "%s\n", read ) ;
		int i, j, k ;
		int len = strlen( read ) ;

		SetPrevAddInfo( -1, -1, -1, -1, -1, 0 ) ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;
		
		overlapCnt = GetOverlapsFromRead( read, strand, barcode, 0, repetitiveData, overlaps ) ;
				
		if ( overlapCnt <= 0 )
			return -1 ;

#ifdef DEBUG
		printf( "geneName: %s\n", geneName ) ;
#endif
		if ( geneName[0] != '\0' )
		{
			k = 0 ;
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				j = 3 ;
				if ( seqs[ overlaps[i].seqIdx ].name[0] >= 'A' && seqs[ overlaps[i].seqIdx ].name[0] <= 'Z' )
				{	
					for ( j = 0 ; j < 3 ; ++j )
						if ( seqs[ overlaps[i].seqIdx ].name[j] != geneName[j] )
							break ;
				}
				if ( j == 3 )
				{
					overlaps[k] = overlaps[i] ;
					++k ;
				}
				/*else
				{
					printf( "haha: %s %s\n%s\n", geneName, seqs[ overlaps[i].seqIdx ].name, read ) ;
				}*/
			}
			overlaps.resize( k ) ;
			overlapCnt = k ;
			
			if ( overlapCnt <= 0 )
				return -1 ;
		}

		std::sort( overlaps.begin(), overlaps.end() ) ;

#ifdef DEBUG
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			printf( "%d: %d %d %s. %d(%d %d). %d %d %d %d. %lf.\n", i, overlaps[i].seqIdx, seqs[ overlaps[i].seqIdx ].consensusLen, 
					seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, 
					seqs[ overlaps[i].seqIdx ].minLeftExtAnchor, seqs[ overlaps[i].seqIdx ].minRightExtAnchor,
					overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd, 
					overlaps[i].similarity ) ; 
			printf( "%s\n", seqs[  overlaps[i].seqIdx ].consensus ) ;
			//if ( !seqs[ overlaps[i].seqIdx ].isRef )
			/*if ( !strcmp( read, "TTACTGTAATATACGATATTTTGACTGGTTATTAAGAGGCGACCCAAGAATCAATACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCT" ) )
			{
				printf( " %s\n",seqs[  overlaps[i].seqIdx ].consensus ) ;
				if ( i == 1 )
				{
					printf( "%d %d\n", seqs[ overlaps[i].seqIdx ].posWeight[103].count[0], 
						seqs[ overlaps[i].seqIdx ].posWeight[103].count[2] ) ;
				}
			}*/
		}
		fflush( stdout ) ;	
#endif		
		// If the read only overlaps with the reference, we will add that to the seq.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( !seqs[ overlaps[i].seqIdx ].isRef )
				break ;
		}
		
		struct _overlap *extendedOverlaps = new struct _overlap[ overlapCnt ];
		struct _overlap *failedExtendedOverlaps = new struct _overlap[ overlapCnt ];
		k = 0 ;
		int ret = -1 ;
		bool addNew = true ;
		int failedExtendedOverlapsCnt = 0 ;
		struct _overlap goodExtendedOverlap ;
		
		goodExtendedOverlap.seqIdx = -1 ;
		if ( i < overlapCnt )
		{
			// Incorporate to existing sequence.
			char *align = new char[3 * len] ;
			char *rcRead = strdup( read ) ;
			ReverseComplement( rcRead, read, len ) ;
			
			char *r ;
			int readInConsensusOffset = 0 ;
			int seqIdx ;
			int tag ;
			bool sortExtendedOverlaps = true ;
			struct _pair *oldMinExtAnchor = new struct _pair[ overlapCnt ] ;

			if ( overlaps[0].strand == 1 )
				r = read ;
			else
				r = rcRead ;
			
#ifdef DEBUG
			if ( overlaps[0].strand == -1 )
				printf( "rc: %s\n", r ) ;
#endif
			/*int tmpCnt = 0 ;
			for ( i = 0 ; i < overlapCnt ; ++i )
				if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].readStart == 0 && overlaps[i].readEnd == 74 && overlaps[i].similarity == 1 )
				{
					++tmpCnt ;
				}
			assert ( tmpCnt <= 1 ) ;*/
			// Filter low similiarty overlaps
			/*if ( similarityThreshold > novelSeqSimilarity )
			{
				// Filter extended sequence.
				int cnt = 0 ;
				for ( i = 0 ; i < overlapCnt ; ++i )
				{
					if ( overlaps[i].similarity >= similarityThreshold )
					{
						overlaps[ cnt ] = overlaps[i] ;
						++cnt ;
					}
				}
				overlapCnt = cnt ;
			}*/
	
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				oldMinExtAnchor[i].a = seqs[ overlaps[i].seqIdx ].minLeftExtAnchor ;
				oldMinExtAnchor[i].b = seqs[ overlaps[i].seqIdx ].minRightExtAnchor ;
				for ( j = 0 ; j < k ; ++j )
				{
					// If is contained in extended sequence.
					int leftRadius = radius ;
					int rightRadius = radius ;
					if ( extendedOverlaps[j].seqStart == 0 )
						leftRadius = 0 ;
					if ( extendedOverlaps[j].seqEnd == seqs[ extendedOverlaps[j].seqIdx ].consensusLen - 1 )
						rightRadius = 0 ;
					if ( overlaps[i].readStart >= extendedOverlaps[j].readStart - leftRadius  
						&& overlaps[i].readEnd <= extendedOverlaps[j].readEnd + rightRadius 
						&& ( overlaps[i].seqStart >= radius 
							|| overlaps[i].seqEnd <= seqs[ overlaps[i].seqIdx ].consensusLen - radius - 1 ) )
						break ;
					
					// Some extended is a subset of this one.
					leftRadius = radius ;
					rightRadius = radius ;
					if ( overlaps[i].seqStart == 0 )
						leftRadius = 0 ;
					if ( overlaps[i].seqEnd == seqs[ overlaps[i].seqIdx ].consensusLen - 1 )
						rightRadius = 0 ;
					if ( extendedOverlaps[j].readStart >= overlaps[i].readStart - leftRadius 
						&& extendedOverlaps[j].readEnd <= overlaps[i].readEnd + rightRadius )
						break ;
				}
				
				struct _seqWrapper &seq = seqs[ overlaps[i].seqIdx ] ; 
				if ( j < k || seq.isRef )
					continue ;
				//if ( !strcmp( read, "TTACTGTAATATACGATATTTTGACTGGTTATTAAGAGGCGACCCAAGAATCAATACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCT" ) && i == 1 )
				//	fprintf( stderr, "hi\n" ) ;
				// Only extend the novel seqs.
				if ( ExtendOverlap( r, len, seq, ( barcode == -1 && !repetitiveData ) ? 
						1.0 : 2.0, align, overlaps[i], extendedOverlaps[k] ) == 1 )
				{
					if ( extendedOverlaps[k].similarity < similarityThreshold )
					{
						if ( ( minKmerCount <= 1 || extendedOverlaps[k].similarity + 0.01 >= similarityThreshold )   
							&& extendedOverlaps[k].readStart == 0 
							&& extendedOverlaps[k].readEnd == len - 1 )
						{
							goodExtendedOverlap = extendedOverlaps[k] ;
						}
						continue ;
					}
					//if ( !strcmp( read, "CCGGCAGCCCCCAGGGAAGGGACTTGAATGGATTGGCTATATCTATTACACTGGGAGCACCATCTACAATCCCTC") )
					//{
					//	fprintf( stderr, "%d %d can be extended.\n", i, overlaps[i].seqIdx ) ;
					//}
					// Double check whether there is subset relationship.
					for ( j = 0 ; j < k ; ++j )
					{
						int leftRadius = radius ;
						int rightRadius = radius ;
						if ( extendedOverlaps[j].seqStart == 0 )
							leftRadius = 0 ;
						if ( extendedOverlaps[j].seqEnd == seqs[ extendedOverlaps[j].seqIdx ].consensusLen - 1 )
							rightRadius = 0 ;
						if ( extendedOverlaps[k].readStart >= extendedOverlaps[j].readStart - leftRadius  
								&& extendedOverlaps[k].readEnd <= extendedOverlaps[j].readEnd + rightRadius 
								&& ( overlaps[i].seqStart > 0 
									|| overlaps[i].seqEnd < seqs[ overlaps[i].seqIdx ].consensusLen - 1 ) )
							break ;

						// Some extended is a subset of this one.
						leftRadius = radius ;
						rightRadius = radius ;
						if ( extendedOverlaps[k].seqStart == 0 )
							leftRadius = 0 ;
						if ( extendedOverlaps[k].seqEnd == seqs[ extendedOverlaps[k].seqIdx ].consensusLen - 1 )
							rightRadius = 0 ;
						if ( extendedOverlaps[j].readStart >= extendedOverlaps[k].readStart - radius 
								&& extendedOverlaps[j].readEnd <= extendedOverlaps[k].readEnd + radius )
							break ;
					}
					if ( j < k )
						continue ;

					// Then check whether the extended porition is a subset of matched portion from other overlaps.
					for ( j = 0 ; j < i ; ++j )
					{
						if ( seqs[ overlaps[j].seqIdx ].isRef )
							continue ;

						if ( extendedOverlaps[k].seqStart == 0 
							&& extendedOverlaps[k].seqEnd == seqs[ extendedOverlaps[k].seqIdx ].consensusLen - 1 )
							continue ;

						if ( extendedOverlaps[k].readStart >= overlaps[j].readStart &&
							extendedOverlaps[k].readEnd <= overlaps[j].readEnd && 
							( overlaps[j].readEnd - overlaps[j].readStart >= 
								extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 10 
							  || overlaps[j].similarity + 0.02 >= extendedOverlaps[k].similarity ) )
						{
							if ( extendedOverlaps[k].readStart > 0 ) //&& extendedOverlaps[k].similarity >= 0.98 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							if ( extendedOverlaps[k].readEnd < len - 1 ) //&& extendedOverlaps[k].similarity >= 0.98 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							break ;
						}
					}
					if ( j < i )
						continue ;

					// The previous overlaps might just failed to extend at one side.
					for ( j = 0 ; j < failedExtendedOverlapsCnt ; ++j )
					{
						if ( extendedOverlaps[k].seqStart == 0 
							&& extendedOverlaps[k].seqEnd == seqs[ extendedOverlaps[k].seqIdx ].consensusLen - 1 )
							continue ;
						
						if ( extendedOverlaps[k].readStart >= failedExtendedOverlaps[j].readStart &&
							extendedOverlaps[k].readEnd <= failedExtendedOverlaps[j].readEnd ) 
							//failedExtendedOverlaps[j].similarity + 0.02 >= extendedOverlaps[k].similarity )
						{
							if ( extendedOverlaps[k].readStart > 0 ) //&& extendedOverlaps[k].similarity >= 0.98 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							if ( extendedOverlaps[k].readEnd < len - 1 ) //&& extendedOverlaps[k].similarity >= 0.98 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							break ;
						}
					}
					if ( j < failedExtendedOverlapsCnt  )
						continue ;
					
					if ( extendedOverlaps[k].readStart > 0 )
					{
						if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor >= extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
							continue ;
					}
					if ( extendedOverlaps[k].readEnd < len - 1 )
					{
						if ( seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor >= extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
							continue ;
					}


					tag = i ;
					++k ;
				}
				else
				{
					failedExtendedOverlaps[ failedExtendedOverlapsCnt ] = extendedOverlaps[k] ;
					++failedExtendedOverlapsCnt ;
				}
			}
				
			if ( k == 1 && 
				extendedOverlaps[0].readStart <= radius && extendedOverlaps[0].readEnd >= len - radius )
			{
				// Further check whehter this could merge two assemblies.
				// This happens when two contigs already overlap with each other.
				for ( i = 0 ; i < overlapCnt ; ++i )
				{
					if ( tag == i )
						continue ;

					struct _seqWrapper &seq = seqs[ overlaps[i].seqIdx ] ; 
					if ( seq.isRef )
						continue ;
				
					if ( ExtendOverlap( r, len, seq, ( barcode == -1 && !repetitiveData ) ? 1.0 : 2.0, 
						align, overlaps[i], extendedOverlaps[k] ) == 1 )
					{
						j = i ;
						++k ;
					}
				}

				if ( k > 2 )
				{
					k = 1 ;
				}
				else if ( k == 2 )
				{
					if ( extendedOverlaps[1].readStart > 0 )
					{
						if ( oldMinExtAnchor[j].a >= extendedOverlaps[1].readEnd - extendedOverlaps[1].readStart + 1 )
							k = 1 ;
					}
					if ( extendedOverlaps[1].readEnd < len - 1 )
					{
						if ( oldMinExtAnchor[j].b >= extendedOverlaps[1].readEnd - extendedOverlaps[1].readStart + 1 )
							k = 1 ;
					}
					
					if ( k == 2 )
					{
						if ( extendedOverlaps[0].seqEnd == seqs[ extendedOverlaps[0].seqIdx ].consensusLen - 1  
								&& extendedOverlaps[1].seqStart == 0  )
						{
							// no need to change.
							sortExtendedOverlaps = false ;
						}
						else if ( extendedOverlaps[0].seqStart == 0 && 
								extendedOverlaps[1].seqEnd == seqs[ extendedOverlaps[1].seqIdx ].consensusLen - 1 )
						{
							// swap 0, 1
							sortExtendedOverlaps = false ;

							struct _overlap tmp = extendedOverlaps[0] ; 
							extendedOverlaps[0] = extendedOverlaps[1] ;
							extendedOverlaps[1] = tmp ;
						}
						else
							k = 1 ;
					}
				}
			}

			if ( similarityThreshold > novelSeqSimilarity )
			{
				// Filter extended sequence.
				int cnt = 0 ;
				for ( i = 0 ; i < k ; ++i )
				{
					if ( extendedOverlaps[i].similarity >= similarityThreshold )
					{
						extendedOverlaps[ cnt ] = extendedOverlaps[i] ;
						++cnt ;
					}
				}
				k = cnt ;
			}

			if ( k == 0 && goodExtendedOverlap.seqIdx != -1 )
			{
				extendedOverlaps[0] = goodExtendedOverlap ;
				k = 1 ;
			}
			
			if ( k > 1 )
			{
				// For merging, if all the overlaps look bad, don't merge them
				for ( i = 0 ; i < k ; ++i )
				{
					if ( extendedOverlaps[i].similarity >= 0.95 )
						break ;
				}

				if ( i >= k )
				{
					int maxtag = 0 ;
					for ( i = 1 ; i < k ; ++i )
						if ( extendedOverlaps[i] < extendedOverlaps[ maxtag ] )
							maxtag = i ;
					extendedOverlaps[0] = extendedOverlaps[maxtag] ;
					k = 1 ;
				}
			}
#ifdef DEBUG
			for ( i = 0 ; i < k ; ++i )
				printf( "extended %d: %d %s. %d. %d %d %d %d %lf\n", i, extendedOverlaps[i].seqIdx, 
						seqs[ extendedOverlaps[i].seqIdx ].name, extendedOverlaps[i].strand, 
						extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, extendedOverlaps[i].seqStart, 
						extendedOverlaps[i].seqEnd, extendedOverlaps[i].similarity ) ; 
			fflush( stdout ) ;
#endif
			

			if ( k > 1 )
			{
				// If we are merging multiple novel seqs, make sure they are different seqs.
				for ( i = 0 ; i < k - 1 ; ++i )
					for ( j = i + 1 ; j < k ; ++j )
					{
						if ( extendedOverlaps[i].seqIdx == extendedOverlaps[j].seqIdx )
						{
							k = 0 ;
							break ;
						}
					}
			}

			/*if ( k == 1 ) 
			{
				// Check wether a match to the annotation fits better, in that case, creating 
				// 	a new sequence is a better.
				for ( i = 0 ; i < overlapCnt ; ++i )
				{
					if ( !seqs[ overlaps[i].seqIdx ].isRef )
						continue ;
					if ( overlaps[i].readStart == 0 && overlaps[i].readEnd == len - 1 
						&& overlaps[i].readStart <= extendedOverlaps[0].readStart 
						&& overlaps[i].readEnd >= extendedOverlaps[0].readEnd 
						&& overlaps[i].similarity > extendedOverlaps[0].similarity 
						&& overlaps[i].matchCnt > extendedOverlaps[0].matchCnt )
					{
						k = 0 ;
					}
				}
			}*/

			if ( k > 1 )
			{
				int eOverlapCnt = k ;
				addNew = false ;		
				// Merge sequences.
				// Reorder the overlaps to the order on the read coordinate.

				if ( sortExtendedOverlaps )
					std::sort( extendedOverlaps, extendedOverlaps + eOverlapCnt, CompSortOverlapsOnReadCoord ) ;
				
#ifdef DEBUG
				for ( i = 0 ; i < k ; ++i )
					printf( "sort extended %d: %d %s. %d. %d %d %d %d\n", i, extendedOverlaps[i].seqIdx, 
							seqs[ extendedOverlaps[i].seqIdx ].name, extendedOverlaps[i].strand, 
							extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, extendedOverlaps[i].seqStart, 
							extendedOverlaps[i].seqEnd ) ; 
#endif			
				// Infer from the gene name to make sure we would not merge two different receptors 
				for ( i = 0 ; i < k ; ++i )
				{
					for ( j = i + 1 ; j < k ; ++j )
					{
						if ( !IsNameCompatible( seqs[ extendedOverlaps[i].seqIdx ].name,
							seqs[ extendedOverlaps[j].seqIdx ].name ) )
						{
							delete[] extendedOverlaps ;
							delete[] failedExtendedOverlaps ;
							free( rcRead ) ;
							delete[] align ;
							delete[] oldMinExtAnchor ;
							return -1 ;
						}
					}
				}
				// Compute the new consensus.
				int sum = 0 ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					sum += seqs[ extendedOverlaps[i].seqIdx ].consensusLen ;
				char *newConsensus = (char *)malloc( sizeof(char) * ( sum + len + 1 ) ) ;
				
				// Compute the location of the seqs in the new merged seq
				int *seqOffset = new int[ eOverlapCnt ] ; 
				int base = 0 ;
				
				if ( extendedOverlaps[0].readStart > 0 )
				{
					for ( i = 0 ; i < eOverlapCnt ; ++i )
						seqOffset[i] = extendedOverlaps[i].readStart ;
				}
				else
				{
					seqOffset[0] = 0 ;
					
					for ( i = 1 ; i < eOverlapCnt ; ++i )
					{
						seqOffset[i] = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen - 1  
							+ ( extendedOverlaps[i].readStart - extendedOverlaps[i - 1].readEnd ) ;	
					}
				}
				
#ifdef DEBUG
				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					printf( "merge %d: %d %d %d %d %d. %d\n", i, extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, 
						extendedOverlaps[i].seqStart, extendedOverlaps[i].seqEnd, seqs[ extendedOverlaps[i].seqIdx ].consensusLen, 
						seqOffset[i] ) ;
				}
#endif
				if ( extendedOverlaps[0].readStart > 0 )
					memcpy( newConsensus, r, len ) ;
				else
					memcpy( newConsensus + extendedOverlaps[0].seqStart, r, len ) ;
				// Copy the original consensus in.
				// The earlier seq has higher weight.
				for ( i = eOverlapCnt - 1 ; i >= 0 ; --i )
				{
					memcpy( newConsensus + seqOffset[i], seqs[ extendedOverlaps[i].seqIdx ].consensus,
						seqs[ extendedOverlaps[i].seqIdx ].consensusLen ) ;	
				}

				int newConsensusLen = 0 ;
				int lastEndExtendedOverlapIdx = eOverlapCnt - 1 ;
				k = 0 ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					if ( seqOffset[i] + seqs[ extendedOverlaps[i].seqIdx ].consensusLen > k )
					{
						k = seqOffset[i] + seqs[ extendedOverlaps[i].seqIdx ].consensusLen ;
						lastEndExtendedOverlapIdx = i ;
					}

				if ( extendedOverlaps[lastEndExtendedOverlapIdx].readEnd < len )
				{
					newConsensusLen = k + ( len - extendedOverlaps[ lastEndExtendedOverlapIdx ].readEnd - 1 ) ;
				}
				else
				{
					newConsensusLen = k ; 
				}
				newConsensus[ newConsensusLen ] = '\0' ;
				//printf( "newconsensus %s %d\n", newConsensus, newConsensusLen ) ;
				
				// Rearrange the memory structure for posWeight.	
				int newSeqIdx = extendedOverlaps[0].seqIdx ;
				k = 0 ;
				for ( i = 1 ; i < eOverlapCnt ; ++i )
					if ( extendedOverlaps[i].seqIdx < newSeqIdx )
					{
						newSeqIdx = extendedOverlaps[i].seqIdx ;
						k = i ;
					}
				SimpleVector<struct _posWeight> &posWeight = seqs[newSeqIdx].posWeight ;
				posWeight.ShiftRight( seqOffset[k] ) ;
				posWeight.ExpandTo( newConsensusLen ) ;
				posWeight.SetZero( 0, seqOffset[k] ) ;
				posWeight.SetZero( seqOffset[k] + seqs[newSeqIdx].consensusLen, 
					newConsensusLen - seqs[ newSeqIdx ].consensusLen -seqOffset[k] ) ;
				

				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					// Though not the most efficient implementation, it seems very straightforward.
					int seqIdx = extendedOverlaps[i].seqIdx ;
					if ( seqIdx == newSeqIdx )
						continue ;

					for ( j = 0 ; j < seqs[ seqIdx ].consensusLen ; ++j )
					{
						int l ;
						posWeight[ seqOffset[i] + j ] += seqs[ seqIdx ].posWeight[j] ;
					}
				}

				// Update the index.
				KmerCode kmerCode( kmerLength ) ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					seqIndex.RemoveIndexFromRead( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus,
						seqs[ extendedOverlaps[i].seqIdx ].consensusLen, extendedOverlaps[i].seqIdx, 0 ) ;
				
				/*if ( seqOffset[0] != 0 )
				{
					seqIndex.UpdateIndexFromRead( kmerCode, seqs[ extendedOverlaps[0].seqIdx ].consensus, 
							seqs[ extendedOverlaps[0].seqIdx].consensusLen, seqOffset[0], 
							extendedOverlaps[0].seqIdx, newSeqIdx ) ; 
				}
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					//printf( "%d->%d.\n", extendedOverlaps[i].seqIdx, newSeqIdx ) ;
					// Don't use the overlapped portion, which will create duplicated entries in the index.
					int start = 0 ;
					if ( seqOffset[i] < seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen - kmerLength + 1 )
					{
						start = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen - kmerLength + 1 -
								seqOffset[i] ;
						printf( "start=%d\n", start ) ;
						seqIndex.RemoveIndexFromRead( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus, 
								start + kmerLength - 1, extendedOverlaps[i].seqIdx, 0 ) ;
					}
					seqIndex.UpdateIndexFromRead( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus, 
							seqs[ extendedOverlaps[i].seqIdx].consensusLen, seqOffset[i], 
							extendedOverlaps[i].seqIdx, newSeqIdx ) ;
				}
				
				// Update the index for the gaps.
				if ( extendedOverlaps[0].readStart > 0 )
				{
					seqIndex.BuildIndexFromRead( kmerCode, r, extendedOverlaps[0].readStart + kmerLength - 1, newSeqIdx ) ;
				}
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( extendedOverlaps[i].readStart > extendedOverlaps[i - 1].readEnd - kmerLength + 1 )
					{
						int start = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen 
									- kmerLength + 1 ;
						int end = seqOffset[i] + kmerLength - 2 ;
						int rstart = extendedOverlaps[i - 1].readEnd - kmerLength + 2 ;
						seqIndex.BuildIndexFromRead( kmerCode, r + rstart, end - start + 1, newSeqIdx, start ) ;
					}
				}
				if ( extendedOverlaps[i - 1].readEnd < len - 1 )
				{
					int start = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen 
						- kmerLength + 1 ;
					int rstart = extendedOverlaps[i - 1].readEnd - kmerLength + 2 ;
					seqIndex.BuildIndexFromRead( kmerCode, r + rstart, len - rstart, newSeqIdx, start ) ;
				}
				*/

				// Update the name.
				// TODO: use array of names.
				sum = 0 ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					sum += strlen( seqs[ extendedOverlaps[i].seqIdx ].name ) ;
				char* nameBuffer = new char[sum + eOverlapCnt + 1 ] ;
				
				strcpy( nameBuffer, seqs[ extendedOverlaps[0].seqIdx ].name ) ;
				sum = strlen( nameBuffer ) ;
				//printf( "%d\n", seqs.size() ) ;
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( strcmp( seqs[ extendedOverlaps[i].seqIdx ].name, seqs[ extendedOverlaps[i - 1].seqIdx ].name ) )
					{
						nameBuffer[ sum ] = '+' ;
						nameBuffer[ sum + 1 ] = '\0' ;
						strcpy( nameBuffer + sum + 1, seqs[ extendedOverlaps[i].seqIdx ].name ) ;
						sum = sum + 1 + strlen( seqs[ extendedOverlaps[i].seqIdx ].name ) ;
					}
				}
				free( seqs[ newSeqIdx ].name ) ;
				seqs[ newSeqIdx ].name = strdup( nameBuffer ) ;
				delete[] nameBuffer ;

				// Relase the memory for merged seqs.
				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					int seqIdx = extendedOverlaps[i].seqIdx ;
					if ( seqIdx == newSeqIdx )
						continue ;
					
					ReleaseSeq( seqIdx ) ;
				}
					
				// Put the new allocated stuff in.
				free( seqs[ newSeqIdx ].consensus ) ;
				seqs[ newSeqIdx ].consensus = newConsensus ;
				seqs[ newSeqIdx ].consensusLen = newConsensusLen ;
				
				// Update the index
				UpdateConsensus( newSeqIdx, false ) ;
				seqIndex.BuildIndexFromRead( kmerCode, newConsensus, newConsensusLen, newSeqIdx ) ;
				
				// Update the anchor requirement.
				seqs[ newSeqIdx ].minLeftExtAnchor = seqs[ extendedOverlaps[0].seqIdx ].minLeftExtAnchor ;
				seqs[ newSeqIdx ].minRightExtAnchor = 
					seqs[ extendedOverlaps[ lastEndExtendedOverlapIdx ].seqIdx ].minRightExtAnchor ;
				
				// either one of the ends of read or seq should be 0.
				readInConsensusOffset = 0 ;
				if ( extendedOverlaps[0].seqStart > 0 )
					readInConsensusOffset = extendedOverlaps[0].seqStart ;

				delete[] seqOffset ;

				seqIdx = newSeqIdx ;
			}
			else if ( k == 1 )
			{
				// Extend a sequence
				addNew = false ;

				seqIdx = extendedOverlaps[0].seqIdx ;
				struct _seqWrapper &seq = seqs[ extendedOverlaps[0].seqIdx ] ;
				
				// Compute the new consensus.
				if ( extendedOverlaps[0].readStart > 0 || extendedOverlaps[0].readEnd < len - 1 )
				{
					char *newConsensus = (char *)malloc( sizeof( char ) * ( 
						( extendedOverlaps[0].readStart + len - 1 -extendedOverlaps[0].readEnd ) + seq.consensusLen + 1 ) ) ;

					if ( extendedOverlaps[0].readStart > 0 )
					{
						// add read[0...readStart-1] to the consensus.
						//for ( i = 0 ; i < extendedOverlaps[0].readStart ; ++i )
						//	newConsensus[i] = r[i] ;
						memcpy( newConsensus, r, extendedOverlaps[0].readStart ) ;
					}
					memcpy( newConsensus + extendedOverlaps[0].readStart, seq.consensus, seq.consensusLen ) ;
					j = extendedOverlaps[0].readStart + seq.consensusLen ;	

					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						for ( i = extendedOverlaps[0].readEnd + 1 ; i < len ; ++i, ++j )
							newConsensus[j] = r[i] ;
					}
					newConsensus[j] = '\0' ;
					int newConsensusLen = strlen( newConsensus ) ;	
					//printf( "new consensus %s %d. %s %d\n", seq.consensus, seq.consensusLen, newConsensus, j ) ;
					
					// Update index 
					int shift = extendedOverlaps[0].readStart ;
					KmerCode kmerCode( kmerLength ) ;
					if ( shift > 0 )
					{
						seqIndex.BuildIndexFromRead( kmerCode, newConsensus, extendedOverlaps[0].readStart + kmerLength - 1, seqIdx ) ;
						seqIndex.UpdateIndexFromRead( kmerCode, seq.consensus, seq.consensusLen, shift, seqIdx, seqIdx ) ; 
					}
					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						int start = extendedOverlaps[0].readStart + extendedOverlaps[0].seqEnd - kmerLength + 2 ;
						seqIndex.BuildIndexFromRead( kmerCode, newConsensus + start , 
							( newConsensusLen - start ), seqIdx, 
							start ) ;
					}
					
					// Rearrange the memory structure for posWeight.
					int expandSize = extendedOverlaps[0].readStart + ( len - 1 - extendedOverlaps[0].readEnd ) ;
					seq.posWeight.ExpandBy( expandSize ) ;
					if ( shift > 0 )
					{
						for ( i = seq.consensusLen - 1 ; i >= 0 ; --i )
							seq.posWeight[i + shift] = seq.posWeight[i] ;
						
						// Lower the weight at the end for original sequence.
						for ( i = 0 ; i < 2 ; ++i )
						{
							if ( i + shift >= len || r[i + shift] == 'N' )
								continue ;
							for ( j = 0 ; j < 4 ; ++j )
								if ( r[i + shift] != numToNuc[j] && seq.posWeight[i + shift].count[j] > 1 )
									--seq.posWeight[i + shift].count[j] ;
						}
						seq.posWeight.SetZero( 0, shift ) ;
					}
					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						int start = extendedOverlaps[0].readStart + seq.consensusLen ;
						seq.posWeight.SetZero( start, len - extendedOverlaps[0].readEnd - 1 ) ;
						
						// Lower the weight at the end for original sequence.
						for ( i = seq.consensusLen - 2 ; i < seq.consensusLen ; ++i )
						{
							int pos = i - extendedOverlaps[0].seqStart ;
							if ( pos < 0 || r[pos] == 'N' )
								continue ;
							for ( j = 0 ; j < 4 ; ++j )
								if ( r[pos] != numToNuc[j] && seq.posWeight[i].count[j] > 1 )
									--seq.posWeight[i].count[j] ;
						}
					}
					
					// Update the anchor requirement
					if ( shift > 0 )
						seq.minLeftExtAnchor = 0 ;
					if ( extendedOverlaps[0].readEnd < len - 1 )
						seq.minRightExtAnchor = 0 ;

					// Adjust the name.
					// Find the possible ref seq
					int refIdx = -1 ;
					for ( i = 0 ; i < overlapCnt ; ++i )
					{
						if ( !seqs[ overlaps[i].seqIdx ].isRef )
							continue ;
						// Use refIdx as the idx in the overlaps list.
						if ( refIdx == -1 || 
							( overlaps[i].readEnd - overlaps[i].readStart > 
								overlaps[refIdx].readEnd - overlaps[refIdx].readStart ) )
						{
							refIdx = i ;
						}

						if ( strstr( seq.name, seqs[ overlaps[i].seqIdx ].name ) != NULL )
						{
							refIdx = i ;
							break ;
						}
					}
					if ( refIdx != -1 )
					{
						// Use refIdx as the idx in the seqs 
						refIdx = overlaps[ refIdx ].seqIdx ;
						if ( strstr( seq.name, seqs[ refIdx ].name ) == NULL )
						{
							char *nameBuffer = new char[ strlen( seqs[ refIdx ].name) + strlen( seq.name ) + 2 ] ;
							if ( extendedOverlaps[0].readStart > 0 )
							{
								sprintf( nameBuffer, "%s+%s", seqs[refIdx].name, seq.name ) ;
							}
							else
							{
								sprintf( nameBuffer, "%s+%s", seq.name, seqs[refIdx].name ) ;
							}
							free( seq.name ) ;
							seq.name = strdup( nameBuffer ) ;
							delete[] nameBuffer ;
						}
					}
					// either one of the ends of read or seq should be 0.
					readInConsensusOffset = 0 ;
					if ( extendedOverlaps[0].seqStart > 0 )
						readInConsensusOffset = extendedOverlaps[0].seqStart ;

					free( seq.consensus ) ;
					seq.consensus = newConsensus ;
					seq.consensusLen = newConsensusLen ;	
					//printf( "new consensus len %d\n", seq.consensusLen ) ;
				}
				else // the read is inside of the seq.
					readInConsensusOffset = extendedOverlaps[0].seqStart ;
			}

			// Update the posweight, assume we already compute the new readStart and shift existing posWeight.
			// seqIdx holds the index that we need to update.
			if ( !addNew )
			{
				struct _seqWrapper &seq = seqs[seqIdx] ;
				//printf( "%d %d. %d %d\n%s\n%s\n", seqIdx, seq.posWeight.Size(), readInConsensusOffset, len, seq.consensus, r) ;
				SimpleVector<int> nPos ;
				for ( i = 0 ; i < len ; ++i )
				{
					if ( r[i] == 'N' )
						continue ;
					++seq.posWeight[i + readInConsensusOffset].count[ nucToNum[ r[i] - 'A' ] ] ;
			
					if ( seq.consensus[i + readInConsensusOffset ] == 'N' )
					{
						nPos.PushBack( i ) ;
					}
				}
				SetPrevAddInfo( seqIdx, 0, len - 1, readInConsensusOffset, readInConsensusOffset + len - 1, overlaps[0].strand ) ;

				int size = nPos.Size() ;
				for ( i = 0 ; i < size ;  )
				{
					for ( j = i + 1 ; j < size ; ++j )
						if ( nPos[j] > nPos[j - 1] + kmerLength - 1 )
							break ;

					// [i,j) holds the N positions that are with kmer-length size.
					int l ;
					// Update the consensus
					for ( l = i ; l < j ; ++l )
						seq.consensus[ nPos[l] + readInConsensusOffset ] = r[ nPos[l] ] ;
					// Update the index
					KmerCode kmerCode( kmerLength ) ;
					int start = nPos[i] - kmerLength + 1 + readInConsensusOffset ;
					if ( start < 0 )
						start = 0 ;
					int end = nPos[j - 1] + kmerLength - 1 + readInConsensusOffset ;
					if ( end >= seq.consensusLen )
						end = seq.consensusLen - 1 ;
					seqIndex.BuildIndexFromRead( kmerCode, seq.consensus + start, end - start + 1, seqIdx, start  ) ;
					i = j ;
				}

				ret = seqIdx ;
			}
			
			free( rcRead ) ;
			delete[] align ;
			delete[] oldMinExtAnchor ;
		}

		k = 0 ;
		// See whether there is a reference seq match is sequence if it does not match any novel seq.
		int anchorSeqIdx = -1 ;
		if ( addNew )
		{
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				if ( seqs[ overlaps[i].seqIdx ].isRef )
				{
					anchorSeqIdx = overlaps[i].seqIdx ;
					break ;
				}
			}
			if ( i >= overlapCnt )
				addNew = false ;

			/*if ( !addNew )
			{
				for ( i = 0 ; i < overlapCnt ; ++i )
					if ( overlaps[i].similarity == 1.0 )
					{
						addNew = true ;
						refSeqIdx = -1 ;
						break ;
					}
			}*/
			/*if ( anchorSeqIdx == -1 )
			{
				for ( i = 0 ; i < overlapCnt ; ++i )
				{
					if ( overlaps[i].matchCnt >= 2 * hitLenRequired 
							&& ( overlaps[i].readStart < 3 
								|| overlaps[i].readEnd >= len - 3 ) )
					{
						anchorSeqIdx = i ;
						break ;
					}
				}
			}*/
		}

		if ( addNew )
		{
			// Add the sequence to SeqSet
			int idx = seqs.size() ;
			struct _seqWrapper ns ;
			
			if ( anchorSeqIdx >=0 )
				ns.name = strdup( seqs[ anchorSeqIdx ].name ) ;
			else
				ns.name = strdup( "unknown" ) ;
			ns.consensus = strdup( read ) ;
			ns.consensusLen = strlen( read ) ;
			if ( overlaps[ anchorSeqIdx ].strand == -1 )
				ReverseComplement( ns.consensus, read, len ) ;
			ns.isRef = false ;
			
			ns.posWeight.Reserve( len ) ;
			ns.posWeight.ExpandTo( len ) ;
			ns.posWeight.SetZero( 0, len ) ;
			for ( i = 0 ; i < len ; ++i )
			{
				//memset( ns.posWeight[i].count, 0, sizeof( ns.posWeight[i].count ) ) ;
				if ( ns.consensus[i] == 'N' )
					continue ;
				
				++ns.posWeight[i].count[ nucToNum[ ns.consensus[i] - 'A' ] ] ;
			}
			//printf( "%d %s %lld\n", ns.posWeight.Size(), ns.consensus, ns.posWeight.BeginAddress() ) ;
			ns.minLeftExtAnchor = ns.minRightExtAnchor = 0 ;
			seqs.push_back( ns ) ;

			// Don't forget to update index.
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, len, idx ) ;			
			
			SetPrevAddInfo( idx, 0, len - 1, 0, len - 1, overlaps[0].strand ) ; 
#ifdef DEBUG
			printf( "add novel seq: %d\n", idx ) ;
#endif
			ret = idx ;
		}

		delete[] extendedOverlaps ;
		delete[] failedExtendedOverlaps ;

		if ( ret == -1 )
		{
			SetPrevAddInfo( -2, -1, -1, -1, -1, 0 ) ; 
			ret = -2 ;
		}

		if ( ret >= 0 && strand == 0 )
			strand = overlaps[0].strand ;

		return ret ;
	}
	
	// Called when we just want to duplicate the add operation 
	//   applied before.
	int RepeatAddRead( char *read )
	{
		if ( prevAddInfo.seqIdx < 0 )
			return prevAddInfo.seqIdx ;

		int i ;
		char *r ;
				
		r = read ;
		if ( prevAddInfo.strand == -1 )
		{
			int len = strlen( read ) ;
			r = strdup( read ) ;
			ReverseComplement( r, read, len ) ;
		}
		
		struct _seqWrapper &seq = seqs[ prevAddInfo.seqIdx ] ;
		//printf( "%d %d. %d %d\n%s\n%s\n", seqIdx, seq.posWeight.Size(), readInConsensusOffset, len, seq.consensus, r) ;
		for ( i = prevAddInfo.readStart ; i <= prevAddInfo.readEnd ; ++i )
		{
			if ( r[i] == 'N' )
				continue ;
			++seq.posWeight[i + prevAddInfo.seqStart].count[ nucToNum[ r[i] - 'A' ] ] ;
		}

		if ( prevAddInfo.strand == -1 )
			free( r ) ;

		return prevAddInfo.seqIdx ;
	}

	void AddAssignedRead( char *read, struct _overlap assign )
	{
		if ( assign.seqIdx == -1 )
			return ;

		char *r = strdup( read ) ;
		int len = strlen( read ) ;
		if ( assign.strand == -1 )
			ReverseComplement( r, read, len ) ;
		
		UpdatePosWeightFromRead( seqs[ assign.seqIdx ].posWeight, assign.seqStart, r ) ;
		free( r ) ;
	}
	

	void UpdateAllConsensus()
	{
		int i ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;
			UpdateConsensus( i, true ) ;
		}
	}

	void UpdateConsensus( int seqIdx, bool updateIndex )
	{
		int i, j ;
		struct _seqWrapper &seq = seqs[ seqIdx ] ;
		SimpleVector<struct _pair> changes ;
		for ( i = 0 ; i < seq.consensusLen ; ++i )
		{
			int max = 0 ;
			int maxTag = 0 ;
			for ( j = 0 ; j < 4 ; ++j )
				if ( seq.posWeight[i].count[j] > max )
				{
					max = seq.posWeight[i].count[j] ;
					maxTag = j ;
				}

			if ( max == 0 ) // A case of N
				continue ;

			if ( nucToNum[ seq.consensus[i] - 'A' ] != maxTag 
				&& seq.posWeight[i].count[  nucToNum[ seq.consensus[i] - 'A' ] ] < max )
			{
				struct _pair np ;
				np.a = i ;
				np.b = maxTag ;
				changes.PushBack( np ) ;
			}
		}

		if ( changes.Size() == 0 )
			return ;
		
		if ( updateIndex )
		{
			// Inefficient implementation, improve in future.
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.RemoveIndexFromRead( kmerCode, seq.consensus, seq.consensusLen, seqIdx, 0 ) ;
		}

		int size = changes.Size() ;
		for ( i = 0 ; i < size ; ++i )
			seq.consensus[ changes[i].a ] = numToNuc[ changes[i].b ] ;

		if ( updateIndex )
		{
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.BuildIndexFromRead( kmerCode, seq.consensus, seq.consensusLen, seqIdx, 0 ) ;
		}
	}
	
	// Remove unneeded entries and rebuild the index.
	void Clean( bool removeRefSeq )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		k = 0 ;
		KmerCode kmerCode( kmerLength ) ;
		seqIndex.Clear() ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL )
			{
				continue ;
			}
			if ( removeRefSeq && seqs[i].isRef )
			{
				//seqIndex.RemoveIndexFromRead( kmerCode, seqs[i].consensus, seqs[i].consensusLen, i, 0 ) ;
				
				ReleaseSeq( i ) ;
				continue ;
			}

			seqs[k] = seqs[i] ;
			//if ( i != k )
			//	seqIndex.UpdateIndexFromRead( kmerCode, seqs[k].consensus, seqs[k].consensusLen, 0, i, k ) ;
			seqIndex.BuildIndexFromRead( kmerCode, seqs[k].consensus, seqs[k].consensusLen, k, 0 ) ;
			++k ;
		}
		SetPrevAddInfo( -1, -1, -1, -1, -1, 0 ) ; 
		seqs.resize( k ) ;
	}

	void ChangeKmerLength( int kl )
	{
		kmerLength = kl ;
		Clean( false ) ;
	}
	
	// Find the seq id this read belongs to.
	int AssignRead( char *read, int strand, int barcode, struct _overlap &assign )
	{
		int i ;

		std::vector<struct _overlap> overlaps ;
		
		int overlapCnt = GetOverlapsFromRead( read, strand, barcode, 0, false, overlaps ) ;
		//printf( "%d %d\n", overlapCnt, mateOverlapCnt ) ;
		//printf( "%d %s\n%d %s\n", overlaps[0].strand, reads[i].seq, mateOverlaps[0].strand, reads[i + 1].seq ) ;
		assign.seqIdx = -1 ;

		if ( overlapCnt == 0 )
		{
			return -1 ;
		}
			
		std::sort( overlaps.begin(), overlaps.end() ) ;

		int len = strlen( read ) ;
		char *rc = new char[len + 1] ;

		ReverseComplement( rc, read, len ) ;

		char *r = read ;
		if ( overlaps[0].strand == -1 )
			r = rc ;

		struct _overlap extendedOverlap ;
		char *align = new char[ 2 * len + 2 ] ;
		/*for ( j = 0 ; j < overlapCnt ; ++j )
		  {
		  printf( "+ %d %d: %d %d %d %lf\n", i, j, overlaps[j].seqIdx, overlaps[j].seqStart, overlaps[j].seqEnd, overlaps[j].similarity) ;
		  }
		  for ( j = 0 ; j < mateOverlapCnt ; ++j )
		  {
		  printf( "- %d %d: %d %d %d %lf\n", i + 1, j, mateOverlaps[j].seqIdx, mateOverlaps[j].seqStart, mateOverlaps[j].seqEnd, mateOverlaps[j].similarity) ;
		  }*/
		int extendCnt = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			//printf( "%d %d: %d-%d %d-%d %lf\n", i, overlaps[i].seqIdx, overlaps[i].readStart, overlaps[i].readEnd,
			//		overlaps[i].seqStart, overlaps[i].seqEnd, overlaps[i].similarity) ;
			if ( ExtendOverlap( r, len, seqs[ overlaps[i].seqIdx ], barcode == -1 ? 1.0 : 2.0, align, 
						overlaps[i], extendedOverlap ) == 1 )
			{
				/*if ( extendCnt == 0 )
				  {
				  extendedOverlap = tmpExtendedOverlap ;
				  ++extendCnt ;
				  }
				  else if ( tmpExtendedOverlap.similarity == extend) */
				if ( extendedOverlap.readStart == 0 && extendedOverlap.readEnd == len - 1 )
					break ;
			}
		}

		if ( i >= overlapCnt )
		{
			delete[] rc ;
			delete[] align ;
			return -1 ;
		}
	
		delete[] rc ;
		delete[] align ;
		assign = extendedOverlap ;
		return assign.seqIdx ;
	}
	

	// Recompute the posweight based on assigned read
	void RecomputePosWeight( std::vector<struct _assignRead> &reads )
	{
		int i, j ;
		int readCnt = reads.size() ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
			seqs[i].posWeight.SetZero( 0, seqs[i].consensusLen ) ;

		for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( reads[i].overlap.seqIdx == -1 )
				continue ;
			
			if ( reads[i].overlap.strand == 1 )
				UpdatePosWeightFromRead( seqs[ reads[i].overlap.seqIdx ].posWeight, reads[i].overlap.seqStart, reads[i].read ) ;
			else
			{
				char *r = strdup( reads[i].read ) ;
				ReverseComplement( r, reads[i].read, strlen( r ) ) ;
				UpdatePosWeightFromRead( seqs[ reads[i].overlap.seqIdx ].posWeight, reads[i].overlap.seqStart, r ) ;
				free( r ) ;
			}
		}

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			struct _seqWrapper &seq = seqs[i] ;
			for ( j = 0 ; j < seq.consensusLen ; ++j )
				if ( seq.consensus[j] != 'N' && seq.posWeight[j].Sum() == 0 )
				{
					seq.posWeight[j].count[ nucToNum[ seq.consensus[j] - 'A' ] ] = 1 ;
				}
		}
	}
	
	// Return:the number of connections made.
	int ExtendSeqFromSeqOverlap( int leastOverlapLen )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int ret = 0 ;
		repeatSimilarity = 0.99 ;

		std::vector<struct _overlap> *adj = new std::vector<struct _overlap>[ seqCnt ] ;
		struct _pair *next = new struct _pair[ seqCnt ] ; // a-index, b-the index in adj[i] 
		struct _pair *prev = new struct _pair[ seqCnt ] ;
		int *containedIn = new int[seqCnt] ;
	
		//BuildSeqOverlapGraph( 100, adj ) ;
		SimpleVector<bool> useInBranch ;
		useInBranch.ExpandTo( seqCnt ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
			useInBranch[i] = true ;
		
		BuildBranchGraph( adj, leastOverlapLen ) ;

		// Keep only the connections that representing overlapping information.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( adj[i][j].readEnd < seqs[i].consensusLen - 1 ||
					adj[i][j].seqStart > 0 )
					adj[i][j].seqIdx = -1 ;
				else if ( adj[i][j].readStart == 0 && adj[i][j].readEnd == seqs[i].consensusLen - 1
					&& adj[i][j].seqStart == 0 && adj[i][j].seqEnd == seqs[ adj[i][j].seqIdx ].consensusLen - 1 
					&& i > adj[i][j].seqIdx )
					adj[i][j].seqIdx = -1 ;
			}
			
			k = 0 ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( adj[i][j].seqIdx != -1 )
				{
					adj[i][k] = adj[i][j] ;
					++k ;
				}
			}
			adj[i].resize( k ) ;
		}

		memset( next, -1, sizeof( struct _pair ) * seqCnt ) ;
		memset( prev, -1, sizeof( struct _pair ) * seqCnt ) ;
		memset( containedIn, -1, sizeof( int ) * seqCnt ) ;
		KmerCode kmerCode( kmerLength ) ;

		// Process the contained in.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			if ( size == 0 )
				continue ;
			
			if ( adj[i][0].readStart == 0 && adj[i][0].readEnd == seqs[i].consensusLen - 1 &&
				containedIn[ adj[i][0].seqIdx ] == -1 )
			{
				int seqIdx = adj[i][0].seqIdx ;
				//printf( "Contain: %d %d %lf\n%s\n%s\n", adj[i][0].readStart, adj[i][0].readEnd, adj[i][0].similarity, 
				//			seqs[i].consensus, seqs[ seqIdx ].consensus ) ;
				
				for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
				{
					seqs[ seqIdx ].posWeight[ adj[i][0].seqStart + j ] += seqs[i].posWeight[j] ;
				}

				seqIndex.RemoveIndexFromRead( kmerCode, seqs[i].consensus, seqs[i].consensusLen, i, 0 ) ;	
				ReleaseSeq( i ) ;

				containedIn[i] = seqIdx ;
				next[i].a = -2 ;
				prev[i].a = -2 ;
				
				++ret ;
			}
		}

		//Clean( true ) ;
		//return ret ;
		
		// Process the partial overlap case
		// Build the path
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			if ( size == 0 )
				continue ;
			k = 0 ;
			if ( containedIn[i] != -1 )
				continue ;

			if ( containedIn[ adj[i][0].seqIdx ] != -1 )
			{
				int father = containedIn[ adj[i][0].seqIdx ] ;
				while ( 1 )
				{
					if ( containedIn[ father ] == -1 )
						break ;
					father = containedIn[ father ] ;
				}
				for ( j = 0 ; j < size ; ++j )
					if ( adj[i][j].seqIdx == father )
						break ;
				if ( j < size )
					k = j ;
				else
					continue ;
			}
			if ( prev[ adj[i][k].seqIdx ].a == -1 )
			{
				next[i].a = adj[i][k].seqIdx ;
				next[i].b = k ;
				prev[ adj[i][k].seqIdx ].a = i ;
				prev[ adj[i][k].seqIdx ].b = k ;
			}
			else if ( prev[ adj[i][k].seqIdx ].a >= 0 )
			{
				next[ prev[ adj[i][k].seqIdx ].a ].a = -2 ;
				prev[ adj[i][k].seqIdx ].a = -2 ;
			}
		}
		
		//for ( i = 0 ; i < seqCnt ; ++i )
		//	printf( "%d next %d prev %d contained %d\n", i, next[i].a, prev[i].a, containedIn[i] ) ;
		
		// Use the paths to merge seqs.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( containedIn[i] == -1 && prev[i].a < 0 && next[i].a >= 0 ) // Head of a chain
			{
				std::vector<int> path ;
				int p = i ;
				int newConsensusLen = 0 ;
				std::vector<int> seqOffset ; 
				int offset = 0 ;

				k = 0 ;
				while ( 1 )
				{
					path.push_back( p ) ;
					++k ;
					seqOffset.push_back( offset ) ;
					
					if ( next[p].a >= 0 )
						offset += adj[p][ next[p].b ].readStart ;
					else
						break ;
					p = next[p].a ;
				}
			
				struct _seqWrapper ns ;
				ns.isRef = false ;

				// Obtain the length after merging.
				newConsensusLen = seqOffset[k - 1] + seqs[ path[k - 1] ].consensusLen ;
				char *newConsensus = (char *)malloc( sizeof( char ) * newConsensusLen + 1 ) ;
				for ( j = 0 ; j < k ; ++j )
					memcpy( newConsensus + seqOffset[j], seqs[ path[j] ].consensus, seqs[ path[j] ].consensusLen ) ;
				newConsensus[ newConsensusLen ] = '\0' ;
				
				ns.consensus = newConsensus ;
				ns.consensusLen = newConsensusLen ;
				ns.posWeight.ExpandTo( newConsensusLen ) ;
				
				// Update posweight
				ns.posWeight.SetZero( 0, newConsensusLen ) ;
				for ( j = 0 ; j < k ; ++j )
				{
					int l, c ;
					int seqIdx = path[j] ;
					for ( l = 0 ; l < seqs[ seqIdx ].consensusLen ; ++l )
						ns.posWeight[l + seqOffset[j] ] += seqs[ seqIdx ].posWeight[l] ; 
				} 
				// Update name
				int sum = 0 ;
				for ( j = 0 ; j < k ; ++j )
					sum += strlen( seqs[ path[j] ].name ) ;
				sum += k ;
				
				ns.name = (char *)malloc( sizeof( char ) * sum ) ;
				strcpy( ns.name, seqs[ path[0] ].name ) ;
				sum = strlen( ns.name ) ;
				for ( j = 1 ; j < k ; ++j )
				{
					ns.name[ sum ] = '+' ;
					ns.name[ sum + 1 ] = '\0' ;
					strcpy( ns.name + sum + 1, seqs[ path[j] ].name ) ;
					sum = sum + 1 + strlen( seqs[ path[j] ].name ) ;
				}

				// Update index	
				//int newSeqIdx = seqs.size() ;
				seqs.push_back( ns ) ;
				for ( j = 0 ; j < k ; ++j )
				{
					//seqIndex.RemoveIndexFromRead( kmerCode, seqs[ path[j] ].consensus, seqs[ path[j] ].consensusLen, path[j], 0 ) ;						
					ReleaseSeq( path[j] ) ;
				}
				//seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, ns.consensusLen, newSeqIdx ) ;
			
				ret += k - 1 ;
			}
		}


		Clean( true ) ;
		
		delete[] adj ;
		delete[] next ;
		delete[] prev ;
		delete[] containedIn ;
		return ret ;
	}

	// Remove the seq that are substring of other seqs.
	int RemoveRedundantSeq()
	{
		int i ;
		int seqCnt = seqs.size() ;
		std::vector<struct _overlap> subsetOf ;
		subsetOf.resize( seqCnt ) ;
		
		for ( i = 0 ; i < seqCnt ; ++i )
			subsetOf[i].seqIdx = -1 ;	
		
		BuildSeqSubstringRelation( subsetOf ) ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( subsetOf[i].seqIdx != -1 )
				ReleaseSeq( i ) ;
		}
		Clean( true ) ;
	
		return seqs.size() ;
	}

	/*int AddRead( char *read )
	{
		int i, j, k ;
		int len = strlen( read ) ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;
		
		overlapCnt = GetOverlapsFromRead( read, overlaps ) ;
		
		if ( overlapCnt == 0 )
			return -1 ;
		
		struct _overlap *extendedOverlaps = new struct _overlap[ overlapCnt ] ;
		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			std::sort( overlaps.begin(), overlaps.end() ) ;
			if ( Extend)
		}
		int eOverlapCnt ;

		delete[] overlapCnt ;
	}*/

	void ResetPosWeight()
	{
		int i ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = seqs[i].posWeight.Size() ;
			seqs[i].posWeight.SetZero( 0, size ) ;
		}
	}
	
	int GetSeqCnt() 
	{
		return seqs.size() ;
	}

	// Return: 2-bit: left bit: V motif, right bit: J motif
	int HasMotif( char *read, int strand )
	{
		if ( strand == 0 )
			return 0 ;

		int i, j, k ;
		char *r = read ;
		char *aa = strdup( read ) ;
		int len = strlen( read ) ;
		if ( strand == -1 )
		{
			r = strdup( read ) ;
			ReverseComplement( r, read, len ) ;
		}
		
		int ret = 0 ;
		for ( k = 0 ; k <= 2 ; ++k )	
		{
			for ( i = k, j = 0 ; i < len ; i += 3, ++j )
				aa[j] = DnaToAa( read[i], read[i + 1], read[i + 2] ) ;

			for ( i = 0 ; i + 2 < j ; ++i ) // YYC
				if ( aa[i] == 'Y' && aa[i + 1] == 'Y' && aa[i + 2] == 'C' )
				{
					ret |= 2 ;
					break ;
				}

			for ( i = 0 ; i + 3 < j ; ++i ) // F/W G*G
				if ( ( aa[i] == 'F' || aa[i] == 'W' ) && 
					aa[i + 1] == 'G' && aa[i + 3] == 'G' )
				{
					ret |= 1 ;
					break ;
				}

		}

		if ( strand == -1 )
		{
			free( r ) ;	
		}
		free( aa ) ;

		return ret ;
	}

	int GetGeneType( char *name )
	{
		int geneType = -1 ;
		switch ( name[3] )
		{
			case 'V': geneType = 0 ; break ;
			case 'D': 
				  if ( name[4] >= '0' && name[4] <= '9' )
					  geneType = 1 ; 
				  else
					  geneType = 3 ;
				  break ;
			case 'J': geneType = 2 ; break ;
			default: geneType = 3 ; break ;
		}
		return geneType ;
	}

	bool IsSameGeneAllele( char *name, char *name2 )
	{
		int i ;
		int ret = true ;
		for ( i = 0 ; name[i] && name2[i] && name[i] != '*' && name2[i] != '*' ; ++i )
		{
			if ( name[i] != name2[i] )
			{
				ret = false ;
				break ;
			}
		}

		return ret ;
	}

	bool IsSameChainType( char *name, char *name2 )
	{
		int i ;
		for ( i = 0 ; name[i] && name2[i] && i < 3 ; ++i )
		{
			if ( name[i] != name2[i] )
				return false ;
		}
		if ( i >= 3 )
			return true ;
		else
			return false ;
	}
	
	// The reference gene may have different length, which makes matchCnt criterion biased
	//   to longer gene, so we want to remove such effect
	// Return: is overlap a is better than b*threshold
	bool IsBetterGeneMatch( struct _overlap &a, struct _overlap &b, double threshold )
	{
		int matchCnt = a.matchCnt ;
		/*if ( a.seqEnd - a.seqStart + 1 > 0.95 * seqs[ a.seqIdx ].consensusLen ) 
			matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
		else if ( GetGeneType( seqs[ a.seqIdx ].name ) == 2 && a.seqEnd >= seqs[ a.seqIdx ].consensusLen - 3 
			&& a.readEnd >= seqs[a.seqIdx].consensusLen )
		{
			matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
		}*/
		int gapAllow = kmerLength + 1 ;
		if ( threshold >= 1 )
			gapAllow = 3 ;
		if ( a.seqIdx == -1 )
			return false ;
		if ( b.seqIdx == -1 )
			return true ;
		int geneType = GetGeneType( seqs[ a.seqIdx ].name ) ;
		if ( geneType == 2 /*&& threshold < 1*/ ) 
		{
			//printf( "hi %lf %lf %.1lf\n", a.similarity, b.similarity, threshold ) ;
			if ( a.seqEnd >= seqs[ a.seqIdx ].consensusLen - gapAllow /*&& a.readEnd >= seqs[a.seqIdx].consensusLen */
					&& b.seqEnd >= seqs[ b.seqIdx ].consensusLen - gapAllow /*&& b.readEnd >= seqs[b.seqIdx].consensusLen*/ ) 
			{
				//matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
				if ( a.similarity - 0.1 > b.similarity && (a.matchCnt > b.matchCnt - 20))
				{
					// Check whether a's high similarity coming from the sequence it matched is a 
					bool directlyBetter = true ;
					if ( a.seqEnd - a.seqStart < b.seqEnd - b.seqStart )
					{
						int mismatchCnt = 0 ;
						int i, j ;
						for ( i = a.seqEnd, j = b.seqEnd ; i >= a.seqStart ; --i, --j )
							if ( seqs[ a.seqIdx ].consensus[i] != seqs[ b.seqIdx ].consensus[j] )
								++mismatchCnt ;
						if ( mismatchCnt <= 1 )
							directlyBetter = false ;
					}
					//printf( "%s %s %d\n", seqs[ a.seqIdx ].name, seqs[b.seqIdx].name, directlyBetter ) ;
					if ( directlyBetter )
						return true ;
				}
				else if ( a.similarity + 0.1 < b.similarity && a.matchCnt <= b.matchCnt - 20 )
					return false ;
			}
			else if (  a.seqEnd >= seqs[ a.seqIdx ].consensusLen - gapAllow && a.readEnd >= seqs[a.seqIdx].consensusLen 
				&& threshold < 1 ) 
			{	
				return true ;
			}
		}
		
		//if ( threshold == 1 )
		//	printf( "%s %d %lf %s %d %lf\n", seqs[ a.seqIdx ].name, a.matchCnt, a.similarity,
		//		 seqs[ b.seqIdx ].name, b.matchCnt, b.similarity ) ;

		if ( a.readStart == b.readStart && a.readEnd == b.readEnd )
		{
			if ( a.similarity > b.similarity )
				return true ;
			else if ( a.similarity < b.similarity )
				return false ;
			else 
			{
				// The assignment is the same, then we should replace the paralog.
				int i ;
				char *name = seqs[b.seqIdx].name ;
				for ( i = 0 ; name[i + 1] ; ++i )
				{
					if ((name[i + 1] == '-' || name[i + 1] == '*') && (name[i] < '0' || name[i] > '9'))
						return true ;
					if (name[i] == 'O' && name[i + 1] == 'R')
						return true ;
				}
			}
		}
		
		// Different alleles may have some effects on the boundary.
		// Then we pick the one with no less coverage portion and higher similiarty.
		if (threshold == 1.0 && IsSameGeneAllele(seqs[a.seqIdx].name, seqs[ b.seqIdx].name))
		{
			if ( ( a.seqEnd - a.seqStart + 1 ) / (double)(seqs[a.seqIdx].consensusLen) 
				>= (b.seqEnd - b.seqStart + 1) / (double)seqs[b.seqIdx].consensusLen
				&& a.similarity > b.similarity)
			{
				return true ;
			}
		}

		if ( matchCnt > b.matchCnt * threshold )
			return true ;
		else if ( threshold < 1.0 && 
			( a.matchCnt + 10 >= b.matchCnt || (
				a.similarity > b.similarity + 0.01 && a.matchCnt + 2 * kmerLength >= b.matchCnt ) ) )
			return true ;
		else
			return false ;
	}
	
	int GetContigIntervals( char *read, SimpleVector<struct _pair> &contigs )
	{
		int i, j ;
		//sprintf( buffer, "%d", len ) ;
		for ( i = 0 ; read[i] ; )
		{
			int NCnt = 0 ;
			for ( j = i + 1 ; read[j] ; ++j )
			{
				//printf( "%c %d %d %d\n", read[j], i, j, NCnt ) ;
				if ( j >= i + gapN && read[j - gapN] == 'N' )
					--NCnt ;
				if ( read[j] == 'N' ) 
					++NCnt ;
				if ( NCnt >= gapN )
					break ;

			}

			struct _pair nc ; 
			nc.a = i ;
			if ( read[j] )
				nc.b = j - gapN ; 
			else
				nc.b = j - 1 ;
			contigs.PushBack( nc ) ;

			if ( !read[j] )
				break ;
			i = j + 1 ;
		}
		return contigs.Size() ;
	}

	int GetContigIdx( int pos, SimpleVector<struct _pair> &contigs )
	{
		int size = contigs.Size() ;
		int i ;
		for ( i = 0 ; i < size ; ++i )
			if ( contigs[i].a <= pos && pos <= contigs[i].b )
				return i ;
		return 0 ;
	}
	
	void ShiftAnnotations( int at, int shift, int baseChange, struct _overlap geneOverlap[4], //struct _overlap cdr[3],
	                        std::vector<struct _overlap> *secondaryGeneOverlaps )
	{
		int i ;
		for ( i = 0 ; i < 4 ; ++i )	
		{
			if ( geneOverlap[i].seqIdx == -1 )
				continue ;

			if ( geneOverlap[i].readStart <= at && at <= geneOverlap[i].readEnd )
			{
				geneOverlap[i].matchCnt += 2 * baseChange ;
				geneOverlap[i].similarity = double( geneOverlap[i].matchCnt ) 
					/ ( geneOverlap[i].readEnd - geneOverlap[i].readStart + 1 + shift
						+ geneOverlap[i].seqEnd - geneOverlap[i].seqStart + 1 ) ;
			}

			if ( geneOverlap[i].readStart >= at )
				geneOverlap[i].readStart += shift ;
			if ( geneOverlap[i].readEnd >= at )
				geneOverlap[i].readEnd += shift ;
			

		}

		/*for ( i = 0 ; i < 3 ; ++i )
		{
			//cdr[2].readStart = cdr3NewStart ;
			//cdr[2].readEnd = cdr3NewEnd ;
			if ( cdr[i].readStart >= at )
				cdr[i].readStart += shift ;
			if ( cdr[i].readEnd >= at )
				cdr[i].readEnd += shift ;
		}
		cdr[2].similarity = 0.01 ;*/

		if ( secondaryGeneOverlaps != NULL )
		{
			std::vector<struct _overlap> &overlaps = *secondaryGeneOverlaps ;
			int size = overlaps.size() ;
			for ( i = 0 ; i < size ; ++i )
			{
				if ( overlaps[i].readStart <= at && at <= overlaps[i].readEnd )
				{
					overlaps[i].matchCnt += 2 * baseChange ;
					overlaps[i].similarity = double( overlaps[i].matchCnt ) 
						/ ( overlaps[i].readEnd - overlaps[i].readStart + 1 + shift +  
								+ overlaps[i].seqEnd - overlaps[i].seqStart + 1 ) ;
				}

				if ( overlaps[i].readStart >= at )
					overlaps[i].readStart += shift ;
				if ( overlaps[i].readEnd >= at )
					overlaps[i].readEnd += shift ;
			}
		}
	}

	
	// Impute the CDR3 by extension to the anchor case.
	int ImputeAnchorCDR3( char *read, char *nr, struct _overlap geneOverlap[4], struct _overlap cdr[3], 
			std::vector<struct _overlap> *secondaryGeneOverlaps )
	{
		int i, j, k ;
		// Locate the insert position
		int insertAt = -1 ; // Every position >= insertAt will be shifted
		int insertLen = -1 ;
		int seqIdx = -1 ;
		int seqStart = -1 ;
		int newStart = cdr[2].readStart ;
		int newEnd = cdr[2].readEnd ;
		
		// Exact one of the anchor should be in the annotation
		bool vInAnchor = false ;
		bool jInAnchor = false ;
		if (  seqs[ geneOverlap[0].seqIdx ].info[2].a >= geneOverlap[0].seqStart
			&& seqs[ geneOverlap[0].seqIdx ].info[2].a + 2 <= geneOverlap[0].seqEnd ) 
				vInAnchor = true ;
		//if ( cdr[2].readStart >= geneOverlap[0].readStart && cdr[2].readStart <= geneOverlap[0].readEnd )
		//	vInAnchor = true ;

		if (  seqs[ geneOverlap[2].seqIdx ].info[2].a >= geneOverlap[2].seqStart
			&& seqs[ geneOverlap[2].seqIdx ].info[2].a + 2 <= geneOverlap[2].seqEnd ) 
				jInAnchor = true ;
		//if ( cdr[2].readEnd >= geneOverlap[0].readStart && cdr[2].readEnd <= geneOverlap[2].readEnd )
		//	jInAnchor = true ;
		SimpleVector<struct _pair> contigs ;
		int contigCnt = GetContigIntervals( read, contigs ) ;
		// Change gap nucleotide from 'N' to 'M'
		for ( i = 0 ; i < contigCnt - 1 ; ++i )
		{
			for ( j = contigs[i].b + 1 ; j < contigs[i + 1].a ; ++j )
				read[j] = 'M' ;
		}
		// Double check whether anchor is in the gap 
		bool vAnchorInGap = false ;
		bool jAnchorInGap = false ;
		if ( vInAnchor ) 
		{
			int dest = geneOverlap[0].readEnd - ( geneOverlap[0].seqEnd - seqs[geneOverlap[0].seqIdx].info[2].a ) ;
			for ( i = geneOverlap[0].readEnd ; i >= dest ; --i )
			{
				if ( read[i] == 'M' )
				{
					vInAnchor = false ;
					vAnchorInGap = true ;
					break ;
				}
			}
		}

		if ( jInAnchor )
		{
			int dest = geneOverlap[2].readStart + ( seqs[geneOverlap[2].seqIdx].info[2].a + 2 - geneOverlap[2].seqStart) ;
			for ( i = geneOverlap[2].readStart ; i <= dest ; ++i )
			{
				if ( read[i] == 'M' )
				{
					jInAnchor = false ;
					jAnchorInGap = true ;
					break ;
				}
			}
		}

		if ( !vInAnchor ) 
		{
			seqIdx = geneOverlap[0].seqIdx ;
			int seqOffset = -1 ; // Locate the assembled CDR3 inside the V gene, where the imputation starts 
			int readOffset = -1 ;
			// seqOffset is the coordinate on the v gene that match with readOffset
			// Locate offset.
			if ( geneOverlap[0].seqEnd < seqs[ seqIdx ].info[2].a )
			{
				// V]...[CDR3]
				int matchLen ;
				int offset = AlignAlgo::LocatePartialSufPrefExactMatch( 
					seqs[ seqIdx ].consensus + seqs[ seqIdx ].info[2].a, seqs[ seqIdx ].consensusLen - seqs[ seqIdx ].info[2].a,
					read + cdr[2].readStart, cdr[2].readEnd - cdr[2].readStart + 1, 5, matchLen ) ;
				if ( offset != -1 )
				{
					seqOffset = offset + seqs[ seqIdx ].info[2].a ;
					readOffset = cdr[2].readStart ;
				}
				/*for ( k = seqs[ seqIdx ].info[2].a ; k + 9 < seqs[ seqIdx ].consensusLen ; ++k )
				{
					for ( i = k, j = cdr[2].readStart ; i < seqs[ seqIdx ].consensusLen, j <= cdr[2].readEnd ; ++i, ++j )
					{
						if ( seqs[ seqIdx ].consensus[i] != read[j] )
							break ;
					}
					if ( i - k >= 9 ) // We got a hit
					{
						seqOffset = k ;
						readOffset = cdr[2].readStart ;
						break ;
					}
				}*/

			}
			else
			{
				if ( vAnchorInGap )
				{
					//[V..NN[CDR3 NN..  V]     CDR3]
					int contigIdx = GetContigIdx( geneOverlap[0].readEnd, contigs ) ;
					readOffset = contigs[contigIdx].a ; 
					seqOffset = geneOverlap[0].seqEnd - (geneOverlap[0].readEnd - readOffset) ;
				}
				else
				{
					// [CDR3 [V..]  ]
					// This should happen in boundary case, otherwise the extension in annotateread should 
					// 	fix it.
					seqOffset = geneOverlap[0].seqStart ;
					readOffset = geneOverlap[0].readStart ;
				}
			}
			if ( seqOffset != -1 )
			{
				// There are could be some overhang nucleotide outside of the annotated CDR3 region.
				bool valid = true ;
				for ( i = seqOffset - 1, j = readOffset - 1 ; i >= seqs[ seqIdx ].info[2].a && j >= 0 ; --i, --j )
				{
					if ( read[j] == 'M' )
						break ;
					if ( seqs[ seqIdx ].consensus[i] != read[j] )
						valid = false ;
				}
				if ( valid )
				{
					// Valid for impute.
					insertAt = j + 1 ;
					insertLen = i - seqs[ seqIdx ].info[2].a + 1 ;
					seqStart = seqs[ seqIdx ].info[2].a ;
					newStart = insertAt ;
					newEnd += insertLen ;
				}
			}
		}
		else if ( !jInAnchor ) 
		{
			// For the j side.
			seqIdx = geneOverlap[2].seqIdx ;
			int seqOffset = -1 ; 
			int readOffset = -1 ;

			if ( geneOverlap[2].seqStart > seqs[ seqIdx ].info[2].a )
			{
				// [CDR3]...[J
				int matchLen ;
				int offset = AlignAlgo::LocatePartialSufSufExactMatch( seqs[ seqIdx ].consensus, seqs[ seqIdx ].info[2].a + 3,
					read + cdr[2].readStart, cdr[2].readEnd - cdr[2].readStart + 1, 5, matchLen ) ;
				if ( offset != -1 )
				{
					seqOffset = offset + matchLen - 1 ;
					readOffset = cdr[2].readEnd ;
				}
				/*for ( k = seqs[ seqIdx ].info[2].a + 2 ; k >= 8 ; --k )
				{
					for ( i = k, j = cdr[2].readEnd ; i >= 0 && j >= 0 ; --i, --j )
					{
						if ( seqs[ seqIdx ].consensus[i] != read[j] )
							break ;
					}
					if ( k - i >= 9 ) // We got a hit
					{
						seqOffset = k ;
						readOffset = cdr[2].readEnd ;
						break ;
					}
				}*/
			}
			else
			{
				if ( jAnchorInGap )
				{
					int contigIdx = GetContigIdx( geneOverlap[2].readStart, contigs ) ;
					readOffset = contigs[contigIdx].b ; 
					seqOffset = geneOverlap[2].seqStart + ( readOffset - geneOverlap[2].readStart ) ;
				}
				else
				{
					readOffset = geneOverlap[2].readEnd ;
					seqOffset = geneOverlap[2].seqEnd ; 
				}
			}

			if ( seqOffset != -1 )
			{
				// There are could be some overhang nucleotide outside of the annotated CDR3 region.
				bool valid = true ;
				for ( i = seqOffset + 1, j = readOffset + 1 ; i <= seqs[ seqIdx ].info[2].a + 2 && read[j] ; ++i, ++j )
				{
					if ( read[j] == 'M' )
						break ;
					if ( seqs[ seqIdx ].consensus[i] != read[j] )
						valid = false ;
				}
				if ( valid )
				{
					// Valid for impute.
					insertAt = j ;
					seqStart = i ;
					insertLen = seqs[ seqIdx ].info[2].a + 2 - seqStart + 1 ;
					newEnd = insertAt + insertLen - 1 ;
				}
			}
		}

		for ( i = 0 ; i < contigCnt - 1 ; ++i )
		{
			for ( j = contigs[i].b + 1 ; j < contigs[i + 1].a ; ++j )
				read[j] = 'N' ;
		}
		if ( insertLen > 0 )
		{
			// Shift the sequence.
			int len = strlen( read ) ; 
			strcpy( nr, read ) ;
			for ( i = len + insertLen ; i >= insertAt + insertLen ; --i )
				nr[i] = nr[i - insertLen] ;
			
			// Put in the imputed seqauence
			for ( i = seqStart, j = insertAt ; j < insertAt + insertLen ; ++i, ++j )
				nr[j] = seqs[ seqIdx ].consensus[i] ;

			cdr[2].readStart = newStart ;
			cdr[2].readEnd = newEnd ;
			cdr[2].similarity = 0.01 ;
			// Shift what we found.
			ShiftAnnotations( insertAt, insertLen, 0, geneOverlap, secondaryGeneOverlaps ) ;
		}
		else if ( insertLen == 0 )
		{
			cdr[2].readStart = newStart ;
			cdr[2].readEnd = newEnd ;
			cdr[2].similarity = 0.5 ;
			return -1 ;	
		}
		if ( insertLen < -1 )
			insertAt = -1 ;
		return insertAt ;

	}

	int ImputeInternalCDR3( char *read, char *nr, struct _overlap geneOverlap[4], struct _overlap cdr[3], 
			std::vector<struct _overlap> *secondaryGeneOverlaps )
	{
		int i, j, k ;
		if ( geneOverlap[0].seqIdx == -1 || geneOverlap[2].seqIdx == -1 )
			return -1 ;
		int vSeqIdx = geneOverlap[0].seqIdx ;
		int jSeqIdx = geneOverlap[2].seqIdx ;
		if ( seqs[ vSeqIdx ].info[2].a == -1 || seqs[ jSeqIdx ].info[2].a == -1 )
			return -1 ;

		SimpleVector<struct _pair> contigs ;
		int contigCnt = GetContigIntervals( read, contigs ) ;
		
		// Make sure there is exactly one gap.
		int gapCnt = 0 ;
		int gapStart, gapEnd ;
		for ( i = 0 ; i < contigCnt - 1 ; ++i )
		{
			if ( contigs[i].b >= cdr[2].readStart && contigs[i].b <= cdr[2].readEnd 
				&& contigs[i + 1].a >= cdr[2].readStart && contigs[i + 1].a <= cdr[2].readEnd )
			{
				gapStart = contigs[i].b + 1 ;
				gapEnd = contigs[i + 1].a - 1 ;
				++gapCnt ;
			}
		}
		if ( gapCnt != 1 )
			return -1 ;
		
		// Decide which gene this gap interrupts
		int vOffset = -1, jOffset = -1 ;
		int vMatchLen, jMatchLen ;
		vOffset = AlignAlgo::LocatePartialSufPrefExactMatch( 
			seqs[ vSeqIdx ].consensus + seqs[ vSeqIdx ].info[2].a, seqs[ vSeqIdx ].consensusLen - seqs[ vSeqIdx ].info[2].a,
			read + gapEnd + 1, cdr[2].readEnd - gapEnd, 5, vMatchLen ) ;
		jOffset = AlignAlgo::LocatePartialSufSufExactMatch( seqs[ jSeqIdx ].consensus, seqs[ jSeqIdx ].info[2].a + 3, 
			read + cdr[2].readStart, gapStart - cdr[2].readStart, 5, jMatchLen ) ;
		
		if ( ( vOffset != -1 && jOffset != -1 ) 
			|| ( vOffset == -1 && jOffset == -1 ) )
			return -1 ;

		struct _pair anchor[2] ; // the coordinate info on the anchor point of the gap :a: seq, b:read 
		int seqIdx = -1 ;
		if ( vOffset != -1 ) // Gap interrupt the V gene
		{
			bool valid = true ;
			struct _seqWrapper &seq = seqs[ vSeqIdx ] ;
			for ( i = seq.info[2].a, j = cdr[2].readStart ; i < seq.consensusLen, j < gapStart ; ++i, ++j )
				if ( seq.consensus[i] != read[j] )
					valid = false ;
			if ( valid == false )
				return -1 ;
			anchor[0].a = i - 1 ;
			anchor[0].b = j - 1 ;
			anchor[1].a = vOffset + seqs[ vSeqIdx ].info[2].a ;
			anchor[1].b = gapEnd + 1 ;
			seqIdx = vSeqIdx ;
		}
		else // jOffset != -1 
		{
			bool valid = true ;
			struct _seqWrapper &seq = seqs[ jSeqIdx ] ;
			for ( i = seq.info[2].a + 2, j = cdr[2].readEnd ; i >= 0 && j > gapEnd ; --i, --j )
				if ( seq.consensus[i] != read[j] )
					valid = false ;
			if ( valid == false )
				return -1 ;
			anchor[0].a = jOffset + jMatchLen - 1 ;
			anchor[0].b = gapStart - 1 ;
			anchor[1].a = i + 1 ;
			anchor[1].b = j + 1 ;
			seqIdx = jSeqIdx ;
		}
		
		// Rearrange the sequence
		// Record where anchor[1] got shifted to
		int shiftAt ;
		int shift ;
		int baseChange = 0 ; // The change t
		for ( j = 0 ; j <= anchor[0].b ; ++j )
			nr[j] = read[j] ;
		if ( anchor[1].a > anchor[0].a ) 
		{
			// Put in the imputed sequence
			for ( i = anchor[0].a + 1 ; i < anchor[1].a ; ++i, ++j )
				nr[j] = seqs[ seqIdx ].consensus[i] ;
			shiftAt = anchor[1].b ;
			shift = j - anchor[1].b ; // shift could be negative is the imputed portion is shorter than gap.
			for ( j = anchor[1].b ; read[j] ; ++j )
				nr[j + shift] = read[j] ;
			baseChange = anchor[1].a - anchor[0].a - 1 ;
		}
		else
		{
			// The two anchor portion overlaps.
			shiftAt = anchor[1].b ;
			int overlap = anchor[0].a - anchor[1].a + 1 ;
			shift = ( anchor[0].b - overlap + 1 ) - anchor[1].b ;
			for ( j = anchor[1].b + overlap ; read[j] ; ++j ) // anchor[1].b+overlap should not excess the read len.
				nr[j + shift] = read[j] ;
			baseChange = -overlap ;
		}
		nr[j + shift] = '\0' ; 
		cdr[2].readEnd += shift ;
		cdr[2].similarity = 0.01 ;
		ShiftAnnotations( shiftAt, shift, baseChange, geneOverlap, secondaryGeneOverlaps ) ;

		return shiftAt ;
	}

	int ImputeCDR3( char *read, char *nr, struct _overlap geneOverlap[4], struct _overlap cdr[3], 
			std::vector<struct _overlap> *secondaryGeneOverlaps )
	{
		if ( cdr[2].seqIdx == -1 || cdr[2].similarity != 0  
				|| geneOverlap[0].seqIdx == -1 || geneOverlap[2].seqIdx == -1 
				|| seqs[ geneOverlap[0].seqIdx ].info[2].a == -1 || seqs[ geneOverlap[2].seqIdx ].info[2].a == -1 
				|| geneOverlap[0].readEnd >= geneOverlap[2].readStart ) 
			return -1 ;
		if ( seqs[ geneOverlap[0].seqIdx ].name[0] != 'T' )
			return -1 ;

		// Exact one of the anchor should be in the annotation
		bool vInAnchor = false ;
		bool jInAnchor = false ;
		if (  seqs[ geneOverlap[0].seqIdx ].info[2].a >= geneOverlap[0].seqStart
				&& seqs[ geneOverlap[0].seqIdx ].info[2].a + 2 <= geneOverlap[0].seqEnd ) 
			vInAnchor = true ;
		//if ( cdr[2].readStart >= geneOverlap[0].readStart && cdr[2].readStart <= geneOverlap[0].readEnd )
		//	vInAnchor = true ;

		if (  seqs[ geneOverlap[2].seqIdx ].info[2].a >= geneOverlap[2].seqStart
				&& seqs[ geneOverlap[2].seqIdx ].info[2].a + 2 <= geneOverlap[2].seqEnd ) 
			jInAnchor = true ;
		//if ( cdr[2].readEnd >= geneOverlap[0].readStart && cdr[2].readEnd <= geneOverlap[2].readEnd )
		//	jInAnchor = true ;
		int ret = -1 ;
		if ( vInAnchor && jInAnchor )	
		{
			int j ;
			for ( j = cdr[2].readStart ; j <= cdr[2].readEnd ; ++j )
				if ( read[j] == 'N' && read[j + 1] == 'N' )
					break ;
			if ( j <= cdr[2].readEnd )
				// There is a gap in side the CDR3 region.
				ret = ImputeInternalCDR3( read, nr, geneOverlap, cdr, secondaryGeneOverlaps ) ;
			else
			{
				ret = ImputeAnchorCDR3( read, nr, geneOverlap, cdr, secondaryGeneOverlaps ) ;
			}
		}
		else if ( vInAnchor || jInAnchor )
		{
			int j ;
			for ( j = cdr[2].readStart ; j <= cdr[2].readEnd ; ++j )
				if ( read[j] == 'N' )
					return -1 ;
			ret = ImputeAnchorCDR3( read, nr, geneOverlap, cdr, secondaryGeneOverlaps ) ;
		}
		
		if ( ret != -1 )
			AnnotateReadDGene( nr, geneOverlap, cdr, secondaryGeneOverlaps ) ;
		return ret ;
	}
	
	// Annotate the D gene
	int AnnotateReadDGene( char *read, struct _overlap geneOverlap[4], struct _overlap cdr[3], 
		std::vector<struct _overlap> *secondaryGeneOverlaps )
	{
		int i, j, k ;
		std::vector< std::vector<struct _overlap> > overlaps ;
		SimpleVector<int> dGeneSeqIdx ;

		int seqCnt = seqs.size() ;
		int anchorSeqIdx = -1 ;
		
		if ( cdr[2].seqIdx == -1 || cdr[2].similarity == 0 )
			return -1 ;

		if ( geneOverlap[0].seqIdx != -1 )
			anchorSeqIdx = geneOverlap[0].seqIdx ;
		else if ( geneOverlap[2].seqIdx != -1 )
			anchorSeqIdx = geneOverlap[2].seqIdx ;
		else
			return -1 ;

		// Only consider IGH, TRB, TRD
		if ( seqs[ anchorSeqIdx ].name[2] != 'H' && seqs[ anchorSeqIdx ].name[2] != 'B' 
			&& seqs[ anchorSeqIdx ].name[2] != 'D' )
			return -1 ;
		
		int maxDlen = 0 ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].isRef && GetGeneType( seqs[i].name ) == 1 
				&& ( seqs[i].name[0] == seqs[anchorSeqIdx].name[0] )
				&& ( seqs[i].name[2] == seqs[anchorSeqIdx].name[2] ) 
				&& ( seqs[i].name[1] == seqs[anchorSeqIdx].name[1] ) )
			{
				if ( seqs[i].consensusLen > maxDlen )
					maxDlen = seqs[i].consensusLen ;
				dGeneSeqIdx.PushBack( i ) ;
			}
		}
		
		seqCnt = dGeneSeqIdx.Size() ;
		std::vector<struct _overlap> dOverlaps ;
		int cdr3Len = cdr[2].readEnd - cdr[2].readStart + 1 ;
		char *align = new char[ cdr3Len ] ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int seqStart ;
			int readStart ;
			int seqIdx = dGeneSeqIdx[i] ;
			int alignScore = AlignAlgo::LocalAlignment( seqs[ seqIdx ].consensus, seqs[ seqIdx ].consensusLen,
				read + cdr[2].readStart, cdr3Len, seqStart, readStart, align ) ;

			if ( alignScore >= 5 * SCORE_MATCH_LOCAL )
			{
				readStart += cdr[2].readStart ;
				int readEnd = readStart - 1 ;
				int seqEnd = seqStart - 1 ;
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( align[j] != EDIT_INSERT )
						++seqEnd ;
					if ( align[j] != EDIT_DELETE )
						++readEnd ;
				}

				//printf( "%s\n", seqs[seqIdx].name ) ;
				//AlignAlgo::VisualizeAlignment( seqs[ seqIdx ].consensus + seqStart, seqEnd - seqStart + 1,
				//	read + readStart, readEnd - readStart + 1, align ) ;
				
				int matchCnt, mismatchCnt, indelCnt ;
				GetAlignStats( align, false, matchCnt, mismatchCnt, indelCnt ) ;

				struct _overlap no ;
				no.seqIdx = seqIdx ;
				no.seqStart = seqStart ; no.seqEnd = seqEnd ;
				no.readStart = readStart ; no.readEnd = readEnd;
				no.matchCnt = 2 * matchCnt ;
				no.similarity = (double)no.matchCnt / ( seqEnd - seqStart + 1 + readEnd - readStart + 1 ) ;
								
				dOverlaps.push_back( no ) ;
			}
		}
		delete[] align ;

		if ( dOverlaps.size() == 0 )
			return -1 ;
		
		int bestTag = 0 ;
		int overlapCnt = dOverlaps.size() ;
		//std::sort( dOverlaps.begin(), dOverlaps.end() ) ;
		for ( i = 1 ; i < overlapCnt ; ++i )
		{
			if ( IsBetterGeneMatch( dOverlaps[i], dOverlaps[bestTag], 1.0) )
			{
				bestTag = i ;
			}
		}
		geneOverlap[1] = dOverlaps[bestTag] ;
		return dOverlaps[bestTag].seqIdx ;
	}

	// Figure out the gene composition for the read. 
	// Return successful or not.
	int AnnotateRead( char *read, int detailLevel, struct _overlap geneOverlap[4], struct _overlap cdr[3], 
		std::vector<struct _overlap> *secondaryGeneOverlaps )
	{
		int i, j, k ;
		
		std::vector< std::vector< struct _overlap > > contigOverlaps ;
		SimpleVector<struct _pair> contigs ; // the range of contigs
		int contigCnt ;

		std::vector<struct _overlap> overlaps ;
		std::vector<struct _overlap> allOverlaps ;
		int overlapCnt ;
		
		char BT = '\0' ; // Bcell, Tcell
		char chain = '\0' ;
		int len = strlen( read ) ;
		//buffer[0] = '\0' ;

		geneOverlap[0].seqIdx = geneOverlap[1].seqIdx = geneOverlap[2].seqIdx = geneOverlap[3].seqIdx = -1 ;
		if ( detailLevel >= 2 )
			cdr[0].seqIdx = cdr[1].seqIdx = cdr[2].seqIdx = -1 ;

		contigCnt = GetContigIntervals( read, contigs ) ;
		char *contigBuffer = new char[len + 1] ;

		int locatePartialMatchMinLen = 8 ; // The minimum length to detect those small hits

		//if ( detailLevel > 0 )
		// Obtain the overlaps for each contig
		contigOverlaps.resize( contigCnt ) ;
		for ( k = 0 ; k < contigCnt ; ++k )
		{
			int contigLen = contigs[k].b - contigs[k].a + 1 ;
			memcpy( contigBuffer, read + contigs[k].a, contigLen ) ;
			contigBuffer[ contigLen ] = '\0' ;
			contigOverlaps[k].clear() ;

			int contigOverlapCnt = GetOverlapsFromRead( contigBuffer, 0, -1, detailLevel == 0 ? 0 : 1, false, contigOverlaps[k] ) ;		
			for ( i = 0 ; i < contigOverlapCnt ; ++i )
			{
				contigOverlaps[k][i].readStart += contigs[k].a ;
				contigOverlaps[k][i].readEnd += contigs[k].a ;
			}
			std::sort( contigOverlaps[k].begin(), contigOverlaps[k].end() ) ;
			
			//std::vector<struct _overlap> &overlaps = contigOverlaps[k] ;
			//for ( i = 0 ; i < contigOverlapCnt ; ++i )
			//	printf( "%d: %d %s %lf %d: %d %d. %d %d\n", k, i, seqs[ overlaps[i].seqIdx].name, overlaps[i].similarity, overlaps[i].matchCnt, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ;
		}
		delete[] contigBuffer ;


		// Extend overlaps that can cross contigs
		int *seqUsed = new int[ seqs.size() ] ;
		if ( detailLevel >= 1 )
		{
			std::vector< std::vector<struct _overlap> > extendedOverlaps = contigOverlaps ; 
			for ( k = 0 ; k < contigCnt ; ++k )	
			{
				memset( seqUsed, -1, sizeof( int ) * seqs.size() ) ;
				int cnt = contigOverlaps[k].size() ;
				std::vector<struct _overlap> &ovs = extendedOverlaps[k] ;
				//for ( i = 0 ; i < cnt ; ++i )
				//	printf( "%d: %d %s %lf %d: %d %d. %d %d\n", k, i, seqs[ ovs[i].seqIdx].name, ovs[i].similarity, ovs[i].matchCnt, ovs[i].readStart, ovs[i].readEnd, ovs[i].seqStart, ovs[i].seqEnd ) ;

				for ( i = 0 ; i < cnt ; ++i )
				{
					if ( seqUsed[ ovs[i].seqIdx ] != -1 )
						continue ;
					
					int effectiveLen = ovs[i].readEnd - ovs[i].readStart + 1 
						+ ovs[i].seqEnd - ovs[i].seqStart +1 ;
					
					for ( j = k - 1 ; j >= 0 ; --j )
					{
						bool extended = false  ;
						int cnt2 = contigOverlaps[j].size() ;
						int l ;
						for ( l = 0 ; l < cnt2 ; ++l )
						{
							struct _overlap &o = contigOverlaps[j][l] ;
							if ( o.seqIdx == ovs[i].seqIdx )
							{
								if ( o.seqEnd < ovs[i].seqStart + 31 
									&& o.readStart <= contigs[j + 1].a + 10 
									&& ovs[i].readEnd <= contigs[j].b - 10 ) 
								{
									ovs[i].readStart = o.readStart ;
									ovs[i].seqStart = o.seqStart ;
									ovs[i].matchCnt += o.matchCnt ;
									effectiveLen += o.readEnd - o.readStart + 1 
										+ o.seqEnd - o.seqStart + 1 ;
									extended = true ;
								}
								//printf( "<=+ %s: %d %d. %d %d\n", seqs[ ovs[i].seqIdx ].name, ovs[i].matchCnt, effectiveLen, o.matchCnt, o.readEnd - o.readStart + 1 + o.seqEnd - o.seqStart + 1 ) ;
								break ;
							}
						}
						
						if ( !extended )
							break ;
					}

					for ( j = k + 1 ; j < contigCnt ; ++j )
					{
						bool extended = false ;
						int cnt2 = contigOverlaps[j].size() ;
						int l ;
						for ( l = 0 ; l < cnt2 ; ++l )
						{
							struct _overlap &o = contigOverlaps[j][l] ;
							if ( o.seqIdx == ovs[i].seqIdx )
							{
								if ( o.seqStart > ovs[i].seqEnd - 31 
									&& o.readEnd <= contigs[j - 1].b - 10 
									&& ovs[i].readStart <= contigs[j].a + 10 ) 
								{
									ovs[i].readEnd = o.readEnd ;
									ovs[i].seqEnd = o.seqEnd ;
									ovs[i].matchCnt += o.matchCnt ;
									effectiveLen += o.readEnd - o.readStart + 1 
										+ o.seqEnd - o.seqStart + 1 ;
									extended = true ;
								}
								//printf( "=>+ %s: %d %d. %d %d\n", seqs[ ovs[i].seqIdx ].name, ovs[i].matchCnt, effectiveLen, o.matchCnt, o.readEnd - o.readStart + 1 + o.seqEnd - o.seqStart + 1 ) ;
								break ;
							}
						}

						if ( !extended )
							break ;
					}
					ovs[i].similarity = (double)ovs[i].matchCnt / effectiveLen ;
				//	printf( "%s: %lf\n", seqs[ ovs[i].seqIdx ].name, ovs[i].similarity ) ;
					seqUsed[ ovs[i].seqIdx ] = i ;
				}
			}

			contigOverlaps = extendedOverlaps ;
		}
		
		// Put the extended version of overlaps into overall overlaps
		for ( k = 0 ; k < contigCnt ; ++k )
		{
			overlaps.insert( overlaps.end(), contigOverlaps[k].begin(), contigOverlaps[k].end() ) ;
		}
		std::sort( overlaps.begin(), overlaps.end() ) ;
		k = 0 ;
		memset( seqUsed, -1, sizeof( int ) * seqs.size() ) ;	
		overlapCnt = overlaps.size() ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			// Remove the hits on the D gene
			if ( GetGeneType( seqs[ overlaps[i].seqIdx ].name ) == 1 )
				continue ;

			// Remove those overlaps that was secondary as well.
			if ( seqUsed[ overlaps[i].seqIdx ] == -1 && overlaps[i].similarity >= 0.8 )
			{
				seqUsed[ overlaps[i].seqIdx ] = k ; // Store the index where this is copied to
				overlaps[k] = overlaps[i] ;
				++k ;
			}
			else if ( seqUsed[overlaps[i].seqIdx] != -1 && GetGeneType( seqs[overlaps[i].seqIdx].name) == 2)
			{
				struct _overlap &baseline = overlaps[ seqUsed[overlaps[i].seqIdx] ] ;
				if ( overlaps[i].matchCnt == baseline.matchCnt && overlaps[i].similarity == baseline.similarity )
				{
					// For equally good hit, we need to break the tie.
					// For simplicity, we only do this for J gene, since it is short.
					for ( j = 0 ; j < k ; ++j )
					{
						if (GetGeneType( seqs[overlaps[j].seqIdx].name ) == 3)
							break ;
					}
					if ( j < k )
					{
						if (overlaps[i].readEnd <= overlaps[j].readStart + 3 )
						{
							if ( baseline.readEnd > overlaps[j].readStart + 3
								|| (ABS(overlaps[i].readEnd - overlaps[j].readStart) <
									ABS(baseline.readEnd - overlaps[j].readStart)))
								baseline = overlaps[i] ;
						}
					}
				}
			}
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		delete[] seqUsed ; 
		
		if ( overlapCnt == 0 )
		{
			//if ( detailLevel >= 2 )
			//	sprintf( buffer + strlen( buffer ), " * * * CDR1(0-0):0.00=null CDR2(0-0):0.00=null CDR3(0-0):0.00=null" ) ;
			return 0 ;
		}
		// Get the coverage of the genes.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *name = seqs[ overlaps[i].seqIdx ].name ;
			if ( BT && name[0] != BT )
				continue ;
			BT = name[0] ;
			if ( chain && !( name[2] == chain 
				|| ( name[2] == 'D' && chain == 'A' ) 
				|| ( name[2] == 'A' && chain == 'D' ) ) )
				continue ;
			chain = name[2] ;

			int geneType = GetGeneType( name ) ;
			//printf( "%d %s %d %lf %d: %d %d. %d %d\n", i, seqs[ overlaps[i].seqIdx].name, overlaps[i].strand, overlaps[i].similarity, overlaps[i].matchCnt, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ;
			/*if ( overlaps[i].similarity < 0.8 )//&& geneType == 0 
				//&& ( overlaps[i].readEnd - overlaps[i].readStart + 1 <= 40 ||
				//	overlaps[i].seqEnd - overlaps[i].seqStart + 1 <= 40 ) ) 
			{
				continue ;
			}*/

			if ( geneType >= 0 && geneOverlap[ geneType ].seqIdx == -1 )
				//IsBetterGeneMatch(overlaps[i], geneOverlap[geneType], 1.0 ) )
				geneOverlap[ geneType ] = overlaps[i] ;

			
			//printf( "%s %d\n", seqs[ overlaps[i].seqIdx ].name, overlaps[i].matchCnt ) ;
			if ( geneType >= 0 && ( IsBetterGeneMatch( overlaps[i], geneOverlap[geneType], 0.95 ) || 
				( geneOverlap[ geneType ].seqIdx != -1 && overlaps[i].similarity - 0.1 > geneOverlap[ geneType ].similarity ) ) 
			   )
			{
				allOverlaps.push_back( overlaps[i] ) ;
			}
			else if ( geneType >= 0 && geneOverlap[ geneType ].seqIdx != -1 &&
				( overlaps[i].readEnd < geneOverlap[geneType].readStart || 
					overlaps[i].readStart > geneOverlap[ geneType ].readEnd )
				&& IsBetterGeneMatch( overlaps[i], geneOverlap[geneType], 0.9 ) )
			{
				// The gene on a different region of the read. Might happen due to false alignment.
				allOverlaps.push_back( overlaps[i] ) ;
			}
		}
		
		// Check whether the match to the constant gene is random.
		if ( geneOverlap[3].seqIdx != -1 && geneOverlap[3].readEnd - geneOverlap[3].readStart + 1 <= len / 2 
			&& geneOverlap[3].readEnd - geneOverlap[3].readStart + 1 <= 50 )
		{
			for ( i = 0 ; i < 3 ; ++i )
			{
				if ( geneOverlap[i].seqIdx >= 0 && 
					( geneOverlap[i].readEnd - 17 > geneOverlap[3].readStart 
						|| geneOverlap[3].readEnd < geneOverlap[i].readEnd ) 
					&& geneOverlap[3].seqStart >= 100 )
				{
					// Filter out all the overlaps from C gene
					geneOverlap[3].seqIdx = -1 ;
					break ;
				}
			}

			if ( i < 3 && detailLevel >= 1 )
			{
				int size = allOverlaps.size() ;
				for ( i = 0, k = 0 ; i < size; ++i )
				{
					int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
					if ( geneType != 3 )
					{
						allOverlaps[k] = allOverlaps[i] ;
						++k ;
					}
				}
				allOverlaps.resize( k ) ;
			}
		}
		// Check inconsistent gene hits.
		if ( detailLevel >= 1 )
		{
			for ( i = 0 ; i < 4 ; ++i )
			{
				if ( i == 1 || geneOverlap[i].seqIdx == -1 )
					continue ;
				for ( j = 0 ; j < 4 ; ++j )
				{
					if ( j == 1 || i == j || geneOverlap[j].seqIdx == -1 )
						continue ;
					if ( ( j < i && geneOverlap[j].readEnd > geneOverlap[i].readEnd ) 
						|| ( j > i && geneOverlap[i].readEnd > geneOverlap[j].readEnd ) )
					{
						int removeType = i ;
						if ( geneOverlap[j].similarity < geneOverlap[i].similarity )
							removeType = j ;
						
						geneOverlap[ removeType ].seqIdx = -1 ;
						int size = allOverlaps.size() ;
						for ( i = 0, k = 0 ; i < size; ++i )
						{
							int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
							if ( geneType != removeType )
							{
								allOverlaps[k] = allOverlaps[i] ;
								++k ;
							}
						}
						allOverlaps.resize( k ) ;

						break ;
					}
				}
			}

			// Remove other secondary hits with inconsistent coordinate.
			int size = allOverlaps.size() ;
			for ( i = 0, k = 0 ; i < size ; ++i )
			{
				int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
				if ( allOverlaps[i].readEnd <= geneOverlap[ geneType ].readStart 
					|| allOverlaps[i].readStart >= geneOverlap[ geneType].readEnd )
					continue ;
				allOverlaps[k] = allOverlaps[i] ;
				++k ;
			}
			allOverlaps.resize( k ) ;
		}
		
		// Extend overlap
		if ( detailLevel >= 2 )
		{
			// Change gap nucleotide from 'N' to 'M'
			for ( i = 0 ; i < contigCnt - 1 ; ++i )
			{
				for ( j = contigs[i].b + 1 ; j < contigs[i + 1].a ; ++j )
					read[j] = 'M' ;
			}
		}

		if ( detailLevel >= 1 )
		{
			char *align = new char[ 2 * len + 2 ] ;
			char *rvr = new char[len + 1] ;
			int size = allOverlaps.size() ;
			for ( i = 0 ; i < size ; ++i )
			{
				// Extend right.
				int seqIdx = allOverlaps[i].seqIdx ;				
				AlignAlgo::GlobalAlignment_OneEnd( seqs[ seqIdx ].consensus + allOverlaps[i].seqEnd + 1, 
					seqs[ seqIdx ].consensusLen - allOverlaps[i].seqEnd - 1, 
					read + allOverlaps[i].readEnd + 1, len - allOverlaps[i].readEnd - 1, 0, align ) ;
				//printf( "%s\n", seqs[ seqIdx ].name ) ;
				//AlignAlgo::VisualizeAlignment( seqs[ seqIdx ].consensus + allOverlaps[i].seqEnd + 1, 
				//	seqs[ seqIdx ].consensusLen - allOverlaps[i].seqEnd - 1, 
				//	read + allOverlaps[i].readEnd + 1, len - allOverlaps[i].readEnd - 1, align ) ;
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( read[ allOverlaps[i].readEnd + 1 ] == 'N' && read[ allOverlaps[i].readEnd + 2] == 'N' )
						break ;
					if ( align[j] == EDIT_MATCH || align[j] == EDIT_MISMATCH )
					{
						++allOverlaps[i].readEnd ;
						++allOverlaps[i].seqEnd ;
						
						if ( align[j] == EDIT_MATCH )
							allOverlaps[i].matchCnt += 2 ;
					}
					else if ( radius > 0 )
					{
						if ( align[j] == EDIT_INSERT )
							++allOverlaps[i].readEnd ;
						else if ( align[j] == EDIT_DELETE )
							++allOverlaps[i].seqEnd ;
					}
					else
						break ;
				}

				// Extend left.
				char *rvs = new char[seqs[ seqIdx ].consensusLen ] ;
				Reverse( rvr, read, allOverlaps[i].readStart ) ;
				Reverse( rvs, seqs[seqIdx].consensus, allOverlaps[i].seqStart ) ;
				//rvr[geneOverlap[i].readStart] = '\0' ;
				//rvs[geneOverlap[i].seqStart] = '\0' ;
					
				AlignAlgo::GlobalAlignment_OneEnd( rvs, allOverlaps[i].seqStart, rvr, allOverlaps[i].readStart, 0, align ) ;
				//AlignAlgo::VisualizeAlignment( rvs, geneOverlap[i].readStart, rvr, rvr[geneOverlap[i].readStart], align ) ;
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( allOverlaps[i].readStart > 1 && read[ allOverlaps[i].readStart - 1 ] == 'N' 
						&& read[ allOverlaps[i].readStart - 2 ] == 'N' )
					{
						break ;
					}
					if ( align[j] == EDIT_MATCH || align[j] == EDIT_MISMATCH )
					{
						--allOverlaps[i].readStart ;
						--allOverlaps[i].seqStart ;
					
						if ( align[j] == EDIT_MATCH )
							allOverlaps[i].matchCnt += 2 ;
					}
					else if ( radius > 0 )
					{
						if ( align[j] == EDIT_INSERT )
							--allOverlaps[i].readStart ;
						else if ( align[j] == EDIT_DELETE )
							--allOverlaps[i].seqStart ;
					}
					else
						break ;
				}
				delete[] rvs ;
				
				allOverlaps[i].similarity = (double)( allOverlaps[i].matchCnt ) / 
						( allOverlaps[i].seqEnd - allOverlaps[i].seqStart + 1 + 
							allOverlaps[i].readEnd - allOverlaps[i].readStart + 1 ) ;
			}

			if ( detailLevel >= 2 )
			{
				// See whether we can improve gene alignment information across contigs anchoring CDR3.
				// There might be short overhang portion near contig break point.
				for ( i = 0 ; i < size ; ++i )
				{
					int seqIdx = allOverlaps[i].seqIdx ;		
					
					int geneType = GetGeneType( seqs[seqIdx].name ) ;
					if ( geneType == 0 && read[ allOverlaps[i].readEnd + 1 ] == 'M' )
					{
						int contigIdx = GetContigIdx( allOverlaps[i].readEnd, contigs ) + 1 ;
						int matchLen = 0 ;
						int geneOffset = AlignAlgo::LocatePartialSufPrefExactMatch(
								seqs[seqIdx].consensus + allOverlaps[i].seqEnd + 1, 
								seqs[seqIdx].consensusLen - allOverlaps[i].seqEnd - 1, 
								read + contigs[contigIdx].a,
								contigs[contigIdx].b - contigs[contigIdx].a + 1,
								locatePartialMatchMinLen, matchLen ) ;
						if ( geneOffset != -1 )
						{
							int tmp = allOverlaps[i].seqEnd - allOverlaps[i].seqStart + 1 +  
								allOverlaps[i].readEnd - allOverlaps[i].readStart + 1 ;
							allOverlaps[i].readEnd = contigs[contigIdx].a + matchLen - 1 ;
							allOverlaps[i].seqEnd = allOverlaps[i].seqEnd + 1 + geneOffset + matchLen - 1 ;
							allOverlaps[i].matchCnt += 2 * matchLen ;
							allOverlaps[i].similarity = (double)(allOverlaps[i].matchCnt ) / (tmp + 2 * matchLen) ;
						}
					}
					else if ( geneType == 2 && allOverlaps[i].readStart > 0
							&& read[allOverlaps[i].readStart - 1] == 'M' )
					{
						int contigIdx = GetContigIdx( allOverlaps[i].readStart, contigs ) - 1 ;
						int matchLen = 0 ;
						
						int geneOffset = AlignAlgo::LocatePartialSufSufExactMatch( 
								seqs[seqIdx].consensus, allOverlaps[i].seqStart,
								read + contigs[contigIdx].a, 
								contigs[contigIdx].b - contigs[contigIdx].a + 1,
								locatePartialMatchMinLen, matchLen ) ;
						if ( geneOffset != -1 )
						{
							int tmp = allOverlaps[i].seqEnd - allOverlaps[i].seqStart + 1 +  
								allOverlaps[i].readEnd - allOverlaps[i].readStart + 1 ;
							allOverlaps[i].readStart = contigs[contigIdx].b - matchLen + 1 ;
							allOverlaps[i].seqStart = geneOffset ;
							allOverlaps[i].matchCnt += 2 * matchLen ;
							allOverlaps[i].similarity = (double)(allOverlaps[i].matchCnt ) / (tmp + 2 * matchLen) ;
						}
		
					}
				}
			}

			std::sort( allOverlaps.begin(), allOverlaps.end() ) ;
			
			for ( i = 0 ; i < 4 ; ++i )
			{
				geneOverlap[i].seqIdx = -1 ;
				geneOverlap[i].matchCnt = -1 ;
			}
			
			for ( i = 0 ; i < size ; ++i )
			{
				//printf( "%d %s %lf %d: %d %d\n", i, seqs[ allOverlaps[i].seqIdx].name, allOverlaps[i].similarity, allOverlaps[i].matchCnt, allOverlaps[i].readStart, allOverlaps[i].readEnd ) ;
				int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
				if ( IsBetterGeneMatch( allOverlaps[i], geneOverlap[ geneType ], 1.0 ) )
				{
					//if ( geneOverlap[ geneType].seqIdx != -1 )
					//	printf( "%s %s\n", seqs[ geneOverlap[ geneType ].seqIdx ].name, seqs[ allOverlaps[i].seqIdx ].name ) ;
					//else
					//	printf( "-1 %s\n", seqs[ allOverlaps[i].seqIdx ].name ) ;
						
					geneOverlap[ geneType ] = allOverlaps[i] ;
				}
			}

			// Rescue constant gene if the anchor is too short
			/*if ( geneOverlap[2].seqIdx != -1 && geneOverlap[3].seqIdx == -1 && 
				geneOverlap[2].readEnd + hitLenRequired - 1 >= len )
			{
				int seqCnt = seqs.size() ;
				int rstart = geneOverlap[2].readEnd + 6 ;
				int plen = len - rstart ;
				if ( plen >= 9 )
				{
					char *tBuffer = new char[len - geneOverlap[2].readEnd + 6] ;
					for ( j = 0 ; j < seqCnt ; ++j )
					{
						int geneType = GetGeneType( seqs[j].name ) ;
						if ( geneType != 3 )
							continue ;
						memcpy( tBuffer, seqs[j].consensus, len - geneOverlap[2].readEnd + 2 ) ;
						tBuffer[ len - geneOverlap[2].readEnd + 2 ] = '\0' ;
						char *p = strstr( tBuffer, read + rstart ) ;
						if ( p == NULL )
							continue ;
						int l ;
						for ( l = rstart ; l >= 0 && p >= tBuffer ; --p, --l )
							if ( read[l] != *p )
							{
								break ;
							}
						if ( p >= tBuffer )
							continue ;

						struct _overlap no ;
						no.seqIdx = j ;
						no.readStart = l ;
						no.readEnd = len - 1 ;
						no.seqStart = 0 ;
						no.seqEnd = no.readEnd - no.readStart ;
						no.matchCnt = 2 * ( no.readEnd - no.readStart ) ;
						no.similarity = 1.0 ;
						geneOverlap[ geneType ] = no ;
						allOverlaps.push_back( no ) ;
					}
					delete[] tBuffer ;
				}
			}*/
			
			// Test whether the V gene's coordinate is wrong if we have good J,C gene alignment.
			if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 && geneOverlap[3].seqIdx != -1 )
			{
				if ( geneOverlap[2].readEnd + 3 >= geneOverlap[3].readStart && 
					geneOverlap[2].readEnd - 3 <= geneOverlap[3].readStart &&
					( geneOverlap[0].readEnd > geneOverlap[2].readStart + 6 || 
					 geneOverlap[0].readEnd + ( seqs[ geneOverlap[0].seqIdx ].consensusLen - geneOverlap[0].seqEnd - 100 ) 
					 	> geneOverlap[2].readStart + 6 ) ) 
				{
					struct _overlap orig = geneOverlap[0] ;
					geneOverlap[0].seqIdx = -1 ;
					geneOverlap[0].matchCnt = -1 ;
					for ( i = 0 ; i < size ; ++i )
					{
						int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
						if ( geneType != 0 )
							continue ;
						if ( allOverlaps[i].readEnd <= geneOverlap[2].readStart + 6 && 
							allOverlaps[i].readEnd + ( seqs[ allOverlaps[i].seqIdx ].consensusLen 
								- allOverlaps[i].seqEnd - 100 )
							        	 <= geneOverlap[2].readStart + 6 
							&& ( geneOverlap[geneType].seqIdx == -1
								|| IsBetterGeneMatch( allOverlaps[i], geneOverlap[ geneType ], 1.0 ) ) )
						{
							geneOverlap[ geneType ] = allOverlaps[i] ;
						}
					}
					//if ( geneOverlap[0].seqIdx == -1 )
					//	geneOverlap[0] = orig ;
				}
				else if ( geneOverlap[2].readEnd + 3 >= geneOverlap[3].readStart && 
						geneOverlap[2].readEnd - 3 <= geneOverlap[3].readStart &&
						( geneOverlap[0].seqEnd + 100 < seqs[ geneOverlap[0].seqIdx ].consensusLen 
							  && geneOverlap[0].readEnd - geneOverlap[0].readStart + 1 <= 50 ) 
					) 
				{
					geneOverlap[0].seqIdx = -1 ;
				}
			}

			delete[] align ;
			delete[] rvr ;
		}


		// Infer CDR1,2,3.
		char *cdr1 = NULL ;
		char *cdr2 = NULL ;
		char *vAlign = NULL ;
		

		if ( detailLevel >= 2 && geneOverlap[0].seqIdx != -1 
			&& ( geneOverlap[2].seqIdx == -1 || geneOverlap[0].readStart < geneOverlap[2].readStart ) )
		{
			// Infer CDR1, 2
			struct _overlap vgene = geneOverlap[0] ; // Overlap with v-gene
			vAlign = new char[ len + seqs[ vgene.seqIdx ].consensusLen + 2 ] ;
			AlignAlgo::GlobalAlignment( seqs[ vgene.seqIdx ].consensus + vgene.seqStart, vgene.seqEnd - vgene.seqStart + 1,
				read + vgene.readStart, vgene.readEnd - vgene.readStart + 1, vAlign ) ;
			//AlignAlgo::VisualizeAlignment( seqs[ vgene.seqIdx ].consensus + vgene.seqStart, vgene.seqEnd - vgene.seqStart + 1,
			//	read + vgene.readStart, vgene.readEnd - vgene.readStart + 1, vAlign ) ;

			// Locate CDR1.
			int cdrIdx ;
			for ( cdrIdx = 0 ; cdrIdx <= 1 ; ++cdrIdx )
			{
				int seqRangeStart = seqs[ vgene.seqIdx ].info[ cdrIdx ].a ; 
				int seqRangeEnd = seqs[ vgene.seqIdx ].info[ cdrIdx ].b ;
				if ( vgene.seqStart <= seqRangeStart && vgene.seqEnd >= seqRangeEnd )
				{
					i = vgene.readStart - 1 ;
					j = vgene.seqStart - 1 ;
					int readRangeStart, readRangeEnd ;
					int matchCnt = 0 ;
					for ( k = 0 ; vAlign[k] != -1 ; ++k )	
					{
						if ( vAlign[k] != EDIT_DELETE )
							++i ;
						if ( vAlign[k] != EDIT_INSERT )
							++j ;
						
						if ( j == seqRangeStart )
							readRangeStart = i ;
						if ( j >= seqRangeStart && vAlign[k] == EDIT_MATCH )
							matchCnt += 2 ;
						if ( j == seqRangeEnd )
						{
							readRangeEnd = i ;
							break ;
						}
					}
					cdr[cdrIdx].seqIdx = vgene.seqIdx ;
					cdr[cdrIdx].readStart = readRangeStart ;
					cdr[cdrIdx].readEnd = readRangeEnd ;
					cdr[cdrIdx].matchCnt = matchCnt ;
					cdr[cdrIdx].similarity = (double)matchCnt / 
						( readRangeEnd - readRangeStart + 1 + seqRangeEnd - seqRangeStart + 1 ) ;
					//printf( "%d %d: %d %d; %d %d\n", cdrIdx, matchCnt, readRangeStart, readRangeEnd, seqRangeStart, seqRangeEnd ) ;

					/*char *r = new char[readRangeEnd - readRangeStart + 2] ;
					memcpy( r, read + readRangeStart, readRangeEnd - readRangeStart + 1 ) ;
					r[  readRangeEnd - readRangeStart + 1 ] = '\0' ;
					if ( cdrIdx == 0 )
						cdr1 = r ;
					else if ( cdrIdx == 1 )
						cdr2 = r ;*/
				}
			}
		}
		
		char *cdr3 = NULL ;
		double cdr3Score = 0 ;
		
		if ( detailLevel >= 2 )
		{
			// Infer CDR3.
			int s, e ;
			int boundS = 0, boundE = len - 2 ; //[boundS, boundE)
			int range = 37 ;
			bool strongLocateS = false ; // whether the locate method is based on IMGT alignment
			bool strongLocateE = false ;
			if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 )
			{
				// The case that we have anchor.
				// Find the motif for anchor
				if ( geneOverlap[2].readEnd > geneOverlap[0].readEnd )
				{
					int startFrame = geneOverlap[0].seqStart % 3 ;
					int ns = geneOverlap[0].readEnd ; //+ ( seqs[ geneOverlap[0].seqIdx ].consensusLen - 1 - geneOverlap[0].seqEnd ) ;
					s = ns - ( ns - geneOverlap[0].readStart + startFrame ) % 3 ;
					s = s + 6 < len ? s + 6 : s ;
					startFrame = ( seqs[ geneOverlap[2].seqIdx ].consensusLen - 1 - geneOverlap[2].seqEnd ) % 3 ;
					//e = geneOverlap[2].readStart + 
					//	( geneOverlap[2].readEnd - geneOverlap[2].readStart + startFrame ) % 3 ;
					e = geneOverlap[2].readStart ;
					e = e - 6 >= 0 ? e - 6 : e ;
					int adjustE = e ;
					int locate = -1 ;
					for ( i = adjustE ; i < geneOverlap[2].readEnd  && i + 11 < len ; ++i )
					{
						// A strong motif, should be used as the anchor.
						if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
									DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
								&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locate = i ;
							break ;
						}
					}

					if ( locate != -1 )
						e = locate ;

					if ( e < s + 12 )
						range += 15 ;

					if ( s - range > boundS )
						boundS = s - range ;
					if ( e + range < boundE )
						boundE = e + range ;

					if ( locate != -1 )
					{
						s = s + ( e - s ) % 3 ;
						if ( s < e - 18 && 
							geneOverlap[0].seqEnd < seqs[ geneOverlap[0].seqIdx ].consensusLen - 31 )
							s = e - 18 ;
					}
					
					int far = false ;
					for ( i = s ; i <= e ; ++i )
						if ( read[i] == 'M' )
							far = true ;
					if ( far )
					{
						/*s = e - 18 ;
						if ( s < 0 )
							s = e % 3 ;
						for ( i = e ; i >= s ; i -= 3 )
							if ( read[i] == 'M' )
							{
								s = i + 3 ;
								break ;
							}
						if ( s >= e )
						{
							s = 0 ;
							e = len ;
							boundS = 1 ;
						}*/

						if ( seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 
							&& geneOverlap[0].seqEnd < seqs[ geneOverlap[0].seqIdx ].info[2].a )
						{
							s = e - 18 ;
						}
						if ( seqs[ geneOverlap[2].seqIdx ].info[2].a != -1 
							&& geneOverlap[2].seqStart > seqs[ geneOverlap[0].seqIdx ].info[2].a )
						{
							e = s + 18 ;
						}
						
					}
				}
				else
				{
					s = 0 ;
					e = len ;
					boundS = 1 ;
				}	
			
			}
			else if ( geneOverlap[2].seqIdx != -1 )
			{
				e = geneOverlap[2].readStart ;
				e = e - 6 >= 0 ? e - 6 : e ;
				s = e - 12 ;
				if ( s - 31 > boundS )
					boundS = s - 31 ;
				int adjustE = e ;
				int locate = -1 ;
				for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
				{
					if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
								DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
							&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
					{
						locate = i ;
						break ;
					}
				}
				if ( locate != -1 )
				{
					e = locate ;
					s = e - 12 ;
					if ( s < 0 )
						s = 0 ;
				}
				/*s = geneOverlap[2].readStart - 3 + ( e - geneOverlap[2].readStart ) % 3 ;
				if ( s < 0 )
					s = 0 ;
				if ( e + 31 < boundE )
					boundE = e + 31 ;*/
			}
			else if ( geneOverlap[0].seqIdx != -1 
				&& geneOverlap[0].seqEnd >= seqs[ geneOverlap[0].seqIdx ].consensusLen - 50 )
			{
				int startFrame = geneOverlap[0].seqStart % 3 ;
				s = geneOverlap[0].readEnd + ( geneOverlap[0].readEnd - geneOverlap[0].readStart - startFrame ) % 3 ;
				s = s + 6 < len ? s + 6 : s ;
				if ( s >= len ) // the case geneOverlap[0].readEnd is len - 1
					s-= 3 ;
				e = s + 12 ;
				if ( s - 31 > boundS )
					boundS = s - 31 ;
				
				int adjustE = e ;
				int locate = -1 ;
				if ( geneOverlap[3].seqIdx != -1 )
					boundE = geneOverlap[3].readStart - 2 ;
				for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
				{
					if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
								DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
							&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
					{
						locate = i ;
						break ;
					}
				}

				if ( locate != -1 )
				{
					e = locate ;
					s = e - 12 ;
					if ( s < 0 )
						s = 0 ;
				}
				/*else
				{
					e = geneOverlap[0].readEnd + ( geneOverlap[0].readEnd - s ) % 3 ;
					if ( e >= len )
						e = len - 1 ;
				}*/
			}
			else
			{
				s = 0 ;
				e = len ;
				boundS = 1 ;
			}
			
			if ( geneOverlap[2].seqIdx != -1 && boundE > geneOverlap[2].readEnd )
				boundE = geneOverlap[2].readEnd ;
			if ( /*geneOverlap[0].seqIdx == -1 &&*/ s >= boundS )
			{	
				for ( i = s ; i >= boundS ; --i )
					if ( read[i] == 'M' )
					{
						boundS = i + 1 ;
						break ;
					}
			}
			if ( /*geneOverlap[2].seqIdx == -1 &&*/ e <= boundE - 1 )
			{
				for ( i = e ; i < boundE ; ++i )
					if ( read[i] == 'M' )
					{
						boundE = i ;
						break ;
					}
			}
			int locateS = -1 ;
			int locateE = -1 ;
			int extendS = -1 ;
			int extendE = -1 ;

			if ( geneOverlap[0].seqIdx != -1 )
			{
				// Try directly using position 312
				int dest = seqs[ geneOverlap[0].seqIdx ].info[2].a  ;
				if ( dest != -1 /*&& dest <= geneOverlap[0].seqEnd*/ )
				{
					// Left to right scan should be more stable.
					if ( vAlign == NULL )
					{
						struct _overlap vgene = geneOverlap[0] ; // Overlap with v-gene
						vAlign = new char[ len + seqs[vgene.seqIdx].consensusLen + 2 ] ;
						AlignAlgo::GlobalAlignment( seqs[ vgene.seqIdx ].consensus + vgene.seqStart, 
								vgene.seqEnd - vgene.seqStart + 1,
								read + vgene.readStart, vgene.readEnd - vgene.readStart + 1, vAlign ) ;
						//AlignAlgo::VisualizeAlignment( seqs[ vgene.seqIdx ].consensus + vgene.seqStart, 
						//		vgene.seqEnd - vgene.seqStart + 1,
						//		read + vgene.readStart, vgene.readEnd - vgene.readStart + 1, vAlign ) ;
					}
					for ( i = geneOverlap[0].readStart - 1, j = geneOverlap[0].seqStart - 1, k = 0 ; 
						vAlign[k] != -1 ; ++k )	
					{
						if ( vAlign[k] != EDIT_DELETE )
							++i ;
						if ( vAlign[k] != EDIT_INSERT )
							++j ;

						if ( j >= dest )
							break ;
					}
					if ( vAlign[k] == -1 )
					{
						--k ;
						if ( vAlign[k] != EDIT_DELETE )
							--i ;
						if ( vAlign[k] != EDIT_INSERT )
							--j ;
					}
					bool ambiguous = false ;
					for ( int l = k ; l >= 0 && l >= k - 6 ; --l )
						if ( vAlign[l] == EDIT_INSERT || vAlign[l] == EDIT_DELETE )
						{
							ambiguous = true ;
							break ;
						}
					if ( k > 0 && !ambiguous )
					{
						if ( j == dest )
						{
							locateS = i ;
							strongLocateS = true ;
						}
						//else if ( j < dest && j + 6 < dest )
						//	locateS = i + ( dest - j ) ;
						else if ( j < dest )
						{
							extendS = i + dest - j + 5 ; 
						}
					}
				}
			}
			
		
			// The YYC motif on V gene, mentioned in TRUST3 paper, but seems not mentioned in IMGT Junction Analysis.
			if ( locateS == -1 )
			{
				if ( s + 8 >= len )
				{
					if ( ( s - len + 9 ) % 3 )
						s = len - 12 + ( s - len + 9 ) % 3  ;
					else
						s = len - 9 ;
				}
			}

			if ( locateS == -1 )
			{
				for ( i = s ; i >= boundS ; i -= 3 )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
							&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
					{
						locateS = i + 6 ;
						break ;
					}
				}
			}
			

			if ( locateS == -1 )
			{
				// Don't follow the frame rule 
				for ( i = s ; i >= boundS ; --i )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
							&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
					{
						locateS = i + 6 ;
						break ;
					}
				}
			}

			if (locateS == -1 && geneOverlap[0].seqIdx != -1 
					&& seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 )
			{
				// Use the gene sequence information to infer the location
				struct _seqWrapper &seq = seqs[ geneOverlap[0].seqIdx ] ;
				for ( i = s ; i >= boundS ; --i )
				{
					if ( DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' ) 
					{
						int matchLen = 0 ; 
						int geneOffset = AlignAlgo::LocatePartialSufPrefExactMatch( 
								seq.consensus + seq.info[2].a, seq.consensusLen - seq.info[2].a,
								read + i + 6, len - ( i + 6 ),
								locatePartialMatchMinLen, matchLen ) ;
						if ( geneOffset != -1 && geneOffset == 0 )
						{
							locateS = i + 6 ;
							strongLocateS = true ;
							break ;
						}
					}

				}
			}

			if ( locateS == -1 )
			{
				// Try the YxC motif
				for ( i = s ; i >= boundS ; i -= 3 )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							//&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
							&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
					{
						locateS = i + 6 ;
						break ;
					}
				}

				if ( locateS == -1 )
				{
					for ( i = s ; i >= boundS ; --i )
					{
						if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
								//&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
								&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
						{
							locateS = i + 6 ;
							break ;
						}
					}

				}

				if ( locateS == -1 && geneOverlap[0].seqIdx != -1 )
				{
					for ( i = s ; i >= boundS ; --i )
					{
						if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y'
								&& read[i + 6] == 'T' && read[i + 7] == 'G' )
						{
							locateS = i + 6 ; 
							break ;
						}
					}
				}
			}

			if ( locateS == -1 && ( geneOverlap[0].seqIdx != -1 || s <= 18 ) )
			{
				for ( i = s ; i >= boundS ; i -= 3 )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
					{
						locateS = i ;
						break ;
					}
				}
			}
			
			
			if ( locateS == -1 && geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 )
			{
				// Expand the search range.
				int newS = e - 12 ; //- ( e - 12 - s ) % 3 ;
				if ( extendS >= 0 && extendS < newS )
					newS = extendS - ( extendS - s ) % 3 ;
				if ( newS > s )
				{
					for ( i = newS ; i > s ; i -= 3 )
						if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
						{
							locateS = i ;
							break ;
						}
					if ( 0 ) //locateS == -1 )
					{
						for ( i = newS ; i > s ; --i )
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
							{
								locateS = i ;
								break ;
							}
					}

					if ( locateS == -1 )
					{
						for ( i = newS ; i > s ; i -= 3 )
						{
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y'
								&& read[i + 6] == 'T' && read[i + 7] == 'G' )
							{
								locateS = i + 6 ; 
								break ;
							}
						}
					}

				}
			}

			if ( 0 ) //locateS == -1 && geneOverlap[0].seqIdx != -1 )
			{
				for ( i = s ; i >= boundS ; --i )
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
					{
						locateS = i ;
						break ;
					}
			}

			if ( locateS == -1 && geneOverlap[0].seqIdx != -1 )
			{
				// YYx motif.
				for ( i = s ; i >= boundS ; --i )
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' )
					{
						locateS = i + 6 ;
						break ;
					}
			}
			
			int adjustE = e  ;
			if ( 1 ) //locateS != -1 )
			{
				// Use sequence to infer locateE ;
				if ( geneOverlap[2].seqIdx != -1 ) 
				{
					char *align ;

					int dest = seqs[ geneOverlap[2].seqIdx ].info[2].a  ;
					if ( dest != -1 /*&& dest >= geneOverlap[2].seqStart*/ )
					{
						// Left to right scan should be more stable.
						struct _overlap jgene = geneOverlap[2] ; // Overlap with j-gene
						align = new char[ len + seqs[ jgene.seqIdx ].consensusLen + 2 ] ;
						AlignAlgo::GlobalAlignment( seqs[ jgene.seqIdx ].consensus + jgene.seqStart, 
								jgene.seqEnd - jgene.seqStart + 1,
								read + jgene.readStart, jgene.readEnd - jgene.readStart + 1, align ) ;
						//AlignAlgo::VisualizeAlignment( seqs[ jgene.seqIdx ].consensus + jgene.seqStart, 
						//		jgene.seqEnd - jgene.seqStart + 1,
						//		read + jgene.readStart, jgene.readEnd - jgene.readStart + 1, align ) ;
						for ( k = 0 ; align[k] != -1 ; ++k )
							;
						for ( i = jgene.readEnd + 1, j = jgene.seqEnd + 1 ; 
								k >= 0 ; --k )	
						{
							if ( align[k] != EDIT_DELETE )
								--i ;
							if ( align[k] != EDIT_INSERT )
								--j ;

							if ( j <= dest )
								break ;
						}

						bool ambiguous = false ;
						int l = k ;
						if ( k == -1 )
						{
							++l ;
							if ( align[0] != EDIT_DELETE )
								++i ;
							if ( align[0] != EDIT_INSERT )
								++j ;
						}

						for ( ; align[l] != -1 && l <= k + 6 ; ++l )
						{
							if ( align[l] == EDIT_INSERT || align[l] == EDIT_DELETE )
							{
								ambiguous = true ;
								break ;
							}
						}
						if ( !ambiguous )
						{
							if ( j == dest )
							{
								locateE = i ;
								strongLocateE = true ;
							}
							else if ( j == dest + 1 && read[ i - (j - dest)] != 'M' )
							{
								locateE = i - ( j - dest ) ;
							}
							else if ( j > dest && j - dest < 10 )
								extendE = i - ( j - dest ) - 5 ;
						}
						delete[] align ;
					}
				}
				if ( locateS != -1 )
					adjustE = e - ( e - locateS ) % 3 ;
				if ( locateE == -1 )
				{
					for ( i = adjustE ; i < boundE && i + 11 < len ; i += 3 )
					{
						if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
									DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
								&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locateE = i ;
							break ;
						}
					}
				}

				if ( locateE == -1 )
				{
					adjustE = e ;
					for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
					{
						if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
									DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
								&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locateE = i ;
							break ;
						}
					}
				}


				if ( locateE == -1 )
				{
					// Use weaker motif.
					if ( locateS != -1 )
					{
						adjustE = e - ( e - locateS ) % 3 ;
						if ( adjustE + 3 < locateS + 18 )
							adjustE = locateS + 15 ;
					}


					for ( i = adjustE ; i < boundE && i + 11 < len ; i += 3 )
					{
						// The GxG motif
						if ( read[i] == 'T' && DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&&  DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locateE = i ;
							break ;
						}
					}


					if (locateE == -1 && geneOverlap[2].seqIdx != -1 
						&& seqs[ geneOverlap[2].seqIdx ].info[2].a != -1 )
					{
						// Use the gene sequence information to infer the location
						struct _seqWrapper &seq = seqs[ geneOverlap[2].seqIdx ] ;
						for ( i = e ; i < boundE ; ++i )
						{
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' 
								||  DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )
							{
								int matchLen = 0 ; 
								int geneOffset = AlignAlgo::LocatePartialSufSufExactMatch( 
										seq.consensus, seq.info[2].a + 1,
										read, i + 1,
										locatePartialMatchMinLen, matchLen ) ;
								if ( geneOffset != -1 )
								{
									locateE = i ;
									strongLocateE = true ;
									break ;
								}
							}
							
						}
					}

					if ( 0 )//locateE == -1 )
					{
						for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
						{
							// The GxG motif
							if ( read[i] == 'T' && DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
									&&  DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
							{
								locateE = i ;
								break ;
							}
						}

					}

					if ( locateE == -1 )
					{
						for ( i = adjustE ; i < boundE && i + 11 < len ; i += 3 )
						{
							if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
										DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
									&& ( DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
										||  DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' ) )
							{
								locateE = i ;
								break ;
							}
						}
					}

					if ( locateE == -1 && e + 40 > len && boundE == len - 2 )
					{
						for ( i = len - 11 ; i < boundE && i + 5 < len ; ++i )
						{
							if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
										DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
									&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' )
							{
								locateE = i ;
								break ;
							}
						}

						if ( locateE == -1 && geneOverlap[2].seqIdx == -1 )
						{
							for ( i = len - 5 - ( len - 5 - locateS ) % 3 ; i < boundE ; i += 3 )
								if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
										DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' ) 
								{
									locateE = i ;
									break ;
								}
						}
					}

					if ( locateE == -1 ) //&& geneOverlap[2].seqIdx != -1 )
					{
						for ( i = adjustE ; i < boundE ; i += 3 )
						{
							if ( i + 5 < boundE )
								continue ;
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' )
							{
								//printf( "%c%c%c=>%c\n", read[i], read[i + 1], read[i + 2],
								//	DnaToAa( read[i], read[i + 1], read[i + 2] ) ) ;
								//printf( "%c%c%c=>%c\n", read[j], read[j + 1], read[j + 2],
								//	DnaToAa( read[j], read[j + 1], read[j + 2] ) ) ;
								locateE = i ;
								break ;
							}
						}
					}
					if ( locateE == -1 ) //&& geneOverlap[2].seqIdx != -1 )
					{
						for ( i = adjustE ; i < boundE ; i += 3 )
						{
							if ( i + 5 < boundE )
								continue ;
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )
							{
								locateE = i ;
								break ;
							}
						}
					}

					if ( 0 ) //locateE == -1 && geneOverlap[2].seqIdx != -1 )
					{
						// frame shift happens or no locateS.
						adjustE = e ; 
						for ( i = adjustE ; i < boundE ; ++i )
						{
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' )
							{
								//printf( "%c%c%c=>%c\n", read[i], read[i + 1], read[i + 2],
								//	DnaToAa( read[i], read[i + 1], read[i + 2] ) ) ;
								//printf( "%c%c%c=>%c\n", read[j], read[j + 1], read[j + 2],
								//	DnaToAa( read[j], read[j + 1], read[j + 2] ) ) ;
								locateE = i ;
								break ;
							}
						}
						if ( locateE == -1 )
						{
							for ( i = e ; i < boundE ; ++i )
							{
								if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )
								{
									locateE = i ;
									break ;
								}
							}
						}
					}
				}
			}
			if ( locateS != -1 && locateE != -1 )
			{
				if ( locateE + 2 - locateS + 1 < 18 )
				{
					if ( geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 )
						locateS = -1 ;
					else if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx == -1  )
						locateE = -1 ;
				}
				else if ( locateE + 2 - locateS + 1 >= 180 && ( geneOverlap[0].seqIdx == -1 || geneOverlap[2].seqIdx == -1 ) ) 
				{
					locateS = locateE = -1 ;
				}
			}
			
			// If there a gap in the middle of CDR3, pick one side.
			if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 
				&& seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 
				&& seqs[ geneOverlap[2].seqIdx ].info[2].a != -1 
				&& locateS != -1 && locateE != -1 )
			{
				for ( i = locateS ; i <= locateE + 2 ; ++i )
				{
					if ( read[i] == 'M' || read[i] == '\0' )
					{
						if ( strongLocateE 
							&& geneOverlap[0].seqEnd < seqs[ geneOverlap[0].seqIdx ].info[2].a )
							locateS = -1 ;
						if ( strongLocateS 
							&& geneOverlap[2].seqStart > seqs[ geneOverlap[2].seqIdx ].info[2].a )
							locateE = -1 ;
					}
					if ( !read[i] )
						break ;
				}
			}

			// Partial CDR3s
			int sContigIdx = GetContigIdx( locateS, contigs ) ;
			int eContigIdx = GetContigIdx( locateE, contigs ) ;
			bool removeLocateS = false ;
			bool removeLocateE = false ;
			if ( locateS == -1 && locateE != -1 && geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 
				//&& locateE + 11 < contigs[ eContigIdx ].b + 1 
				&& locateE > 15 + contigs[ eContigIdx ].a && locateE <= 60 + contigs[ eContigIdx ].a ) 
			{
				if ( strongLocateE || ( locateE + 11 < len
					&& ( DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'W' ||
					   DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'F' )
					 && DnaToAa( read[locateE + 3], read[ locateE + 4], read[ locateE + 5 ] ) == 'G' 
					 && DnaToAa( read[locateE + 9], read[ locateE + 10], read[ locateE + 11 ] ) == 'G' ) )
				{
					locateS = locateE % 3 ;			
					s = locateS ;
					e = locateE + 2 ;
					
					if ( e - s + 1 >= 18 )
					{
						bool flag = false ;
						for ( i = s ; i <= s + 9 && e - i + 1 >= 18 ; i += 3 )
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' )
							{
								locateS = i ;
								flag = true ;
								break ;
							}
						if ( !flag )
							removeLocateS = true ;
					}
					else
						locateS = -1 ;
				}
			}
			else if ( locateS != -1 && locateE == -1 && geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx == -1 
				//&& locateS - 6 > contigs[ sContigIdx ].a
				&& locateS + 18 < contigs[ sContigIdx ].b + 1 && locateS + 2 + 60 > contigs[ sContigIdx ].b + 1 )
			{
				if ( strongLocateS || ( locateS - 6 >= 0 
					&& DnaToAa( read[locateS], read[ locateS + 1], read[ locateS + 2 ] ) == 'C' 
					&& DnaToAa( read[locateS - 3], read[ locateS - 2], read[ locateS - 1 ] ) == 'Y'  
					&& DnaToAa( read[locateS - 6], read[ locateS - 5], read[ locateS - 4 ] ) == 'Y' ) )
				{
					locateE = ( contigs[ sContigIdx ].b + 1) - 3 - ( ( contigs[ sContigIdx ].b + 1 ) - 3 - locateS ) % 3 ;
				
					s = locateS ;
					e = locateE + 2 ;
						
					if ( e - s + 1 >= 18 )
					{
						bool flag = false ;
						for ( i = e ; i >= e - 9 && i - s + 1 >= 18 ; i -= 3 )
							if ( DnaToAa( read[i - 2], read[i - 1], read[i] ) == 'W' 
								|| DnaToAa( read[i - 2], read[i - 1], read[i] ) == 'F' )
							{
								locateE = i - 2 ;
								flag = true ;
								break ;
							}
						if ( !flag )
							removeLocateE = true ;
					}
					else
						locateE = -1 ;
				}
			}


			// Try to see whether extremely short anchor works using the sequence around CDR anchor.
			// For V gene.
			sContigIdx = GetContigIdx( locateS, contigs ) ; 
			eContigIdx = GetContigIdx( locateE, contigs ) ;
			bool forcePartial = false ; // Let the CDR3 be partial.
			if ( locateS != -1 && locateS <= 18 && geneOverlap[0].seqIdx == -1 )
			{
				// So far, just ignore indel.
				int bestMatchCnt = 0 ;
				int bestHitLen = 0 ;
				int readStart = 0 ;
				SimpleVector<struct _pair> bestTags ;
				int seqCnt = seqs.size() ; 
				for ( i = 0 ; i < seqCnt ; ++i )
				{
					if ( GetGeneType( seqs[i].name ) != 0 || seqs[i].info[2].a == -1 )
						continue ;
					
					struct _seqWrapper &seq = seqs[i] ;
					int matchCnt = 0 ;
					int hitLen = 0 ;

					int geneOffset = seqs[i].info[2].a ; // The position of the gene match with locateS
					// Locate the offset if the match site is not a motif  
					if ( 1 )//DnaToAa( read[ locateS ], read[ locateS + 1], read[ locateS + 2 ] ) != 'C' )
					{
						int matchLen = 0 ;
						geneOffset = AlignAlgo::LocatePartialSufPrefExactMatch( 
								seq.consensus + seq.info[2].a, seq.consensusLen - seq.info[2].a,
								read + locateS, len - locateS,
								locatePartialMatchMinLen, matchLen ) ;
						if ( geneOffset == -1 )
							geneOffset = seq.info[2].a ;
						else
							geneOffset += seq.info[2].a ;
					}

					for ( k = locateS - 1, j = geneOffset - 1 ; k >= 0 && j >=0 ; --k, --j )
					{
						if ( read[k] == 'M' )
							break ;
						if ( seqs[i].consensus[j] == read[k] )
							++matchCnt ;
						++hitLen ;
					}
					int tmp = k + 1 ;

					// The other end, just extend to a mismatch point.
					for ( k = locateS, j = geneOffset ; k < len && j < seq.consensusLen ; ++k, ++j )
					{
						if ( seq.consensus[j] != read[k] )
							break ;
						++matchCnt ;
						++hitLen ;
					}
					
					if ( matchCnt > bestMatchCnt )
					{
						bestMatchCnt = matchCnt ;
						bestHitLen = hitLen ; 
						bestTags.Clear() ;
						struct _pair np ;
						np.a = i ;
						np.b = geneOffset ;
						readStart = tmp ;
						bestTags.PushBack( np ) ;
					}
					else if ( matchCnt == bestMatchCnt )
					{
						struct _pair np ;
						np.a = i ;
						np.b = geneOffset ;
						bestTags.PushBack( np ) ;
					}
				}
				
				int anchorSeqIdx = -1 ;
				int anchorType = -1 ;
				if ( geneOverlap[2].seqIdx != -1 )
				{
					anchorSeqIdx = geneOverlap[2].seqIdx ;
					anchorType = 2 ;
				}
				else if ( geneOverlap[3].seqIdx != -1 )
				{
					anchorSeqIdx = geneOverlap[3].seqIdx ;	
					anchorType = 3 ;
				}
				//printf( "%d %d %d\n", bestMatchCnt, bestHitLen, bestTags.Size() ) ;
				if ( bestHitLen > 9 && bestMatchCnt / (double)bestHitLen >= 0.9 )
				{
					int size = bestTags.Size() ;
					bool start = false ;
					for ( i = 0 ; i < size ; ++i )
					{
						struct _overlap no ;
						no.seqIdx = bestTags[i].a ;
						no.readStart = readStart ;
						no.readEnd = readStart + bestHitLen - 1 ;
						no.seqStart = bestTags[i].b - ( locateS - readStart ) ;
						no.seqEnd = no.seqStart + bestHitLen - 1 ;
						no.matchCnt = 2 * bestMatchCnt ;
						no.similarity = bestMatchCnt / (double)bestHitLen ;
						if ( anchorSeqIdx != -1 )
						{
							if ( no.readEnd > geneOverlap[ anchorType ].readStart 
								|| !IsSameChainType( seqs[ no.seqIdx ].name, seqs[ anchorSeqIdx ].name ) )
								continue ;
						}

						if ( !start )
						{
							geneOverlap[0] = no ;
							if ( seqs[ bestTags[i].a ].info[2].a != bestTags[i].b ) 
							{
								// adjust locateS, if we can match the point
								int diff = bestTags[i].b - seqs[ bestTags[i].a ].info[2].a ; 
								if ( locateS - diff >= no.readStart && locateS + diff <= no.readEnd )
								{
									locateS -= diff ;
									removeLocateS = false ;
								}
							}
							if ( removeLocateS &&  // The CDR3 is likely to be partial
								seqs[ bestTags[i].a ].info[2].a != bestTags[i].b )
								forcePartial = true ;
							start = true ;
						}
						allOverlaps.push_back( no ) ;
					}
					removeLocateS = false ;
				}
			}
			
			// For J gene
			if ( locateE != - 1 )
			{
				int distToEnd = contigs[ eContigIdx ].b - locateE ;
				if ( distToEnd <= 18 && geneOverlap[2].seqIdx == -1 )
				{
					// So far, just ignore indel.
					int bestMatchCnt = 0 ;
					SimpleVector<struct _pair> bestTags ;
					int seqCnt = seqs.size() ;
					int bestHitLen = 0 ;
					int readEnd = 0 ;
					for ( i = 0 ; i < seqCnt ; ++i )
					{
						if ( GetGeneType( seqs[i].name ) != 2 || seqs[i].info[2].a == -1 )
							continue ;
						struct _seqWrapper &seq = seqs[i] ;

						int geneOffset = seq.info[2].a ; // The coordinate of the gene match with locateE 
						int matchCnt = 0 ;
						int hitLen = 0 ;
						if ( locateE < len ) /*locateE + 2 >= len || 
							( DnaToAa( read[ locateE ], read[ locateE + 1], read[ locateE + 2 ] ) != 'F' 
							  && DnaToAa( read[ locateE ], read[ locateE + 1], read[ locateE + 2 ] ) != 'W' ) )*/ 
						{
							int matchLen = 0 ;
							geneOffset = AlignAlgo::LocatePartialSufSufExactMatch( 
								seq.consensus, seq.info[2].a + 1,
								read, locateE + 1,
								locatePartialMatchMinLen, matchLen ) ;
							if ( geneOffset == -1 )
								geneOffset = seq.info[2].a ;
							else
								geneOffset += matchLen - 1 ; 
							/*for ( ; geneOffset >= 8 ; --geneOffset )
							{
								cnt = 0 ;
								for ( k = locateE, j = geneOffset ; k >= 0 && j >= 0 ; --k, --j )
								{
									if ( seq.consensus[j] != read[k] )
										break ;
									++cnt ;
								}
								if ( cnt > 9 )
									break ;
							}
							if ( cnt <= 9 )
							{
								geneOffset = seq.info[2].a ;
							}*/
						}
						
						for ( k = locateE + 1, j = geneOffset + 1 ; k < len && j < seqs[i].consensusLen ; ++k, ++j )
						{
							if ( read[k] == 'M' )
								break ;
							if ( seq.consensus[j] == read[k] )
								++matchCnt ;
							++hitLen ;
						}
						int tmp = k - 1 ;
						for ( k = locateE, j = geneOffset ; k >= 0 && j >= 0 ; --k, --j )
						{
							if ( seq.consensus[j] != read[k] )
								break ;
							++matchCnt ;
							++hitLen ;
						}
						
						struct _pair np ;
						if ( matchCnt > bestMatchCnt )
						{
							bestMatchCnt = matchCnt ;
							bestHitLen = hitLen ;
							bestTags.Clear() ;
							readEnd = tmp ;
							np.a = i ;
							np.b = geneOffset ;
							bestTags.PushBack( np ) ;
						}
						else if ( matchCnt == bestMatchCnt )
						{
							np.a = i ;
							np.b = geneOffset ;
							bestTags.PushBack( np ) ;
						}
					}

					int anchorSeqIdx = -1 ;
					int anchorType = -1 ;
					if ( geneOverlap[0].seqIdx != -1 )
					{
						anchorSeqIdx = geneOverlap[0].seqIdx ;
						anchorType = 0 ;
					}
					else if ( geneOverlap[3].seqIdx != -1 )
					{
						anchorSeqIdx = geneOverlap[3].seqIdx ;	
						anchorType = 3 ;
					}
					if ( bestHitLen > 9 && bestMatchCnt / (double)bestHitLen >= 0.9 )
					{
						int size = bestTags.Size() ;
						bool start = false ;
						for ( i = 0 ; i < size ; ++i )
						{
							struct _overlap no ;
							no.seqIdx = bestTags[i].a ;
							no.readStart = readEnd - bestHitLen + 1 ;
							no.readEnd = readEnd ;
							no.seqEnd = bestTags[i].b + ( readEnd - locateE ) ;
							no.seqStart = no.seqEnd - bestHitLen + 1 ;
							no.matchCnt = 2 * bestMatchCnt ;
							no.similarity = bestMatchCnt / (double)bestHitLen ;
							if ( anchorSeqIdx != -1 )
							{
								if ( ( anchorType == 0 && no.readStart < geneOverlap[ anchorType ].readEnd ) 
										|| !IsSameChainType( seqs[ no.seqIdx ].name, 
												seqs[ anchorSeqIdx ].name ) )
									continue ;
							}
							if ( !start )
							{
								geneOverlap[2] = no ;

								if ( seqs[ bestTags[i].a ].info[2].a != bestTags[i].b ) 
								{
									// adjust locateS, if we can match the point
									int diff = bestTags[i].b - seqs[ bestTags[i].a ].info[2].a ; 
									if ( locateE - diff >= no.readStart && locateE + diff <= no.readEnd )
									{
										locateE -= diff ;
										removeLocateE = false ;
									}
								}

								if ( removeLocateE &&  // The CDR3 is likely to be partial
										seqs[ bestTags[i].a ].info[2].a != bestTags[i].b )
									forcePartial = true ;
								start = true ;
							}
							allOverlaps.push_back( no ) ;
						}
						removeLocateE = false ;
					}
				}
			}

			if ( removeLocateS )
				locateS = -1 ;
			if ( removeLocateE )
				locateE = -1 ;
			
			sContigIdx = GetContigIdx( locateS, contigs ) ; 
			eContigIdx = GetContigIdx( locateE, contigs ) ;
			if ( locateS != -1 && locateE != -1 && locateE + 2 - locateS + 1 >= 18 )
			{
				s = locateS ;
				e = locateE + 2 ;
				
				/*cdr3 = new char[e - s + 2  + 1 ] ;
				memcpy( cdr3, read + s, e - s + 1 ) ;
				cdr3[e - s + 1] = '\0' ;*/

				cdr[2].seqIdx = 0 ;
				cdr[2].readStart = s ;
				cdr[2].readEnd = e ;

				int leftCnt = 0 ;
				int rightCnt = 0 ;

				// Use the anchor motif to score the cdr
				if ( geneOverlap[0].seqIdx != -1 && seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 )
				{
					// Because the YYC motif is not always true on V gene, 
					// use the sequence from reference to test
					char *ref = seqs[ geneOverlap[0].seqIdx ].consensus ;
					int offset = seqs[ geneOverlap[0].seqIdx ].info[2].a ;
					if ( locateS - 6 > 0 )
						if ( DnaToAa( read[locateS - 6], read[ locateS - 5], read[ locateS - 4 ] ) 
							== DnaToAa( ref[offset - 6], ref[ offset - 5], ref[ offset - 4 ] )  )
						{
							cdr3Score += 100.0 / 6 ;
							++leftCnt ;
						}
					if ( locateS - 3 > 0 )
						if ( DnaToAa( read[locateS - 3], read[ locateS - 2], read[ locateS - 1 ] ) 
							== DnaToAa( ref[offset - 3], ref[ offset - 2], ref[ offset - 1 ] )  )
						{
							cdr3Score += 100.0 / 6 ;
							++leftCnt ;
						}
					if ( DnaToAa( read[locateS], read[ locateS + 1], read[ locateS + 2] ) 
							== DnaToAa( ref[offset], ref[ offset + 1], ref[ offset + 2 ] )  )
					{
						cdr3Score += 100.0 / 6 ;
						++leftCnt ;
					}
						
				}
				else
				{
					if ( locateS - 6 > 0 )
						if ( DnaToAa( read[locateS - 6], read[ locateS - 5], read[ locateS - 4 ] ) == 'Y' )
						{
							cdr3Score += 100.0 / 6 ;
							++leftCnt ;
						}
					if ( locateS - 3 > 0 )
						if ( DnaToAa( read[locateS - 3], read[ locateS - 2], read[ locateS - 1 ] ) == 'Y' )
						{
							cdr3Score += 100.0 / 6 ;
							++leftCnt ;
						}

					if ( DnaToAa( read[locateS], read[ locateS + 1], read[ locateS + 2 ] ) == 'C' )
					{
						cdr3Score += 100.0 / 6 ;
						++leftCnt ;
					}
				}
				if ( locateE + 2 < len && ( DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'W' ||
						DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'F' ) )
				{
					cdr3Score += 100.0 / 6 ;
					++rightCnt ;
				}

				if ( locateE + 5 < len )
					if ( DnaToAa( read[locateE + 3], read[ locateE + 4], read[ locateE + 5 ] ) == 'G' )
					{
						cdr3Score += 100.0 / 6 ;
						++rightCnt ;
					}
				if ( locateE + 11 < len )
					if ( DnaToAa( read[locateE + 9], read[ locateE + 10], read[ locateE + 11 ] ) == 'G' )
					{
						cdr3Score += 100.0 / 6 ;
						++rightCnt ;
					}

				if ( s < 0 )  
				{
					// The CDR3 anchor based on alignment of V,J genes could create some boundary cases when indel. 
					//    Since these come from alignment, it should happend to the case of missing V/J gene.
					s = e % 3 ;
					cdr[2].readStart = s ;
					cdr3Score = 0 ;
				}
				if ( e >= len )
				{
					e = len - 1 - ( len - s ) % 3 ; 
					cdr[2].readEnd = e ;
					cdr3Score = 0 ;
				}
				
				if ( cdr3Score < 99 && ( ( leftCnt < 3 && geneOverlap[0].seqIdx == -1 ) 
						|| ( rightCnt < 3 && geneOverlap[2].seqIdx == -1 ) ) )
					cdr3Score = 0 ;
				else if ( e + 6 >= len && locateE + 2 < len &&
						!( DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'W' ||
							DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'F' ) )
					cdr3Score = 0 ;
				else if ( cdr3Score < 99 && 
					( geneOverlap[0].seqIdx != -1 && geneOverlap[0].seqStart > 100 && geneOverlap[0].readStart > 100 ) 
					&& ( !strongLocateS || leftCnt < 3 ) )
					cdr3Score = 0 ;
				else if ( cdr3Score < 99 &&
						geneOverlap[0].seqIdx != -1 && ( !strongLocateS || leftCnt < 3 ) &&
						( GetContigIdx( geneOverlap[0].readEnd, contigs ) == GetContigIdx(s, contigs)) && 
						( ( seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 
							&& geneOverlap[0].seqEnd + ( s - geneOverlap[0].readEnd ) + 5 < 
								seqs[ geneOverlap[0].seqIdx ].info[2].a ) 
						|| ( seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 
							&& geneOverlap[0].seqEnd + ( s - geneOverlap[0].readEnd ) + 51 < 
								seqs[ geneOverlap[0].seqIdx ].consensusLen ) ) )
					cdr3Score = 0 ;
				else if ( cdr3Score < 99 &&
						geneOverlap[2].seqIdx != -1 && ( !strongLocateE || rightCnt < 3 ) &&
						( GetContigIdx( geneOverlap[2].readStart, contigs ) == GetContigIdx(e, contigs)) && 
						( seqs[ geneOverlap[2].seqIdx ].info[2].a != -1 
							&& geneOverlap[2].seqStart + ( (e - 2) - geneOverlap[2].readStart ) - 5 > 
								seqs[ geneOverlap[2].seqIdx ].info[2].a ) ) 
					cdr3Score = 0 ;
				else if ( geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 && s >= geneOverlap[2].readStart )
					cdr3Score = 0 ;
				else if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx == -1 && e <= geneOverlap[0].readEnd )
					cdr3Score = 0 ;
				else if ( geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 )
				{
					// Check the contigs containing s.
					for ( i = 0 ; i < contigCnt ; ++i )
						if ( s <= contigs[i].b )
							break ;
					if ( i >= contigCnt || s - 50 >= contigs[i].a )
						cdr3Score = 0 ;
				}
				else if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx == -1 )
				{
					for ( i = contigCnt - 1 ; i >= 0 ; --i )	
						if ( e >= contigs[i].a )
							break ;
					if ( i < 0 || ( e + 50 <= contigs[i].b && rightCnt < 3 ) )
						cdr3Score = 0 ;
				}
				else if ( forcePartial )
					cdr3Score = 0 ;
				
				// Now consider whether the gaps could create some false positive score CDR3.
				if ( cdr3Score > 0 )
				{
					int nCnt = 0 ;
					// Gap on the boundary, then we need to adjust the CDR3 region.
					//   this could happen when doing alignment, it does not check the overhang nucleotide.
					if ( read[s] == 'M' )
					{
						while ( read[s] == 'M' && s <= e )
							s += 3 ;
						cdr[2].readStart = s ;
						cdr3Score = 0 ;
						if ( s >= e )
						{
							cdr[2].seqIdx = -1 ;
							cdr[2].readStart = cdr[2].readEnd = -1 ;
						}
					}
					if ( read[e] == 'M' )
					{
						while ( read[e] == 'M' && e >= s )
							e -= 3 ;
						cdr[2].readEnd = e ;
						cdr3Score = 0 ;
						if ( e <= s )
						{
							cdr[2].seqIdx = -1 ;
							cdr[2].readStart = cdr[2].readEnd = -1 ;
						}
					}
					// Gap in the middle.
					for ( i = s ; i <= e ; ++i )
					{
						if ( read[i] == 'N' )
						{
							++nCnt ;
							if ( nCnt >= 2 )
							{
								cdr3Score = 0 ;
								break ;
							}
						}
						else if ( read[i] == 'M' )
						{
							cdr3Score = 0 ;
							break ;
						}
					}
				}

				// If the V, J gene already contain the anchor and but the anchor motif is in the gap then we 
				//    need to force the CDR3 to be partial.
				if ( geneOverlap[0].seqIdx != -1 && seqs[ geneOverlap[0].seqIdx ].info[2].a != -1 )
				{
					if ( geneOverlap[0].seqEnd >= seqs[ geneOverlap[0].seqIdx ].info[2].a + 2 
						&& s > geneOverlap[0].readEnd )
						cdr3Score = 0 ;
				}
				if ( geneOverlap[2].seqIdx != -1 && seqs[ geneOverlap[2].seqIdx ].info[2].a != -1 )
				{
					if ( geneOverlap[2].seqStart <= seqs[ geneOverlap[2].seqIdx ].info[2].a  
						&& e < geneOverlap[2].readStart )
						cdr3Score = 0 ;
				}


				if ( cdr3Score < 100 )
				{
					for ( i = 1 ; i < contigCnt ; ++i )
					{
						/*if ( s >= contigs[i].a && s < contigs[i].a + 18 
							&& geneOverlap[0].seqIdx != -1 && geneOverlap[0].readEnd <= contigs[i - 1].b 
							&& DnaToAa( read[s], read[s + 1], read[s + 2] ) != 'C' )
						{
							cdr3Score = 0 ;
							break ;
						}
						else if ( s >= contigs[i].a && s >= contigs[i].a + 18 && s <= contigs[i].b 
							&& geneOverlap[0].seqIdx != -1 && geneOverlap[0].readEnd <= contigs[i - 1].b 
							&& leftCnt < 3 )
						{
							cdr3Score = 0 ;
							break ;
						}*/

						if ( s >= contigs[i].a && s <= contigs[i].b )
						{
							if ( geneOverlap[0].seqIdx != -1 && geneOverlap[0].readEnd <= contigs[i - 1].b 
								&& leftCnt < 3 && !strongLocateS )
							{
								// Check whether it matches the sequence of the V gene.
								int matchCnt = 0 ;
								int hitLen = 0 ;
								int seqIdx = geneOverlap[0].seqIdx ;
								if ( seqs[ seqIdx ].info[2].a != -1 )
								{
									int rightMatch = 0 ;
									for ( j = s , k = seqs[ seqIdx ].info[2].a ; 
											read[j] && seqs[ seqIdx ].consensus[k] ; ++j, ++k )
									{
										//if ( read[j] != seqs[seqIdx].consensus[k] )
										//	break ;

										if ( read[j] == seqs[ seqIdx ].consensus[k] )
										{
											++rightMatch ;
											if ( (double)rightMatch / ( j - s + 1 ) >= 0.75 )
											{
												matchCnt = rightMatch ;
												hitLen = j - s + 1 ;
											}
										}
									}

									for ( j = s - 1, k = seqs[ seqIdx ].info[2].a - 1; 
											j >= 0 && k >= 0 ; --j, --k ) 
									{
										if ( read[j] == 'M' )
											break ;

										if ( read[j] == seqs[ seqIdx ].consensus[k] )
											++matchCnt ;
										++hitLen ;
									}
								}
								// hitLen by default is 0.
								double similarity = 0.9 ;
								if ( seqs[ seqIdx ].name[0] == 'I' )
									similarity = 0.8 ;
								if ( hitLen <= 9 || (double)matchCnt / hitLen < similarity )
									cdr3Score = 0 ;
								break ;
							}
							break ;
						}
					}
					for ( i = contigCnt - 2 ; i > 0 ; --i )
					{
						/*if ( e <= contigs[i].b && e > contigs[i].b - 18 
								&& geneOverlap[2].seqIdx != -1 && geneOverlap[2].readStart >= contigs[i + 1].a 
								&& DnaToAa( read[e - 2], read[e - 1], read[e] ) != 'W' &&
								  DnaToAa( read[e - 2], read[e - 1], read[e] ) != 'F' )
						{
							cdr3Score = 0 ;
							break ;
						}
						else if ( e <= contigs[i].b && e > contigs[i].b - 18 && e >= contigs[i].a  
								&& geneOverlap[2].seqIdx != -1 && geneOverlap[2].readStart >= contigs[i + 1].a
								&& rightCnt < 3 ) 
						{
							cdr3Score = 0 ;
							break ;
						}*/
						if ( e >= contigs[i].a && e <= contigs[i].b )
						{
							if ( geneOverlap[2].seqIdx != -1 && geneOverlap[2].readStart >= contigs[i + 1].a 
								&& rightCnt < 3 && !strongLocateE )
							{
								int matchCnt = 0 ;
								int hitLen = 0 ;
								int seqIdx = geneOverlap[2].seqIdx ;
								if ( seqs[ seqIdx ].info[2].a != -1 )
								{
									int leftMatch = 0 ;
									for ( j = e, k = seqs[ seqIdx ].info[2].a + 2 ; j >= 0 && k >= 0 ; --j, --k ) 
									{
										//if ( read[j] != seqs[seqIdx].consensus[k] )
										//	break ;

										if ( read[j] == seqs[ seqIdx ].consensus[k] )
										{
											++leftMatch ;
											if ( (double)leftMatch / ( e - j + 1 ) >= 0.75 )
											{
												matchCnt = leftMatch ;
												hitLen = e - j + 1 ;
											}
										}
									}
									for ( j = e + 1, k = seqs[ seqIdx ].info[2].a + 3 ; 
										read[j] && seqs[ seqIdx ].consensus[k] ; ++j, ++k )
									{
										if ( read[j] == 'M' )
											break ;

										if ( read[j] == seqs[ seqIdx ].consensus[k] )
											++matchCnt ;
										++hitLen ;

									}
								}
								double similarity = 0.9 ;
								if ( seqs[ seqIdx ].name[0] == 'I' )
									similarity = 0.8 ;
								if ( hitLen <= 9 || (double)matchCnt / hitLen < similarity )
									cdr3Score = 0 ;
								break ;
							}
							break ;
						}
					}
					
				}
			}
			// Partial CDR3s
			else if ( locateS == -1 && locateE != -1 && geneOverlap[2].seqIdx != -1 
				&& ( geneOverlap[0].seqIdx == -1 
					|| GetContigIdx( geneOverlap[0].readStart, contigs ) != GetContigIdx( geneOverlap[2].readStart, contigs  ) )
				&& locateE > 15 + contigs[ eContigIdx ].a && locateE <= 60 + contigs[ eContigIdx ].a ) 
			{
				if ( strongLocateE || ( locateE + 11 < len  
					 && ( DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'W' ||
					    DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'F' )
					 && DnaToAa( read[locateE + 3], read[ locateE + 4], read[ locateE + 5 ] ) == 'G' 
					 && DnaToAa( read[locateE + 9], read[ locateE + 10], read[ locateE + 11 ] ) == 'G' ) )
				{
					while ( locateE + 2 >= len )
						locateE -= 3 ;
					locateS = contigs[ eContigIdx ].a + ( locateE - contigs[ eContigIdx ].a ) % 3 ;			
					cdr3Score = 0 ;

					s = locateS ;
					e = locateE + 2 ;
					
					if ( e - s + 1 >= 18 )
					{
						for ( i = s ; i <= s + 9 && e - i + 1 >= 18 ; i += 3 )
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' )
							{
								s = i ;
								break ;
							}
						
						/*cdr3 = new char[e - s + 2  + 1 ] ;
						memcpy( cdr3, read + s, e - s + 1 ) ;
						cdr3[e - s + 1] = '\0' ;*/
						if ( s + 2 < geneOverlap[2].readStart )
						{
							cdr[2].seqIdx = 0 ;
							cdr[2].readStart = s ;
							cdr[2].readEnd = e ;
						}
					}
				}
			}
			else if ( locateS != -1 && locateE == -1 && geneOverlap[0].seqIdx != -1 //&& geneOverlap[2].seqIdx == -1
				&& ( geneOverlap[2].seqIdx == -1 || 
					GetContigIdx( geneOverlap[0].readStart, contigs ) != GetContigIdx( geneOverlap[2].readStart, contigs ) )  
				&& locateS + 18 < contigs[ sContigIdx ].b + 1 && locateS + 2 + 60 > contigs[ sContigIdx ].b + 1 )
			{
				if ( strongLocateS || ( locateS - 6 >= 0  
					&& DnaToAa( read[locateS], read[ locateS + 1], read[ locateS + 2 ] ) == 'C' 
					&& DnaToAa( read[locateS - 3], read[ locateS - 2], read[ locateS - 1 ] ) == 'Y'  
					&& DnaToAa( read[locateS - 6], read[ locateS - 5], read[ locateS - 4 ] ) == 'Y' ) )
				{
					while ( locateS < 0 )
						locateS += 3 ;

					locateE = contigs[ sContigIdx ].b - 2 - ( contigs[ sContigIdx ].b - 2 - locateS ) % 3 ;
					cdr3Score = 0 ;
				
					s = locateS ;
					e = locateE + 2 ;
						
					if ( e - s + 1 >= 18 )
					{
						for ( i = e ; i >= e - 9 && i - s + 1 >= 18 ; i -= 3 )
							if ( DnaToAa( read[i - 2], read[i - 1], read[i] ) == 'W' 
								|| DnaToAa( read[i - 2], read[i - 1], read[i] ) == 'F' )
							{
								e = i ;
								break ;
							}

						/*cdr3 = new char[e - s + 2  + 1 ] ;
						memcpy( cdr3, read + s, e - s + 1 ) ;
						cdr3[e - s + 1] = '\0' ;*/
						if ( e - 2 > geneOverlap[0].readEnd )
						{
							cdr[2].seqIdx = 0 ;
							cdr[2].readStart = s ;
							cdr[2].readEnd = e ;
						}
					}
				}
			}

			cdr[2].similarity = cdr3Score / 100.0 ;
				
		}
		
		if ( detailLevel >= 2 && cdr[2].similarity > 0 )
		{
			AnnotateReadDGene( read, geneOverlap, cdr, secondaryGeneOverlaps ) ;
		}
		if ( vAlign != NULL )
			delete[] vAlign ;

		if ( detailLevel >= 2 )
		{
			// Change gap nucleotide back from 'M' to 'N'
			for ( i = 0 ; i < contigCnt - 1 ; ++i )
			{
				//printf( "%d %d %d %d\n", contigs[i].a, contigs[i].b, contigs[i + 1].a, contigs[i + 1].b ) ;
				for ( j = contigs[i].b + 1 ; j < contigs[i + 1].a ; ++j )
					read[j] = 'N' ;
			}

		}
		// Compute the name
		if ( secondaryGeneOverlaps != NULL )
		{
			for ( i = 0 ; i < 4 ; ++i )
			{	
				if ( i == 1 ) // skip the D gene.
					continue ;
				if ( geneOverlap[i].seqIdx != -1 )
				{
					//int offset = strlen( buffer ) ;
					int seqIdx = geneOverlap[i].seqIdx ;
					/*sprintf( buffer + offset, " %s(%d):(%d-%d):(%d-%d):%.2lf",
					  seqs[ seqIdx ].name, seqs[ seqIdx ].consensusLen,
					  geneOverlap[i].readStart, geneOverlap[i].readEnd, 
					  geneOverlap[i].seqStart, geneOverlap[i].seqEnd, geneOverlap[i].similarity * 100 ) ;*/

					// Output the secondary assignment	
					int size = allOverlaps.size() ;
					int reportCnt = 0 ;
					SimpleVector<int> usedSeqIdx ;
					usedSeqIdx.Reserve( 5 ) ;
					for ( j = 0 ; j < size ; ++j )
					{
						if ( GetGeneType( seqs[ allOverlaps[j].seqIdx ].name ) != i )
							continue ;
						int l ;
						int seqIdx2 = allOverlaps[j].seqIdx ;
						if ( seqIdx2 == seqIdx ||
								!IsBetterGeneMatch( allOverlaps[j], geneOverlap[i], 0.95 )
								/*|| ( allOverlaps[j].readStart < geneOverlap[i].readStart - 20 
								  || allOverlaps[j].readStart > geneOverlap[i].readStart +  20 ) 
								  || ( allOverlaps[j].readEnd < geneOverlap[i].readEnd - 20 
								  || allOverlaps[j].readEnd > geneOverlap[i].readEnd +  20 )  */
						   )
							continue ;

						if ( IsSameGeneAllele( seqs[ seqIdx ].name, seqs[ seqIdx2 ].name ) )
							continue ;
						for ( l = 0 ; l < usedSeqIdx.Size() ; ++l )
						{
							if ( IsSameGeneAllele( seqs[ usedSeqIdx[l] ].name, seqs[ seqIdx2 ].name ) )
								break ;
						}
						if ( l < usedSeqIdx.Size() )
							continue ;

						++reportCnt ;
						/*sprintf( buffer + strlen( buffer ), ",%s(%d):(%d-%d):(%d-%d):%.2lf",
						  seqs[ seqIdx2 ].name, seqs[ seqIdx2 ].consensusLen,
						  allOverlaps[j].readStart, allOverlaps[j].readEnd, 
						  allOverlaps[j].seqStart, allOverlaps[j].seqEnd, allOverlaps[j].similarity * 100 ) ; */
						secondaryGeneOverlaps->push_back( allOverlaps[j] ) ;
						usedSeqIdx.PushBack( allOverlaps[j].seqIdx ) ;
						if ( reportCnt >= 2 )
							break ;
					}
				}
				/*else
				  {
				  sprintf( buffer + strlen( buffer ), " *" ) ;
				  }*/
			}
		}
		
		/*if ( cdr1 == NULL)
			sprintf( buffer + strlen( buffer), " CDR1(0-0):0.00=null" ) ;
		else
			sprintf( buffer + strlen( buffer), " CDR1(%d-%d):%.2lf=%s", cdr[0].readStart, cdr[0].readEnd, 
				cdr[0].similarity * 100, cdr1 ) ;
		
		if ( cdr2 == NULL)
			sprintf( buffer + strlen( buffer), " CDR2(0-0):0.00=null" ) ;
		else
			sprintf( buffer + strlen( buffer), " CDR2(%d-%d):%.2lf=%s", cdr[1].readStart, cdr[1].readEnd, 
				cdr[1].similarity * 100, cdr2 ) ;

		if ( cdr3 == NULL)
			sprintf( buffer + strlen( buffer), " CDR3(0-0):0.00=null" ) ;
		else
			sprintf( buffer + strlen( buffer), " CDR3(%d-%d):%.2lf=%s", cdr[2].readStart, cdr[2].readEnd, cdr3Score, cdr3 ) ;*/
		
		if ( cdr1 != NULL )
			delete[] cdr1 ;
		if ( cdr2 != NULL )
			delete[] cdr2 ;
		if ( cdr3 != NULL )
			delete[] cdr3 ;
		return 1 ;
	}
	
	//  How similar the CDR3 is to the germline sequence.
	double GetCDR3Similarity( char *seq, struct _overlap geneOverlap[4], struct _overlap cdr[3] )
	{
		int i ;
		if ( cdr[2].similarity <= 0 ) // partial CDR3.
			return 0 ;
		if ( geneOverlap[0].seqIdx == -1 || geneOverlap[2].seqIdx == -1 )
			return 0 ;
		int seqIdx = geneOverlap[0].seqIdx ;
		int hasD = 0 ;
		
		if ( geneOverlap[0].readEnd < cdr[2].readStart || geneOverlap[0].readStart > cdr[2].readStart )
			return 0 ;
		if ( geneOverlap[2].readStart > cdr[2].readEnd || geneOverlap[2].readEnd < cdr[2].readEnd )
			return 0 ;
		
		if ( seqs[seqIdx].name[2] == 'H' || seqs[seqIdx].name[2] == 'B' 
			|| seqs[seqIdx].name[2] == 'D' )
		{
			// Contains D gene
			if ( geneOverlap[1].seqIdx == -1 )
				return 0 ;
			hasD = 1 ;
		}
		
		int matchCnt = 0, mismatchCnt = 0, indelCnt = 0 ;
		int len = strlen( seq ) ;
		char *align ;
		align = new char[ 3 * len + 102] ;
		len = 0 ;
		int seqStart, seqEnd, readStart, readEnd ;
		for ( i = 0 ; i < 3 ; ++i )
		{
			if ( hasD == 0 && i == 1 )
				continue ;
			struct _overlap gene = geneOverlap[i] ;
			if ( i == 0 )
			{
				readStart = cdr[2].readStart ;
				readEnd = gene.readEnd ;
				
				if ( seqs[gene.seqIdx].info[2].a != -1 )
					seqStart = seqs[gene.seqIdx].info[2].a ;
				else
					seqStart = gene.seqEnd - (readEnd - readStart) ;
				seqEnd = gene.seqEnd ; 
			}
			else if ( i == 1 )
			{
				readStart = gene.readStart ;
				readEnd = gene.readEnd ;
				seqStart = gene.seqStart ;
				seqEnd = gene.seqEnd ;
			}
			else if ( i == 2 )
			{
				readStart = gene.readStart ;
				readEnd = cdr[2].readEnd ;

				seqStart = gene.seqStart ;
				if ( seqs[gene.seqIdx].info[2].a != -1 )
					seqEnd = seqs[gene.seqIdx].info[2].a ;
				else
					seqEnd = gene.seqStart + (readEnd - readStart) ;
			}

			if ( readEnd - readStart < 0 || seqEnd - seqStart < 0 )
			{
				// Annotation coordinate is wrong. Force the similarity to be 0.
				matchCnt = 0 ;
				break ;
			}
			
			AlignAlgo::GlobalAlignment( seqs[ gene.seqIdx ].consensus + seqStart,
					seqEnd - seqStart + 1,
					seq + readStart - cdr[2].readStart, readEnd - readStart + 1, align ) ;

			GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;
			len += ( seqEnd - seqStart + 1 ) ;
		}

		delete[] align ;
		if ( len == 0 )
			return 0 ;
		else
			return (double)matchCnt / len ;
	}

	int GetEqualSecondaryGeneOverlap( struct _overlap &primaryOverlap, int geneType, std::vector<struct _overlap> *secondaryGeneOverlaps, 
		std::vector<int> &secondaryOverlapIdx )
	{
		int i ;
		int seqIdx = primaryOverlap.seqIdx ;
		if ( seqIdx == -1 || secondaryGeneOverlaps == NULL )
			return 0 ;
			
		std::vector<struct _overlap> &overlaps = *secondaryGeneOverlaps ;
		secondaryOverlapIdx.resize( 0 ) ;
		int size = overlaps.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			struct _overlap &overlap = overlaps[i] ;
			
			if ( GetGeneType( seqs[ overlap.seqIdx ].name ) != geneType )
				continue ;
			if ( primaryOverlap.similarity == overlap.similarity &&
				primaryOverlap.readEnd - primaryOverlap.readStart  == overlap.readEnd - overlap.readStart )
				// &&primaryOverlap.seqEnd - primaryOverlap.seqStart == overlap.seqEnd - overlap.seqStart )
				secondaryOverlapIdx.push_back( i ) ;
		}
		return secondaryOverlapIdx.size() ;
	}
	
	// Convert the annotation information to string.
	void AnnotationToString( char *read, struct _overlap geneOverlap[4], struct _overlap cdr[3],
	                std::vector<struct _overlap> *secondaryGeneOverlaps, bool outputGeneAlignment, char *buffer ) 
	{
		int i, j ;
		buffer[0] = '\0' ;
		char *r = strdup( read ) ; // another buffer

		// Output the V, J, C information
		for ( i = 0 ; i < 4 ; ++i )
		{
			//if ( i == 1 )
			//	continue ;
			if ( geneOverlap[i].seqIdx != -1 )
			{
				int seqIdx = geneOverlap[i].seqIdx ;
				sprintf( buffer + strlen( buffer ), " %s(%d):(%d-%d):(%d-%d):%.2lf",
						seqs[ seqIdx ].name, seqs[ seqIdx ].consensusLen,
						geneOverlap[i].readStart, geneOverlap[i].readEnd, 
						geneOverlap[i].seqStart, geneOverlap[i].seqEnd, geneOverlap[i].similarity * 100 ) ;

				if ( secondaryGeneOverlaps != NULL )
				{
					std::vector<struct _overlap> &overlaps = *secondaryGeneOverlaps ;
					int size = overlaps.size() ;
					for ( j = 0 ; j < size ; ++j )
					{
						seqIdx = overlaps[j].seqIdx ;
						if ( GetGeneType( seqs[ seqIdx ].name ) != i )
							continue ;
						sprintf( buffer + strlen( buffer ), ",%s(%d):(%d-%d):(%d-%d):%.2lf",
							seqs[ seqIdx ].name, seqs[ seqIdx ].consensusLen,
							overlaps[j].readStart, overlaps[j].readEnd, 
							overlaps[j].seqStart, overlaps[j].seqEnd, overlaps[j].similarity * 100 ) ;
					}
				}
			}
			else
			{
				sprintf( buffer + strlen( buffer ), " *" ) ;
			}
		}

		// Output CDR1,2,3 information.
		//fprintf( stderr, "ERROR %s\n", read ) ;
		for ( i = 0 ; i < 3 ; ++i )
		{
			//printf( "%d %d\n", cdr[i].readStart, cdr[i].readEnd ) ;
			if ( cdr[i].seqIdx != -1 )
			{
				int len = cdr[i].readEnd - cdr[i].readStart + 1 ;
				memcpy( r, read + cdr[i].readStart, len ) ;
				//if ( len >= strlen( read ) || cdr[i].readStart < 0 || cdr[i].readStart + len - 1 >= strlen( read ))
				//	fprintf( stderr, "ERROR %s\n", read ) ;
			        r[len] = '\0' ;
				
				sprintf( buffer + strlen( buffer), " CDR%d(%d-%d):%.2lf=%s", i + 1, cdr[i].readStart, cdr[i].readEnd, 
						cdr[i].similarity * 100, r ) ;
			}
			else
				sprintf( buffer + strlen( buffer), " CDR%d(0-0):0.00=null", i + 1 ) ;
		}

		if ( outputGeneAlignment )
		{
			for ( i = 0 ; i < 4 ; ++i )
			{
				//if ( i == 1 )
				//	continue ;

				char *align = GetGeneOverlapAlignment( read, geneOverlap[i] ) ;
				if ( align != NULL )
				{
					int k, l ;
					j = geneOverlap[i].readStart ;
					k = geneOverlap[j].seqStart ;
					
					// Convert the alignment results to another type of string, but showed 
					//   how to obtain the read sequence from reference sequence.
					// =: match. [ACGT]: mismatch. [\]: deletion. [acgt]: insertion. 
					for ( l = 0 ; align[l] != -1 ; ++l )
					{
						if ( align[l] == EDIT_MATCH )
						{
							align[l] = '=' ;
							++j ; ++k ;
						}
						else if ( align[l] == EDIT_MISMATCH )
						{
							align[l] = read[j] ;
							++j ; ++k ;
						}
						else if ( align[l] == EDIT_DELETE )
						{
							align[l] = '\\' ;
							++k ;
						}
						else if ( align[l] == EDIT_INSERT )
						{
							align[l] = read[j] - 'A' + 'a' ;
							++j ;
						}
					}
					align[l] = '\0' ;

					sprintf( buffer + strlen( buffer ), " %s", align ) ;
					delete[] align ;
				}
				else
				{
					sprintf( buffer + strlen( buffer ), " *" ) ;
				}
			}
		}
		free( r ) ;
	}

	char *GetGeneOverlapAlignment( char *read, struct _overlap gene )
	{
		if ( gene.seqIdx == -1 )
			return NULL ;

		int len = strlen( read ) ;
		char *align ;
		align = new char[ len + seqs[gene.seqIdx].consensusLen + 2 ] ;
		AlignAlgo::GlobalAlignment( seqs[ gene.seqIdx ].consensus + gene.seqStart,
				gene.seqEnd - gene.seqStart + 1,
				read + gene.readStart, gene.readEnd - gene.readStart + 1, align ) ;
		return align ;
	}
	
	// Use the refSet to annotate current set.
	void Annotate( SeqSet &refSet )
	{
		int i ;
		char *buffer = new char[10240] ;
		int seqCnt = seqs.size() ;
		struct _overlap geneOverlap[4];
		struct _overlap cdr[3] ;
		for ( i = 0 ; i < seqCnt  ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;
		
			free( seqs[i].name ) ;
			refSet.AnnotateRead( seqs[i].consensus, 2, geneOverlap, cdr, NULL ) ;
			seqs[i].name = strdup( buffer ) ;
		}

		delete[] buffer ;
	}
	
	
	void BreakFalseAssembly( std::vector<struct _Read> reads )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		std::vector<struct _overlap> *adj = new std::vector<struct _overlap>[ seqCnt ] ;
		//BuildBranchGraph( adj, 31 ) ;
		
		int readCnt = reads.size() ;
		
		// Mark whether to trust the branch overlap portion.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			std::vector<struct _overlap> overlaps ;
			int overlapCnt = GetOverlapsFromRead( reads[i].seq, 1, -1, 1, false, overlaps ) ;
			if ( overlapCnt == 0 )
				continue ;

			std::sort( overlaps.begin(), overlaps.end() ) ;
			
			struct _overlap extendedOverlap ;
			int len = strlen( reads[i].seq ) ;
			char *rc = strdup( reads[i].seq ) ;
			ReverseComplement( rc, reads[i].seq, len ) ;
			char *r = reads[i].seq ;
			if ( overlaps[0].strand == -1 )
				r = rc ;
			char *align = new char[ 2 * len + 2 ] ;	
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( ExtendOverlap( r, len, seqs[ overlaps[j].seqIdx ], 1.0, align, 
					overlaps[j], extendedOverlap ) == 1 )
				{
					int seqIdx = extendedOverlap.seqIdx ;
					int adjCnt = adj[ seqIdx ].size() ;
					for ( k = 0 ; k < adjCnt ; ++k )
					{
						if ( extendedOverlap.seqStart < adj[ seqIdx ][k].readStart &&
							extendedOverlap.seqEnd > adj[ seqIdx ][k].readEnd )
						{
							adj[ seqIdx ][k].similarity = 0 ;
						}
					}
					break ;
				}
			}
			free( rc ) ;
			delete[] align ;
		}

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( adj[i][j].similarity < repeatSimilarity )
					adj[i][j].similarity = 0 ;

				//printf( "%d branch %d %d: %d %d %lf\n", i, j, adj[i][j].seqIdx, 
				//	adj[i][j].readStart, adj[i][j].readEnd, adj[i][j].similarity ) ;
			}
		}

		// Break up the Seq 
		int newSeqId = seqCnt ;
		KmerCode kmerCode( kmerLength ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int adjCnt = adj[i].size() ;
			// Find the regions.
			SimpleVector<struct _pair> breakCoords ;
			breakCoords.Reserve( adjCnt ) ;
			for ( k = 0 ; k < adjCnt ; ++k )
			{
				if ( adj[i][k].similarity <= 0 )
					continue ;
				
				struct _pair nb ;
				nb.a = adj[i][k].readStart ;
				nb.b = adj[i][k].readEnd ;
				if ( nb.a == 0 || nb.b == seqs[i].consensusLen - 1 )
					continue ;

				breakCoords.PushBack( nb ) ;
			}
			std::sort( breakCoords.BeginAddress(), breakCoords.EndAddress(), CompSortPairAInc ) ;	
			int size = breakCoords.Size() ;

			if ( size == 0 )
				continue ;
			k = 1 ;
			for ( j = 1 ; j < size ; ++j )
			{
				if ( breakCoords[j].a <= breakCoords[k - 1].b )
				{
					if ( breakCoords[j].b > breakCoords[k - 1].b )
						breakCoords[k - 1].b = breakCoords[j].b ;
				}
				else
				{
					breakCoords[k] = breakCoords[j] ;
					++k ;
				}
			}
			breakCoords.Resize( k ) ;
			
			int start, end ;
			//for ( j = 0 ; j < k ; ++j )
			//	printf( "%d %s breakCoords %d: %d %d\n", i, seqs[i].name, j, breakCoords[j].a, breakCoords[j].b ) ;
			for ( j = 0 ; j <= k ; ++j )
			{
				if ( j == 0 )
					start = 0 ;
				else
					start = breakCoords[j - 1].a ;

				if ( j == k )
					end = seqs[i].consensusLen - 1 ;
				else
					end = breakCoords[j].b ;

				struct _seqWrapper ns ;
				ns.consensus = (char *)malloc( sizeof( char ) * ( end - start + 2 ) ) ;
				memcpy( ns.consensus, seqs[i].consensus + start, end - start + 1 ) ;
				ns.consensus[ end - start + 1 ] = '\0' ;
				ns.consensusLen = end - start + 1 ;
				ns.name = strdup( seqs[i].name ) ;
				ns.isRef = false ;
				
				ns.posWeight.Reserve( end - start + 1 ) ;
				int l ;
				for ( l = start ; l <= end ; ++l )
					ns.posWeight.PushBack( seqs[i].posWeight[l] ) ;	

				seqs.push_back( ns ) ;
				printf( "Break %d to %d.\n", i, seqs.size() - 1 ) ;
			}
			// Clean up original seq
			seqIndex.RemoveIndexFromRead( kmerCode, seqs[i].consensus, seqs[i].consensusLen, i, 0 ) ;
			free( seqs[i].consensus ) ;
			free( seqs[i].name ) ;
			seqs[i].consensus = seqs[i].name = NULL ;
			seqs[i].posWeight.Release() ;
		}
		//Clean( true ) ;	
	}

	// The readId is for the first mate, the other mate readId is readId+1.
	//   The main program should determine which one to use.
	void UpdateMateAdjGraph( int from, int fromStart, int fromEnd, int to, int toStart, int toEnd, 
			int readId, std::vector<struct _overlap> *mateAdj )
	{
		int j ;
		int size = mateAdj[from].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			if ( mateAdj[from][j].seqIdx == to )
			{
				if ( fromStart < mateAdj[from][j].readStart )
					mateAdj[from][j].readStart = fromStart ;
				if ( fromEnd > mateAdj[from][j].readEnd )
					mateAdj[from][j].readEnd = fromEnd ;
				if ( toStart < mateAdj[from][j].seqStart )
					mateAdj[from][j].seqStart = toStart ;
				if ( toEnd > mateAdj[from][j].seqEnd )
					mateAdj[from][j].seqEnd = toEnd ;
				++mateAdj[from][j].matchCnt ;
				mateAdj[from][j].info->PushBack( readId ) ;
				break ;
			}
		}
		if ( j >= size )
		{
			struct _overlap na ;
			na.seqIdx = to ;
			na.readStart = fromStart ;
			na.readEnd = fromEnd ;
			na.seqStart = toStart ;
			na.seqEnd = toEnd ;
			na.matchCnt = 1 ;
			na.similarity = 0 ;
			na.info = new SimpleVector<int> ;
			na.info->PushBack( readId ) ;
			mateAdj[from].push_back( na ) ;
		}
	}

	// For assemblies with no extension at all, allow gaps, and put that information in branch graph.
	void AugmentBranchGraphByMate( std::vector<struct _overlap> *branchAdj, 
			std::vector<struct _overlap> *prevAdj, std::vector<struct _overlap> *nextAdj, 
			std::vector<struct _assignRead> &reads )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( prevAdj[i].size() + nextAdj[i].size() != 1 )
				continue ;
			struct _seqWrapper &seq = seqs[i] ;
			/*for ( j = 0 ; j < seq.consensusLen ; ++j )
			{
				if ( seq.posWeight[j].Sum() != 1 )
					break ;
			}
			if ( j < seq.consensusLen )
				continue ;*/		
			int bsize = branchAdj[i].size() ;
			int readIdx, mateReadIdx ;
			char *fr ; // first read
			char *sr ; // second read;
			int flen, slen ; // their lengths

			// For overlaps
			int offset = -1 ;
			int overlapSize = -1 ;
			int bestMatchCnt = -1 ;

			if ( prevAdj[i].size() == 1 )
			{
				// Extend left
				for ( j = 0 ; j < bsize ; ++j )
				{
					if ( branchAdj[i][j].seqIdx == prevAdj[i][0].seqIdx )
						break ;
				}

				if ( j < bsize )
					continue ;

				//if ( prevAdj[i][0].info->Size() == 0 )
				//	continue ;
				int size = prevAdj[i][0].info->Size() ;
				for ( j = 0 ; j < size ; ++j )
				{
					readIdx = prevAdj[i][0].info->Get(j) ;
					mateReadIdx = readIdx + 1 ;
					if ( reads[ readIdx ].overlap.seqIdx != i )
					{
						int tmp = readIdx ;
						readIdx = mateReadIdx ;
						mateReadIdx = tmp ;
					}

					flen = strlen( reads[ mateReadIdx ].read  ) ;
					slen = strlen( reads[ readIdx ].read ) ;
					//if ( slen != seq.consensusLen || flen != seqs[ prevAdj[i][0].seqIdx ].consensusLen  )
					//	continue ;
					if ( reads[ readIdx ].overlap.seqStart > 3 )
						continue ;
					
					fr = strdup( reads[ mateReadIdx ].read ) ;
					sr = strdup( reads[ readIdx ].read ) ;
					int minOverlap = ( flen + slen) / 20 ;
					if ( minOverlap >= 31 )
						minOverlap = 31 ;

					if ( reads[ mateReadIdx ].overlap.strand == -1 )
						ReverseComplement( fr, reads[ mateReadIdx ].read, flen ) ;
					if ( reads[ readIdx ].overlap.strand == -1 )
						ReverseComplement( sr, reads[ readIdx ].read, slen ) ;
					overlapSize = AlignAlgo::IsMateOverlap( fr, flen, sr, slen, minOverlap, offset, bestMatchCnt ) ;
				
					if ( overlapSize == -1 )
					{
						free( fr ) ; free( sr ) ;
					}
					else
						break ;
				}
			}
			else
			{
				// Extend right
				for ( j = 0 ; j < bsize ; ++j )
				{
					if ( branchAdj[i][j].seqIdx == nextAdj[i][0].seqIdx )
						break ;
				}
				if ( j < bsize )
					continue ;

				//if ( nextAdj[i][0].info->Size() == 0 )
				//	continue ;
				int size = nextAdj[i][0].info->Size() ;
				for ( j = 0 ; j < size ; ++j )
				{
					readIdx = nextAdj[i][0].info->Get(j) ;
					mateReadIdx = readIdx + 1 ;

					if ( reads[ readIdx ].overlap.seqIdx != i )
					{
						int tmp = readIdx ;
						readIdx = mateReadIdx ;
						mateReadIdx = tmp ;
					}

					flen = strlen( reads[ readIdx ].read ) ;
					slen = strlen( reads[ mateReadIdx ].read ) ;
					//if ( flen != seq.consensusLen /*|| slen != seqs[ nextAdj[i][0].seqIdx ].consensusLen */)
					if ( reads[ readIdx ].overlap.seqEnd < seq.consensusLen - 4 )
						continue ;
					int minOverlap = ( flen + slen) / 20 ;
					if ( minOverlap >= 31 )
						minOverlap = 31 ;

					fr = strdup( reads[ readIdx ].read ) ;
					sr = strdup( reads[ mateReadIdx ].read ) ;

					if ( reads[ readIdx ].overlap.strand == -1 )
						ReverseComplement( fr, reads[ readIdx ].read, flen ) ;
					if ( reads[ mateReadIdx ].overlap.strand == -1 )
						ReverseComplement( sr, reads[ mateReadIdx ].read, slen ) ;
					overlapSize = AlignAlgo::IsMateOverlap( fr, flen, sr, slen, minOverlap, offset, bestMatchCnt ) ;
					
					if ( overlapSize == -1 )
					{
						free( fr ) ; free( sr ) ;
					}
					else
						break ;
				}
			}
			

			if ( overlapSize == -1 ) // Ambiguous overlap or no overlap
			{
				//free( fr ) ; free( sr ) ;
				continue ;
			}
			bestMatchCnt *= 2 ;
			if ( prevAdj[i].size() == 1 )
			{
				struct _overlap no ;
				no.seqIdx = prevAdj[i][0].seqIdx ;
				no.readStart = reads[ readIdx ].overlap.seqStart ;
				no.readEnd = no.readStart + overlapSize - 1 ;
				no.seqStart = reads[ mateReadIdx ].overlap.seqStart + offset ;
				no.seqEnd = no.seqStart + overlapSize - 1 ;
				no.matchCnt = bestMatchCnt ;
				no.similarity = bestMatchCnt / (2.0 * overlapSize ) ;
				branchAdj[i].push_back( no ) ;
				
				//printf( "<=%d %d\n%s\n%s\n\n", no.seqIdx, i, fr, sr ) ;
			}
			else
			{
				struct _overlap no ;
				no.seqIdx = nextAdj[i][0].seqIdx ;
				no.readStart = reads[ readIdx ].overlap.seqStart + offset ;
				no.readEnd = no.readStart + overlapSize - 1 ;
				no.seqStart = reads[ mateReadIdx ].overlap.seqStart ;
				no.seqEnd = no.seqStart + overlapSize - 1 ;
				no.matchCnt = bestMatchCnt ;
				no.similarity = bestMatchCnt / (2.0 * overlapSize) ;
				branchAdj[i].push_back( no ) ;
				//printf( ">=%d %d\n%s\n%s\n\n", i, no.seqIdx, fr, sr ) ;
			}
			free( fr ) ; free( sr ) ;
		}
	}
	
	// return: 0: no extension. 1: end-to-inside extension. 2: end-to-end extension
	int GetExtendSeqCoord( int from, struct _overlap mateInfo, int direction, 
			std::vector<struct _overlap> *branchAdj, bool aggressive, struct _overlap &coord )
	{
		int i, k ;
		int adjCnt = branchAdj[from].size() ;	
		int to = mateInfo.seqIdx ;

		coord.seqIdx = -1 ;
		int overhangSize = 5 ; // Allowing the end to have this much random extension in the raw assembly.
		
		for ( i = 0 ; i < adjCnt ; ++i )
		{
			if ( direction == 1 )
			{
				if ( branchAdj[from][i].seqIdx == to 
						&& branchAdj[from][i].readEnd >= seqs[from].consensusLen - overhangSize ) 
						//&& branchAdj[from][i].seqStart <= 5 - 1 )
					break ;
			}
			else if ( direction == -1 )
			{
				if ( branchAdj[from][i].seqIdx == to 
						&& branchAdj[from][i].readStart <= overhangSize - 1 ) 
						//&& branchAdj[from][i].seqEnd <= seqs[to].consensusLen - 5 )
					break ;
			}
		}


		if ( i >= adjCnt )
			return 0 ;
		k = i ;
		
		if ( direction == 1 && mateInfo.seqEnd <= branchAdj[from][k].seqEnd )
			return 0 ;
		else if ( direction == -1 && mateInfo.seqStart >= branchAdj[from][k].seqStart )
			return 0 ;

		/*for ( i = 0 ; i < adjCnt ; ++i )
		{
		}*/
		
		coord.seqIdx = to ;
		coord.matchCnt = branchAdj[from][k].readEnd - branchAdj[from][k].readStart + 1 ; // Record the length of the branch overlap.
		int ret = 1 ;
		
		if ( direction == 1 )
		{
			coord.readStart = 0 ;
			coord.readEnd = branchAdj[from][k].readEnd ;
		}
		else
		{
			coord.readStart = branchAdj[from][k].readStart ;
			coord.readEnd = seqs[from].consensusLen - 1 ;
		}

		if ( direction == 1 )
		{
			coord.seqStart = branchAdj[from][k].seqEnd + 1 ;
			if ( aggressive )
				coord.seqEnd = seqs[to].consensusLen - 1 ;
			else
				coord.seqEnd = mateInfo.seqEnd ;

			if ( branchAdj[from][k].seqStart <= overhangSize - 1 )
				ret = 2 ;
		}
		else
		{
			if ( aggressive )
				coord.seqStart = 0 ;
			else
			{
				//printf( "%d %d: %d %d\n", from, to, mateInfo.seqStart, branchAdj[from][k].seqStart ) ;
				coord.seqStart = mateInfo.seqStart ;
			}
			coord.seqEnd = branchAdj[from][k].seqStart - 1 ;

			if ( branchAdj[from][k].seqEnd >= seqs[ branchAdj[from][k].seqIdx ].consensusLen - overhangSize )
				ret = 2 ;
		}

		return ret ;
	}
	
	bool CanGapExtend( int from, struct _overlap mateInfo, int direction, std::vector<struct _overlap> *branchAdj )	
	{
		int i, j ;
		int size = branchAdj[from].size() ;
		if ( branchAdj != NULL )
		{
			for ( i = 0 ; i < size ; ++i )
				if ( branchAdj[from][i].seqIdx == mateInfo.seqIdx )
				{
					int bs, be, ms, me ;
					bs = branchAdj[from][i].seqStart ;
					be = branchAdj[from][i].seqEnd ;
					ms = mateInfo.seqStart ;
					me = mateInfo.seqEnd ;

					if ( bs <= ms && be >= me )
						return false ;
					if ( ms <= bs && me >= be ) 
						return false ;
					if ( be <= ms || bs >= me )
						continue ;

					if ( bs <= ms && be <= me && be - ms + 1 >= 17 )
						return false ;
					else if ( bs >= ms && be >= me && me - bs + 1 >= 17 )
						return false ;
				}
		}

		if ( direction == -1 )
		{
			if ( mateInfo.readStart < 50 )
				return true ;

			for ( i = 0 ; i < 3 ; ++i )	
				if ( seqs[from].info[i].a != -1 )
					break ;
			if ( i < 3 )
			{
				int to = mateInfo.seqIdx ;
				for ( j = 0 ; j < 3 ; ++j )
					if ( seqs[to].info[j].b + 3 >= mateInfo.seqEnd )
						break ;

				if ( j < i )
					return true ;
			}
		}
		else
		{
			if ( mateInfo.readEnd >= seqs[ from ].consensusLen - 50 )
				return true ;

			for ( i = 2 ; i >= 0 ; --i )
				if ( seqs[from].info[i].a != -1 )
					break ;
			if ( i >= 0 )
			{
				int to = mateInfo.seqIdx ;
				for ( j = 2 ; j >= 0 ; --j )
					if ( seqs[to].info[j].a >= 0 && seqs[to].info[j].a - 3 <= mateInfo.seqStart )
						break ;

				if ( i < j )
					return true ;
			}
		}

		return false ;
	}

	// Test whether it is reasonable to extend 
	int GetGapExtendSeqCoord( int from, struct _overlap mateInfo, int direction, struct _overlap &coord )
	{
		if ( direction == -1 )
		{
			coord = mateInfo ; 
			coord.readStart = 0 ;
			coord.readEnd = seqs[from].consensusLen - 1 ; 
			coord.matchCnt = 0 ;
			return 1 ;
		}
		else 
		{
			coord = mateInfo ; 
			coord.readStart = 0 ;
			coord.readEnd = seqs[from].consensusLen - 1 ; 
			coord.matchCnt = 0 ;
			return 1 ;
		}
		return 0 ;
	}
	
	// Extend seq where one mate is aligned the other mate is missing, and the mates overlap with each other.
	void ExtendSeqFromMissingOverlapMate( std::vector< struct _assignRead> reads )
	{
		int i, j, k ;
		int readCnt = reads.size() ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < readCnt ; ++i )
		{
			bool paired = false ;
			if ( i < readCnt - 1 && !strcmp( reads[i].id, reads[i + 1].id ) )
				paired = true ;
			if ( !paired )
				continue ;

			if ( ( reads[i].overlap.seqIdx != -1 && reads[i + 1].overlap.seqIdx != -1 )
				|| ( reads[i].overlap.seqIdx == -1 && reads[i + 1].overlap.seqIdx == -1 ) )
			{
				++i ;
				continue ;
			}
			
			int anchorId = i ;
			int anchorSeqIdx = reads[i].overlap.seqIdx ;
			int hangId = i + 1 ;
			int asLen, hsLen ;
			if ( reads[i].overlap.seqIdx == -1 )		
			{
				anchorId = i + 1 ;
				anchorSeqIdx = reads[i + 1].overlap.seqIdx ;
				hangId = i ;
			}
		
			struct _seqWrapper &seq = seqs[ anchorSeqIdx ] ;
			//printf( "%d %s\n%s\n", anchorSeqIdx, seq.consensus, reads[ anchorId ].read ) ;
			if ( ( reads[ anchorId ].overlap.strand == 1 && strcmp( seq.consensus, reads[ anchorId ].read ) ) 
				|| ( reads[ anchorId ].overlap.strand == -1 && !IsReverseComplement( seq.consensus, reads[ anchorId ].read ) ) )
			{
				++i ;
				continue ;
			}
			asLen = strlen( reads[ anchorId ].read ) ;
			hsLen = strlen( reads[ hangId ].read ) ;
			
			char *r = strdup( reads[ hangId ].read ) ;
			char *fr, *sr ;
			int direction = 0 ;
			int flen, slen ;
			if ( reads[ anchorId ].overlap.strand == 1 )
			{
				// Extend towards right
				ReverseComplement( r, reads[ hangId ].read, hsLen ) ;

				fr = reads[ anchorId ].read ;
				flen = asLen ;
				sr = r ;
				slen = hsLen ;
				direction = 1 ;
			}
			else
			{
				fr = r ;
				flen = hsLen ;
				sr = seq.consensus ;
				slen = asLen ;
				direction = -1 ;
			}
			
			int offset ;
			int matchCnt ;
			int minOverlap = ( flen + slen) / 20 ;
			if ( minOverlap >= 31 )
				minOverlap = 31 ;
			//printf( "overlap test %d %d %d\n%s\n%s\n", flen, slen, minOverlap, fr, sr ) ;
			if ( AlignAlgo::IsMateOverlap( fr, flen, sr, slen, minOverlap, 
				offset, matchCnt ) == -1 )
			{
				++i ;
				free( r ) ;
				continue ;
			}

			if ( flen - offset >= slen ) // sr is contained in fr
			{
				++i ;
				free( r ) ;
				continue ;
			}

			int newConsensusLen = offset + slen ; 
			char *newConsensus = (char*)malloc( sizeof( char ) * ( newConsensusLen + 1 ) ) ;
			memcpy( newConsensus, fr, offset ) ;
			memcpy( newConsensus + offset, sr, slen ) ;
			newConsensus[newConsensusLen] = '\0' ;
			
			// Update the pos weight.
			if ( direction == 1 )
			{
				// Append
				seq.posWeight.ExpandTo( newConsensusLen ) ;
				seq.posWeight.SetZero( flen, newConsensusLen - flen ) ;
				UpdatePosWeightFromRead( seq.posWeight, offset, sr ) ;
				/*for ( j = 0 ; j < slen ; ++j )
				{
					if ( sr[j] == 'N' )
						continue ;
					++seq.posWeight[j + offset].count[ nucToNum[ sr[j] - 'A' ] ] ;
				}*/
			}
			else
			{
				seq.posWeight.ShiftRight( offset ) ;
				seq.posWeight.SetZero( 0, offset ) ;
				UpdatePosWeightFromRead( seq.posWeight, 0, fr ) ;
				/*for ( j = 0 ; j < flen ; ++j )
				{
					if ( fr[j] == 'N' )
						continue ;
					++seq.posWeight[j].count[ nucToNum[ fr[j] - 'A' ] ] ;
				}*/
			}
			
			free( r ) ;
			//printf( "pass\n%s\n%s\n%d %d\n", seq.consensus, newConsensus, seq.consensusLen, newConsensusLen ) ;
			//fflush( stdout ) ;
			free( seq.consensus ) ;
			seq.consensus = newConsensus ;
			seq.consensusLen = newConsensusLen ;
			seq.minLeftExtAnchor = seq.minRightExtAnchor = 0 ;
			++i ;
		}
	}

	// Use this set of reads to extend,rearrange the seq 
	void ExtendSeqFromReads( std::vector<struct _assignRead> &reads, int leastOverlapLen, SeqSet &refSet )
	{
		int i, j, k ;
		int readCnt = 0 ;
		int seqCnt = seqs.size() ;
		double backupNovelSeqSimilarity = novelSeqSimilarity ;

		novelSeqSimilarity = 1.00 ;
		int ret = 0 ;

		std::vector<struct _overlap> *branchAdj = new std::vector<struct _overlap>[ seqCnt ] ;

		// Mate-pair support graph. In here, start,end represnts the rough range with read support
		//	and matchCnt represent how many mates support this connection.
		std::vector<struct _overlap> *nextAdj = new std::vector<struct _overlap>[ seqCnt ] ; 
		std::vector<struct _overlap> *prevAdj = new std::vector<struct _overlap>[ seqCnt ] ; 
		readCnt = reads.size() ;
		
		SimpleVector<bool> useInBranch ;
		useInBranch.ExpandTo( seqCnt ) ;
		useInBranch.SetZero( 0, seqCnt ) ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			seqs[i].info[0].a = seqs[i].info[0].b = i ;
		}

		// Directly use overlapped mate pairs for extension on 
		//   singleton assembly where the other mate has no perfect alignment.
		// Since these reads will not be applied on more sophisticated extension, it is fine 
		//   to make it an independent component.
		// We are merging reads at beginning, so there is no need for this procedure now.
		std::sort( reads.begin(), reads.end(), CompSortAssignedReadById ) ;
		//ExtendSeqFromMissingOverlapMate( reads ) ;
		
		// Then do more sophisticated extension
		// Build the mate adj graph.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			bool paired = false ;
			if ( i < readCnt - 1 && !strcmp( reads[i].id, reads[i + 1].id ) )
				paired = true ; 
			
			if ( paired )
			{
				if ( reads[i].overlap.seqIdx == -1 || reads[i + 1].overlap.seqIdx == -1 
					|| reads[i].overlap.strand == reads[i + 1].overlap.strand 
					|| ( reads[i].overlap.similarity < 1 && reads[i + 1].overlap.similarity < 1 ) 
					//|| reads[i].overlap.similarity < 0.95 || reads[i + 1].overlap.similarity < 0.95
					|| reads[i].overlap.seqIdx == reads[i + 1].overlap.seqIdx )			
				{
					++i ;
					continue ;
				}
				int from, to ;
				int fromStart, fromEnd, toStart, toEnd ;
				struct _overlap extendedOverlap = reads[i].overlap ;
				struct _overlap mateExtendedOverlap = reads[i + 1].overlap ;
				
				bool validNext = true ;
				bool validPrev = true ;
				if ( extendedOverlap.strand == 1 )
				{
					from = extendedOverlap.seqIdx ;
					fromStart = extendedOverlap.seqStart ;
					fromEnd = extendedOverlap.seqEnd ;

					to = mateExtendedOverlap.seqIdx ;
					toStart = mateExtendedOverlap.seqStart ;
					toEnd = mateExtendedOverlap.seqEnd ;

					if ( reads[i].overlap.similarity < 1 )
						validNext = false ;
					if ( reads[i + 1].overlap.similarity < 1 )
						validPrev = false ;
				}
				else
				{
					to = extendedOverlap.seqIdx ;
					toStart = extendedOverlap.seqStart ;
					toEnd = extendedOverlap.seqEnd ;

					from = mateExtendedOverlap.seqIdx ;
					fromStart = mateExtendedOverlap.seqStart ;
					fromEnd = mateExtendedOverlap.seqEnd ;

					if ( reads[i + 1].overlap.similarity < 1 )
						validNext = false ;
					if ( reads[i].overlap.similarity < 1 )
						validPrev = false ;
				}
			
				useInBranch[from] = true ;
				useInBranch[to] = true ;
				if ( validNext )
					UpdateMateAdjGraph( from, fromStart, fromEnd, to, toStart, toEnd, i, nextAdj ) ;
				if ( validPrev )
					UpdateMateAdjGraph( to, toStart, toEnd, from, fromStart, fromEnd, i, prevAdj ) ;
				
				//printf( "%d(%d %d) %d(%d %d)\n", 
				//	from, fromStart, fromEnd,
				//	to, toStart, toEnd ) ;
			}
			else
			{
				// Single-end reads, not sure how to use them for now.
				// 	may be look for fragmented secondary alignment?
			}
			
			if ( paired )
				++i ;
		
		}

		// Sort the edge by the weight.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			std::sort( prevAdj[i].begin(), prevAdj[i].end() ) ;
			std::sort( nextAdj[i].begin(), nextAdj[i].end() ) ;
		}

		// Briefly annotate the reads
		char *buffer = new char[10000] ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( !useInBranch[i] )
			{
				for ( k = 0 ; k < 3 ; ++k )
				{
					seqs[i].info[k].a = -1 ;
					seqs[i].info[k].b = -1 ;
					seqs[i].info[k].c = -1 ;
				}
				continue ;
			}
			struct _overlap geneOverlap[4] ;
			struct _overlap cdr[3] ;
			refSet.AnnotateRead( seqs[i].consensus, 0, geneOverlap, cdr, NULL ) ;
			for ( j = 0 ; j < 4 ; ++j )
			{
				k = j ;
				if ( j > 1 )
					--k ;
				if ( j == 1 )
					continue ;
				else if ( geneOverlap[j].seqIdx == -1 )	
				{
					seqs[i].info[k].a = -1 ;
					seqs[i].info[k].b = -1 ;
					seqs[i].info[k].c = -1 ;
					continue ;
				}
				seqs[i].info[k].a = geneOverlap[j].readStart ;
				seqs[i].info[k].b = geneOverlap[j].readEnd ;
				seqs[i].info[k].c = geneOverlap[j].seqIdx ;
			}

		}
		delete[] buffer ;
		
		int backupHLR = hitLenRequired ;
		hitLenRequired = leastOverlapLen ;
		BuildBranchGraph( branchAdj, leastOverlapLen, prevAdj, nextAdj ) ;
		hitLenRequired = backupHLR ;

		
		// Each Seq will be extend to left and right one step.
		SimpleVector<int> toRemoveSeqIdx ;
		toRemoveSeqIdx.Reserve( seqCnt ) ;
		
		// Find which seq id extend to.
		SimpleVector<struct _pair> matePrevNext ; 
		SimpleVector<struct _pair> matePrevNextType ; // 0: no extension(same with mateprevnext=-1); 
							      // 1: overlapped extension; 2: gapped extension 
		matePrevNext.ExpandTo( seqCnt ) ;
		matePrevNextType.ExpandTo( seqCnt ) ;
#ifdef DEBUG
		for ( i = 0 ; i < seqCnt ; ++i )
			printf( "%d (V:%d-%d) (J:%d-%d) (C:%d-%d)\n%s\n", i, 
				seqs[i].info[0].a, seqs[i].info[0].b,
				seqs[i].info[1].a, seqs[i].info[1].b,
				seqs[i].info[2].a, seqs[i].info[2].b,
				seqs[i].consensus ) ;
#endif

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevAdjCnt = prevAdj[i].size() ;
			int prevTag = -1 ;
			int max = -1 ;
			for ( j = 0 ; j < prevAdjCnt ; ++j )
			{
				if ( prevAdj[i][j].seqIdx == i )
					continue ;

				if ( prevAdj[i][j].matchCnt > max )
				{
					prevTag = j ;
					max = prevAdj[i][j].matchCnt ;
				}
				else if ( prevAdj[i][j].matchCnt >= max * 0.9 )
				{
					if ( prevAdj[i][j].seqEnd - prevAdj[i][j].seqStart > 
							prevAdj[i][ prevTag ].seqEnd - prevAdj[i][prevTag].seqStart )
						prevTag = j ;
				}
#ifdef DEBUG				
				printf( "<= %d: %d %d. %d %d %d %d\n", i, prevAdj[i][j].seqIdx, prevAdj[i][j].matchCnt, 
						prevAdj[i][j].readStart, prevAdj[i][j].readEnd,
						prevAdj[i][j].seqStart, prevAdj[i][j].seqEnd ) ;
#endif
			}

			int nextAdjCnt = nextAdj[i].size() ;
			int nextTag = -1 ;
			max = -1 ;
			for ( j = 0 ; j < nextAdjCnt ; ++j )
			{
				if ( nextAdj[i][j].seqIdx == i )
					continue ;

				if ( nextAdj[i][j].matchCnt > max )
				{
					nextTag = j ;
					max = nextAdj[i][j].matchCnt ;
				}
				else if ( nextAdj[i][j].matchCnt >= max * 0.9 )
				{
					if ( nextAdj[i][j].seqEnd - nextAdj[i][j].seqStart > 
							nextAdj[i][ nextTag ].seqEnd - nextAdj[i][nextTag].seqStart )
						nextTag = j ;
				}

#ifdef DEBUG
				printf( "=> %d: %d %d. %d %d %d %d\n", i, nextAdj[i][j].seqIdx, nextAdj[i][j].matchCnt, 
						nextAdj[i][j].readStart, nextAdj[i][j].readEnd,
						nextAdj[i][j].seqStart, nextAdj[i][j].seqEnd ) ;
#endif
			}
			matePrevNext[i].a = prevTag ;
			matePrevNext[i].b = nextTag ;
			matePrevNextType[i].a = matePrevNextType[i].b = 0 ;
		}
		//AugmentBranchGraphByMate( branchAdj, matePrevNext, prevAdj, nextAdj, reads ) ;

		// Figure out the seqs whose extension will create repeat sequence
		// 	e.g.: end-to-end extension from two anchor seq.
		SimpleVector<struct _pair> extensionType ;
		SimpleVector<int> uniqueSuccessorOf ;
		extensionType.ExpandTo( seqCnt ) ;
		uniqueSuccessorOf.ExpandTo( seqCnt ) ;
		
		// Filter out the connection that the heaviest mate pair edge has no sequence overlap with 
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;
			
			struct _overlap leftExtend, rightExtend ;
			leftExtend.seqIdx = -1 ;
			rightExtend.seqIdx = -1 ;
			extensionType[i].a = extensionType[i].b = 0 ;
			matePrevNextType[i].a = matePrevNextType[i].b = -1 ;
			if ( prevTag >= 0 )
			{
				extensionType[i].a = GetExtendSeqCoord( i, prevAdj[i][ prevTag ], -1, branchAdj, false, leftExtend ) ;
				if ( leftExtend.seqIdx != -1 )
					matePrevNextType[i].a = 1 ;
				else if ( CanGapExtend( i, prevAdj[i][ prevTag ], -1, branchAdj ) )
				{
					//printf( "hi %d<N=%d\n", prevAdj[i][prevTag], i ) ;
					//|| GetGapExtendSeqCoord( i, prevAdj[i][ prevTag ], -1, leftExtend ) )
					matePrevNextType[i].a = 2 ;
				}
				else
				{
					// Try to use the secondary mate edge
					int size = prevAdj[i].size() ;
					int threshold = prevAdj[i][0].matchCnt * 0.5 ; // prevTag>=0 makes sure size>=0	
					int found = 0 ;
					for ( j = 0 ; j < size ; ++j )
					{
						if ( prevAdj[i][j].matchCnt < threshold )
							break ;

						extensionType[i].a = GetExtendSeqCoord( i, prevAdj[i][j], -1, branchAdj, false, leftExtend ) ;
						if ( leftExtend.seqIdx != -1 )
						{
							matePrevNextType[i].a = 1 ;
							matePrevNext[i].a = j ;
							found = 1 ;
							break ;
						}
						else if ( CanGapExtend( i, prevAdj[i][j], -1, branchAdj ) )
						{
							matePrevNextType[i].a = 2 ;
							matePrevNext[i].a = j ;
							found = 1 ;
							break ;
						}
					}
					if ( !found )
						matePrevNext[i].a = -1 ;
				}
					
			}
			if ( nextTag >= 0 )
			{
				extensionType[i].b = GetExtendSeqCoord( i, nextAdj[i][ nextTag ], 1, branchAdj, false, rightExtend ) ;
				if ( rightExtend.seqIdx != -1 )
					matePrevNextType[i].b = 1 ;
				else if ( CanGapExtend( i, nextAdj[i][nextTag], 1, branchAdj )  )
				{
					// || GetGapExtendSeqCoord( i, nextAdj[i][ nextTag ], 1, rightExtend ) )
					//printf( "hi %d=N>%d\n", i, nextAdj[i][nextTag] ) ;
					matePrevNextType[i].b = 2 ;
				}
				else
				{
					// Try to use the secondary mate edge
					int size = nextAdj[i].size() ;
					int threshold = nextAdj[i][0].matchCnt * 0.5 ;	
					int found = 0 ;
					for ( j = 0 ; j < size ; ++j )
					{
						if ( nextAdj[i][j].matchCnt < threshold )
							break ;

						extensionType[i].b = GetExtendSeqCoord( i, nextAdj[i][j], 1, branchAdj, false, rightExtend ) ;
						if ( leftExtend.seqIdx != -1 )
						{
							matePrevNextType[i].b = 1 ;
							matePrevNext[i].b = j ;
							found = 1 ;
							break ;
						}
						else if ( CanGapExtend( i, nextAdj[i][j], 1, branchAdj ) )
						{
							matePrevNextType[i].b = 2 ;
							matePrevNext[i].b = j ;
							found = 1 ;
							break ;
						}
					}
					if ( !found )
						matePrevNext[i].b = -1 ;
				}
			}
		}

		// Rescue some partial extension. e.g.: 
		// A |;  A<=B. A could not extend right due to no overlap for the heaviest mate, then we will put A=>B back.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;

			if ( prevTag >= 0 )
			{
				int seqIdx = prevAdj[i][ prevTag ].seqIdx ;
				if ( ( matePrevNext[ seqIdx ].b == -1 || matePrevNextType[ seqIdx ].b == 2 )
					&& extensionType[i].a == 2 )
				{
					// rescue the connection.
					int size = nextAdj[ seqIdx ].size() ;
					for ( j = 0 ; j < size ; ++j )
						if ( nextAdj[ seqIdx ][j].seqIdx == i )
						{
							struct _overlap tmp ;
							extensionType[ seqIdx ].b = GetExtendSeqCoord( seqIdx, nextAdj[ seqIdx ][j], 
										1, branchAdj, false, tmp ) ;
							if ( extensionType[ seqIdx ].b == 2 )
							{
								matePrevNext[ seqIdx ].b = j ;
								matePrevNextType[ seqIdx ].b = 1 ;
							}
							break ;
						}
				}
			}

			if ( nextTag >= 0 )
			{
				int seqIdx = nextAdj[i][ nextTag ].seqIdx ;
				if ( ( matePrevNext[ seqIdx ].a == -1 || matePrevNextType[ seqIdx ].a == 2 )
					&& extensionType[i].b == 2 )
				{
					int size = prevAdj[seqIdx].size() ;
					for ( j = 0 ; j < size; ++j )
						if ( prevAdj[ seqIdx ][j].seqIdx == i )
						{
							struct _overlap tmp ;
							extensionType[ seqIdx ].a = GetExtendSeqCoord( seqIdx, prevAdj[ seqIdx ][j], 
										-1, branchAdj, false, tmp ) ;
							if ( extensionType[ seqIdx ].a == 2 )
							{
								matePrevNext[ seqIdx ].a = j ;
								matePrevNextType[ seqIdx ].a = 1 ;
							}
							break ;
						}
				}
			}
		}

		// For gapped extension, they might have short overlap.
		SimpleVector<struct _pair> shortOverlapSeqIdx ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int minOverlap = 10 ;
			if ( matePrevNextType[i].a == 2 )
			{
				int offset = -1 ;
				int bestMatchCnt = -1 ;
				int prevSeqIdx = prevAdj[i][ matePrevNext[i].a ].seqIdx ;
				int overlapSize = AlignAlgo::IsMateOverlap( seqs[prevSeqIdx].consensus, seqs[ prevSeqIdx ].consensusLen, 
							seqs[i].consensus, seqs[i].consensusLen, 
							minOverlap, offset, bestMatchCnt, true ) ;
				if ( overlapSize >= 0 )
				{
					matePrevNextType[i].a = 1 ;
					extensionType[i].a = 2 ;

					// Add that information to the branch graph.
					struct _overlap o ;
					o.readStart = 0 ;
					o.readEnd = overlapSize - 1 ; 
					o.seqIdx = prevSeqIdx ;
					o.seqStart = offset ;
					o.seqEnd = seqs[ prevSeqIdx ].consensusLen - 1 ;
					o.matchCnt = 2 * bestMatchCnt ;
					o.similarity = (double)bestMatchCnt / overlapSize ;  
					branchAdj[i].push_back( o ) ;
					
					struct _pair np ;
					np.a = i ;
					np.b = -1 ;
					shortOverlapSeqIdx.PushBack( np ) ;
				}

			}
			if ( matePrevNextType[i].b == 2 )	
			{
				int offset = -1 ;
				int bestMatchCnt = -1 ;
				int nextSeqIdx = nextAdj[i][ matePrevNext[ i ].b ].seqIdx ;
				int overlapSize = AlignAlgo::IsMateOverlap( seqs[i].consensus, seqs[i].consensusLen, 
						seqs[ nextSeqIdx ].consensus, seqs[nextSeqIdx ].consensusLen,               
						minOverlap, offset, bestMatchCnt, true );
				if ( overlapSize >= 0 ) 
				{
					matePrevNextType[i].b = 1 ;
					extensionType[i].b = 2 ;
					
					// Add that information to the branch graph.
					struct _overlap o ;
					o.readStart = offset ;
					o.readEnd = seqs[i].consensusLen - 1 ;
					o.seqIdx = nextSeqIdx ;
					o.seqStart = 0 ;
					o.seqEnd = overlapSize - 1  ;
					o.matchCnt = 2 * bestMatchCnt ;
					o.similarity = (double)bestMatchCnt /  overlapSize ;
					branchAdj[i].push_back( o ) ;
					
					struct _pair np ;
					np.a = i ;
					np.b = 1 ;
					shortOverlapSeqIdx.PushBack( np ) ;
				}
			}
		}
		// Recompute extensionType for those small overlap seqs
		k = shortOverlapSeqIdx.Size() ;
		for ( i = 0 ; i < k ; ++i )
		{
			int seqIdx = shortOverlapSeqIdx[i].a ;
			struct _overlap extend ;
			if ( shortOverlapSeqIdx[i].b == -1 )
			{
				int prevTag = matePrevNext[ seqIdx ].a ;
				extensionType[ seqIdx ].a = GetExtendSeqCoord( seqIdx, prevAdj[ seqIdx ][ prevTag ], -1, 
					branchAdj, false, extend ) ;
			}
			else
			{
				int nextTag = matePrevNext[ seqIdx ].b ; 
				extensionType[ seqIdx ].b = GetExtendSeqCoord( seqIdx, nextAdj[ seqIdx ][ nextTag ], 1, 
					branchAdj, false, extend ) ;
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			uniqueSuccessorOf[i] = -1 ;
			
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;
			//if ( i == 47 )
			//	fprintf( stderr, "hi %d %d %d %d\n", prevTag, nextTag,
			//		prevAdj[i][prevTag].seqIdx, nextAdj[i][nextTag].seqIdx ) ;
			//if ( prevTag >= 0 && nextTag >= 0 )
			//	continue ;
			/*if ( nextTag >= 0 )
			{
				// Keep the second part.
				if ( extensionType[i].b == 2 )
				{
					int seqIdx = nextAdj[i][nextTag].seqIdx ;
					if ( matePrevNext[ seqIdx ].a >= 0 
						&& prevAdj[ seqIdx ][ matePrevNext[ seqIdx ].a ].seqIdx == i 
						&& extensionType[ seqIdx ].a == 2 )
					{
						directlyFilter[i] = true ;
					}
				}
			}*/

			if ( prevTag >= 0 )
			{
				if ( extensionType[i].a == 2 || matePrevNextType[i].a == 2 )
				{
					int seqIdx = prevAdj[i][ prevTag ].seqIdx ;
					if ( matePrevNext[seqIdx].b >= 0 
						&& nextAdj[ seqIdx ][ matePrevNext[ seqIdx ].b ].seqIdx == i 
						&& ( extensionType[ seqIdx ].b == 2 
							|| matePrevNextType[seqIdx].b == 2 ) )
					{
						uniqueSuccessorOf[i] = seqIdx ;
					}
				}
			}
		}
		

		// Filter some middle part, if its two extension are connected.
		/*for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ; 
			if ( directlyFilter[i] || prevTag < 0 || nextTag < 0 || extensionType[i].a != 2 
				|| extensionType[i].b != 2 )
				continue ;
			
			int prevSeqIdx = prevAdj[i][ prevTag ].seqIdx ;
			int nextSeqIdx = nextAdj[i][ nextTag ].seqIdx ;
			if ( ( matePrevNext[ prevSeqIdx ].b >= 0 
					&& nextAdj[ prevSeqIdx ][ matePrevNext[ prevSeqIdx ].b ].seqIdx == nextSeqIdx ) 
				|| ( matePrevNext[ nextSeqIdx ].a >= 0 
					&& prevAdj[ nextSeqIdx ][ matePrevNext[ nextSeqIdx ].a ].seqIdx == prevSeqIdx ) )
				directlyFilter[i] = true ;
		}*/

		// Do the extesion
		SimpleVector<int> chain ;
		SimpleVector<int> offset ;
		SimpleVector<struct _pair> range ; // the range of a seq that needs to be copied in the new seq.
		SimpleVector<int> origRangeB ;
		SimpleVector<struct _pair> shiftSeq ; // where the seq moved to after extesion. a-new seqIdx, b-offset in new seq.
					              // Note that the ending of a seq might be trimmed when extension,
						      //  the b(offset) will regard those trimmed part existing, and hence
						      //  starts a bit earlier than the portion that is actually used.
		SimpleVector<struct _pair> gapPos ;   // a: start of the gap, b: length of the gap.
		chain.Reserve( seqCnt ) ;	    
		offset.Reserve( seqCnt ) ;
		range.Reserve( seqCnt ) ;
		origRangeB.Reserve( seqCnt ) ; // the range.b that are not adjusted.
		gapPos.Reserve( 10 ) ;

		shiftSeq.ExpandTo( seqCnt ) ; 
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			shiftSeq[i].a = i ;
			shiftSeq[i].b = 0 ;
		}
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			//if ( inCnt[i] > 0 )
			//	continue ;
			// Each seq will try to extend to each direction once.
			if ( uniqueSuccessorOf[i] != -1 )
			{
				toRemoveSeqIdx.PushBack( i ) ;
				continue ;
			}
			
			int first = i ;
			int firstPrevTag = matePrevNext[i].a ;
			int last = i ;
			int lastNextTag = matePrevNext[i].b ;
			
			chain.Clear() ;
			chain.PushBack( i ) ;
			while ( 1 )
			{
				if ( lastNextTag >= 0 && uniqueSuccessorOf[ nextAdj[last][ lastNextTag ].seqIdx ] == last )
				{
					last = nextAdj[last][ lastNextTag ].seqIdx ;
					lastNextTag = matePrevNext[last].b ;
					chain.PushBack( last ) ;
				}
				else
					break ;
			}
			int chainSize = chain.Size() ;
			// Compute the length of new seq.
			int newConsensusLen = 0 ;
			offset.Clear() ;
			range.Clear() ; 
			origRangeB.Clear() ;
			gapPos.Clear() ;
			struct _overlap leftMostExtend, rightMostExtend ;
			leftMostExtend.seqIdx = -1 ;
			rightMostExtend.seqIdx = -1 ;

			offset.ExpandTo( chainSize ) ;
			range.ExpandTo( chainSize ) ;
			origRangeB.ExpandTo( chainSize ) ;
			for ( j = 0 ; j < chainSize ; ++j )
			{
				struct _overlap leftExtend, rightExtend ;
				int prevTag = matePrevNext[ chain[j] ].a ;
				int nextTag = matePrevNext[ chain[j] ].b ;
				leftExtend.seqIdx = -1 ;
				rightExtend.seqIdx = -1 ;
				
				if ( prevTag >= 0 && matePrevNextType[ chain[j] ].a == 1 )
				{
					bool aggressive = true ;
					if ( j == 0 )
					{
						aggressive = false ;
					}
					GetExtendSeqCoord( chain[j], prevAdj[ chain[j] ][ prevTag ], -1, branchAdj, aggressive, leftExtend ) ;
				}
				if ( nextTag >= 0 && matePrevNextType[ chain[j] ].b == 1 )
				{
					bool aggressive = true ;
					if ( j == chainSize - 1 )
					{
						aggressive = false ;
						if ( seqs[ chain[j] ].info[2].c == -1 
							&& seqs[ nextAdj[ chain[j] ][nextTag].seqIdx ].info[2].c != -1 )
						{
							if ( nextAdj[ chain[j] ][nextTag].seqEnd < 
									seqs[ nextAdj[ chain[j] ][nextTag].seqIdx ].info[2].a )
							{
								int size = nextAdj[ chain[j] ].size() ;
								for ( k = 0 ; k < size ; ++k )
								{
									if ( k == nextTag || 
										seqs[ nextAdj[ chain[j] ][k].seqIdx ].info[2].c == -1 )
										continue ;

									if ( seqs[ nextAdj[ chain[j] ][k].seqIdx ].info[2].c 
											== seqs[ nextAdj[ chain[j] ][ nextTag ].seqIdx ].info[2].c
											&& nextAdj[ chain[j] ][k].seqEnd > 
											seqs[ nextAdj[ chain[j] ][k].seqIdx ].info[2].a )
									{
										aggressive = true ;
										break ;
									}
								}
							}
							else
								aggressive = true ;
						}
					}
					GetExtendSeqCoord( chain[j], nextAdj[ chain[j] ][ nextTag ], 1, branchAdj, aggressive, rightExtend ) ;
				}
				if ( matePrevNextType[ chain[j] ].a == 2 )
				{
					GetGapExtendSeqCoord( chain[j], prevAdj[ chain[j] ][ prevTag ], -1 ,leftExtend ) ;
				}
				if ( matePrevNextType[ chain[j] ].b == 2 )
				{
					GetGapExtendSeqCoord( chain[j], nextAdj[ chain[j] ][ nextTag ], 1, rightExtend ) ;
				}
				
				if ( j == 0 )
				{
					if ( leftExtend.seqIdx != -1 )
					{
						newConsensusLen += leftExtend.seqEnd - leftExtend.seqStart + 1 ;
						if ( matePrevNextType[ chain[j] ].a == 2 )
						{
							struct _pair ng ;
							ng.a = newConsensusLen ;
							ng.b = gapN ;
							gapPos.PushBack( ng ) ;

							newConsensusLen += gapN ;
						}
						leftMostExtend = leftExtend ;
					}
				}
				offset[j] = newConsensusLen ;
				if ( leftExtend.seqIdx != -1 )
					range[j].a = leftExtend.readStart ;
				else
					range[j].a = 0 ;

				if ( rightExtend.seqIdx != -1 )
					range[j].b = rightExtend.readEnd ;
				else
					range[j].b = seqs[ chain[j] ].consensusLen - 1 ;
				origRangeB[j] = range[j].b ;

				if ( j < chainSize - 1 )
				{
					// For the middle ones, we only keep the left overlap.
					/*struct _overlap tmp ;
					GetExtendSeqCoord( chain[j + 1], prevAdj[ chain[j + 1] ][ matePrevNext[ chain[j + 1] ].a ], 
						-1, branchAdj, tmp ) ;
					
					if ( tmp.seqIdx != -1 )
					{
						range[j].b = tmp.seqEnd ;
						if ( range[j].b < range[j].a )
							range[j].b = range[j].a - 1 ;
					}
					else
					{
						// should never get here.
						chainSize = j + 1 ;
					}*/
					range[j].b -= rightExtend.matchCnt ;
					if ( range[j].b < range[j].a )
						range[j].b = range[j].a - 1 ;
				} 

				newConsensusLen += range[j].b - range[j].a + 1 ;

				if ( matePrevNextType[ chain[j] ].b == 2 )
				{
					struct _pair ng ;
					ng.a = newConsensusLen ;
					ng.b = gapN ;
					gapPos.PushBack( ng ) ;
					
					newConsensusLen += gapN ;
				}
				if ( j == chainSize - 1 && rightExtend.seqIdx != -1 )
				{
					newConsensusLen += rightExtend.seqEnd - rightExtend.seqStart + 1 ;
					rightMostExtend = rightExtend ;
				}
			}
			if ( newConsensusLen == seqs[i].consensusLen ) // no extension.
			{
				continue ;
			}
			//printf( "%d %d %d\n", chainSize, seqs[i].consensusLen, newConsensusLen ) ;	
			char *newConsensus = ( char * )malloc( sizeof( char ) * ( newConsensusLen + 1 ) ) ;
			// Put in the seqs.
			if ( leftMostExtend.seqIdx != -1 )
				memcpy( newConsensus, seqs[ leftMostExtend.seqIdx ].consensus + leftMostExtend.seqStart,
					leftMostExtend.seqEnd - leftMostExtend.seqStart + 1 ) ;
			for ( j = 0 ; j < chainSize ; ++j )
			{
				memcpy( newConsensus + offset[j], seqs[ chain[j] ].consensus + range[j].a,
					range[j].b - range[j].a + 1 ) ;
			}
			if ( rightMostExtend.seqIdx != -1 )
			{
				int lastOffset =  offset[j - 1] + range[j - 1].b - range[j - 1].a + 1 ;
				if ( matePrevNextType[ chain[chainSize - 1] ].b == 2 )
					lastOffset += gapN ;
				memcpy( newConsensus + lastOffset,
					seqs[ rightMostExtend.seqIdx ].consensus + rightMostExtend.seqStart,
					rightMostExtend.seqEnd - rightMostExtend.seqStart + 1 ) ;
			}
			// Fill the gaps with Ns
			int gapCnt = gapPos.Size() ;
			for ( j = 0 ; j < gapCnt ; ++j )
			{
				int l ;
				for ( l = gapPos[j].a ; l < gapPos[j].a + gapPos[j].b ; ++l )
				{
					newConsensus[l] = 'N' ;
				}
			}
			
			newConsensus[newConsensusLen] = '\0' ;

			struct _seqWrapper ns ;
			ns.isRef = false ;
			ns.consensus = newConsensus ;
			ns.consensusLen = newConsensusLen ;
			ns.barcode = seqs[i].barcode ;
			
			ns.name = strdup( seqs[i].name ) ;
			ns.posWeight.ExpandTo( newConsensusLen ) ;
			ns.posWeight.SetZero( 0, newConsensusLen ) ;

			// Assign the posWeight of the core part.	
			int newSeqIdx = seqs.size() ;
			for ( j = 0 ; j < chainSize ; ++j )
			{
				int l ;
				for ( l = range[j].a ; l <= origRangeB[j] && offset[j] + l - range[j].a < newConsensusLen ; ++l )
				{
					ns.posWeight[ offset[j] + l - range[j].a ] += seqs[ chain[j] ].posWeight[l] ;
				}
		
				seqs[ chain[j] ].info[0].b = newSeqIdx ; // They are mapped to this new index.
			}
			ns.info[0].a = ns.info[1].b = newSeqIdx ;
			ns.info[1].a = chain[0] ;
			ns.info[1].b = chain[ chainSize - 1 ] ;

			// Adjust the posWeight for the overhang part
			if ( leftMostExtend.seqIdx != -1 )
			{
				int from = leftMostExtend.seqIdx ;
				int to = chain[0] ;
				
				int size = nextAdj[from].size() ;
				for ( j = 0 ; j < size ; ++j )
					if ( nextAdj[from][j].seqIdx == to )
						break ;
				
				if ( j < size )
				{
					// Should always get here
					SimpleVector<int> &readIdx = *( nextAdj[from][j].info ) ; 
					size = readIdx.Size() ;
					int l ;
					for ( l = 0 ; l < size ; ++l )
					{
						int ridx ;
						if ( reads[ readIdx[l] ].overlap.seqIdx == from )
							ridx = readIdx[l] ;
						else
							ridx = readIdx[l] + 1 ;
						
						if ( reads[ ridx ].overlap.seqEnd > leftMostExtend.seqEnd + leftMostExtend.matchCnt )
							continue ;

						int m, rm ;
						for ( m = reads[ ridx ].overlap.seqStart, rm = 0 ; 
							m <= reads[ ridx ].overlap.seqEnd ; ++m, ++rm )
						{
							if ( reads[ ridx ].read[rm] != 'N' )
							{
								if ( m - leftMostExtend.seqStart >= 0 &&
									m - leftMostExtend.seqStart < newConsensusLen )
									++ns.posWeight[m - leftMostExtend.seqStart].
										count[ nucToNum[ reads[ ridx ].read[rm] - 'A' ] ] ;
								
								if ( shiftSeq[from].b + m >= 0 && 
									shiftSeq[from].b + m < seqs[ shiftSeq[from].a ].consensusLen )
									--seqs[ shiftSeq[from].a ].posWeight[ shiftSeq[from].b + m ].
												count[ nucToNum[ reads[ ridx ].read[rm] - 'A' ] ] ;
							}
						}
					}
				}
			}

			if ( rightMostExtend.seqIdx != -1 )
			{
				int from = chain[ chainSize - 1 ] ;
				int to = rightMostExtend.seqIdx ;
				
				int size = nextAdj[from].size() ;
				for ( j = 0 ; j < size ; ++j )
					if ( nextAdj[from][j].seqIdx == to )
						break ;
				
				if ( j < size )
				{
					// Should always get here
					SimpleVector<int> &readIdx = *( nextAdj[from][j].info ) ; 
					size = readIdx.Size() ;
					int l ;
					int lastOffset = offset[chainSize - 1] + range[chainSize - 1].b - range[chainSize - 1].a + 1 ; 
					if ( matePrevNextType[ chain[ chainSize - 1] ].b == 2 )
						lastOffset += gapN ;
					for ( l = 0 ; l < size ; ++l )
					{
						int ridx ;
						if ( reads[ readIdx[l] ].overlap.seqIdx == from )
							ridx = readIdx[l] + 1 ;
						else
							ridx = readIdx[l] ;
						
						if ( reads[ ridx ].overlap.seqStart < rightMostExtend.seqStart - rightMostExtend.matchCnt )
							continue ;

						int m ;
						int rm ;
						char *s = strdup( reads[ ridx ].read ) ;
						if ( reads[ ridx ].overlap.strand == -1 )
						{
							// should always get here
							ReverseComplement( s, reads[ ridx ].read, strlen( reads[ ridx ].read ) ) ;
						}
						for ( m = reads[ ridx ].overlap.seqStart, rm = 0 ; 
							m <= reads[ ridx ].overlap.seqEnd ; ++m, ++rm )
						{
							if ( s[rm] != 'N' )
							{
								int adjustM = m - rightMostExtend.seqStart + lastOffset ;  
								//printf( "%d %d %d. %c %c\n", m, rm, adjustM, ns.consensus[adjustM], reads[ ridx ].read[rm] ) ;
								if ( adjustM >= 0 && adjustM < newConsensusLen )
									++ns.posWeight[adjustM].count[ nucToNum[ s[rm] - 'A' ] ] ;
								if ( shiftSeq[to].b + m >= 0 && 
									shiftSeq[to].b + m < seqs[ shiftSeq[to].a ].consensusLen )
									--seqs[ shiftSeq[to].a ].posWeight[ shiftSeq[to].b + m].
												count[ nucToNum[ s[rm] - 'A' ] ] ;
							}
						}
						free( s ) ;
					}
				}

			}
			
			// For other not updated region, just assign a number there.
			for ( j = 0 ; j < newConsensusLen ; ++j )
				if ( newConsensus[j] != 'N' && 
					ns.posWeight[j].Sum() == 0 )
					ns.posWeight[j].count[ nucToNum[ newConsensus[j] - 'A' ] ] = 1 ;

			// Update the shift information
			for ( j = 0 ; j < chainSize ; ++j )
			{
				shiftSeq[ chain[j] ].a = seqs.size() ;
				shiftSeq[ chain[j] ].b = offset[j] - range[j].a ;
			}

#ifdef DEBUG	
			if ( leftMostExtend.seqIdx != -1 )
				printf( "left 0: %d %s\n", leftMostExtend.seqIdx, seqs[ leftMostExtend.seqIdx ].consensus ) ;

			for ( j = 0 ; j < chainSize ; ++j )
				printf( "chain %d: %d %s\n", j + 1, chain[j], seqs[ chain[j] ].consensus ) ;
			if ( rightMostExtend.seqIdx != -1 )
				printf( "right %d: %d %s\n", j + 1, rightMostExtend.seqIdx, seqs[ rightMostExtend.seqIdx ].consensus ) ;
			printf( "%d new %s\n", i, newConsensus) ;
			fflush( stdout ) ;
#endif 

			ns.minLeftExtAnchor = ns.minRightExtAnchor = 0 ;
			seqs.push_back( ns ) ;
			toRemoveSeqIdx.PushBack( i ) ;
		}
		
		int size ;
		delete[] branchAdj ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			size = nextAdj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( nextAdj[i][j].info != NULL )
					delete nextAdj[i][j].info ;
			}
			
			size = prevAdj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( prevAdj[i][j].info != NULL )
					delete prevAdj[i][j].info ;
			}

		}
		delete[] nextAdj ;
		delete[] prevAdj ;
		
		size = toRemoveSeqIdx.Size() ;
		for ( i = 0 ; i < size ; ++i )
			ReleaseSeq( toRemoveSeqIdx[i] ) ;	
		
		// Recompute the seq id the read is assigned to.
		/*for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( reads[i].seqIdx != -1 )
			{
				reads[i].seqIdx = shiftSeq[ reads[i].seqIdx ].a ;
				reads[i].seqStart += shiftSeq[ reads[i].seqIdx ].b ;
				reads[i].seqEnd += shiftSeq[ reads[i].seqIdx ].b ;
			}
		}*/

		// Recompute the posweight that becomes negative 
		// since the alignment might be changed when adding and after adding states
		seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL )
				continue ;

			int s = -1, e = 0 ;
			for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
			{
				int sum = 0 ;
				for ( k = 0 ; k < 4 ; ++k )
				{
					if ( seqs[i].posWeight[j].count[k] < 0 )
						seqs[i].posWeight[j].count[k] = 0 ;
					sum += seqs[i].posWeight[j].count[k] ;
				}	
				if ( sum > 0 )
				{
					if ( s == -1 )
						s = j ;
					e = j ;
				}

				if ( sum == 0 && seqs[i].consensus[j] != 'N' )
				{
					seqs[i].posWeight[j].count[ nucToNum[ seqs[i].consensus[j] - 'A' ] ] = 1 ;
				}
			}

			if ( s + 10 > e ) // too short
			{
				ReleaseSeq( i ) ;
				continue ;
			}
			
			if ( s > 0 || e < seqs[i].consensusLen - 1 ) // trim.
			{
				for ( j = s ; j <= e ; ++j )
				{
					seqs[i].posWeight[j - s] = seqs[i].posWeight[j] ;
					seqs[i].consensus[j - s] = seqs[i].consensus[j] ;
				}
				seqs[i].consensusLen = e - s + 1 ;
				seqs[i].consensus[ e - s + 1 ] = '\0' ;
			}
		}

		Clean( true ) ;
		novelSeqSimilarity = backupNovelSeqSimilarity ;
	}

	void Output( FILE *fp, std::vector<std::string> *barcodeIntToStr = NULL )
	{
		int i, j, k ;
		int size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;
			
			if ( barcodeIntToStr == NULL || seqs[i].barcode == -1 || seqs[i].barcode >= barcodeIntToStr->size() )
				fprintf( fp, ">assemble%d %s\n%s\n", i, seqs[i].name, seqs[i].consensus ) ;
			else
				fprintf( fp, ">%s_%d %s\n%s\n", barcodeIntToStr->at( seqs[i].barcode ).c_str(), i, seqs[i].name, 
							seqs[i].consensus ) ;
			
			for ( k = 0 ; k < 4 ; ++k )
			{
				for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
					fprintf( fp, "%d ", seqs[i].posWeight[j].count[k] ) ;
				fprintf( fp, "\n" ) ;
			}
		}
	}

	void OutputRef( FILE *fp )
	{
		int i, j, k ;
		int size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( !seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;

			fprintf( fp, ">%s %d %d %d %d %d\n%s\n", seqs[i].name, 
					seqs[i].info[0].a, seqs[i].info[0].b, 
					seqs[i].info[1].a, seqs[i].info[1].b, 
					seqs[i].info[2].a, seqs[i].consensus ) ;
		}
	}

	char *GetSeqName( int seqIdx )
	{
		return seqs[ seqIdx ].name ;
	}
	
	char *GetSeqConsensus( int seqIdx )
	{
		return seqs[ seqIdx ].consensus ;
	}

	void SetSeqConsensus( int seqIdx, char *nc )
	{
		int len = strlen( nc ) ;
		struct _seqWrapper &seq = seqs[ seqIdx ] ;
		free( seq.consensus ) ;
		seq.consensus = strdup( nc ) ;
		seq.consensusLen = len ;	
		seq.posWeight.ExpandTo( len ) ;
		seq.posWeight.SetZero( 0, len ) ;
		UpdatePosWeightFromRead( seq.posWeight, 0, nc ) ;
	}

	int GetSeqConsensusLen( int seqIdx )
	{
		return seqs[ seqIdx ].consensusLen ;
	}

	int GetSeqWeightSum( int seqIdx )
	{
		int ret = 0 ;
		int i ;
		for ( i = 0 ; i < seqs[ seqIdx ].consensusLen ; ++i )
			ret += seqs[seqIdx].posWeight[i].Sum() ;
		return ret ;
	}

	int GetConsensusWeightSumRange( int seqIdx, int start, int end )
	{
		int ret = 0 ;
		int i ;
		for ( i = start ; i <= end ; ++i )
			ret += seqs[seqIdx].posWeight[i].count[ nucToNum[seqs[seqIdx].consensus[i] - 'A' ] ] ;
		return ret ;
	}
	
	void SetIsLongSeqSet( bool in )
	{
		isLongSeqSet = in ;
	}

	void SetBarcodeFromSeqName(std::map<std::string, int> &barcodeStrToInt)
	{
		int seqCnt = seqs.size() ;
		int i, j ;
		char *buffer = new char[10001] ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int len = strlen( seqs[i].name ) ;
			for ( j = len - 1 ; j >= 0 && seqs[i].name[j] != '_' ; --j )	
				;
			strcpy( buffer, seqs[i].name ) ;	
			if ( j >= 0 )
				buffer[j] = '\0' ;

			std::string s(buffer) ;
			int barcode = -1 ;
			if ( barcodeStrToInt.find( s ) != barcodeStrToInt.end() )
				barcode = barcodeStrToInt[s] ;
			else
			{
				barcode = barcodeStrToInt.size() ;
				barcodeStrToInt[s] = barcode ;
			}
			seqs[i].barcode = barcode ;
		}
		delete[] buffer ;
	}
} ;


#endif
