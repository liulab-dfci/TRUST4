// The data structure holds the set of sequences (can be "assembled" from several reads)
#ifndef _MOURISL_SEQSET_HEADER
#define _MOURISL_SEQSET_HEADER

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <queue>

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
} ;

struct _hit
{
	struct _indexInfo indexHit ;
	int readOffset ;
	int strand ; // -1: different strand, 1: same strand. When strand==-1, the readOffset is the offset in the rcSeq.
	
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

	SimpleVector<struct _pair> hitCoords ;

	bool operator<( const struct _overlap &b ) const
	{
		// The overlap with more matched bases should come first.
		if ( matchCnt > b.matchCnt + 2 || matchCnt < b.matchCnt - 2 )
			return matchCnt > b.matchCnt ;
		else if ( similarity != b.similarity )
			return similarity > b.similarity ; 
		else if ( readEnd - readStart != b.readEnd - b.readStart )
			return readEnd - readStart > b.readEnd - b.readStart ;
		else if ( seqIdx != b.seqIdx )
			return seqIdx < b.seqIdx ;
		else if ( strand !=  b.strand )
			return strand < b.strand ;
		else if ( readStart != b.readStart )
			return readStart < b.readStart ;
		else if ( readEnd != b.readEnd )
			return readEnd < b.readEnd ;
		else if ( seqStart != b.seqStart )
			return seqStart < b.seqStart ;
		else 
			return seqEnd < b.seqEnd ; 

		return false ;
	}
} ;

class SeqSet
{
private:
	std::vector<struct _seqWrapper> seqs ;
	KmerIndex seqIndex ;
	int kmerLength ;
	int minHitRequired ;
	int radius ;
	
	// Some threshold
	double novelSeqSimilarity ;
	double refSeqSimilarity ;
	double repeatSimilarity ; // e.g., the repeat when building the branch graph.

	struct _overlap prevAddInfo ; 

	static bool CompSortPairBInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.b < p2.b ;
	}
	
	static bool CompSortPairAInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.a < p2.a ;
	}

	static bool CompSortOverlapsOnReadCoord( const struct _overlap &a, const struct _overlap &b )
	{
		return a.readStart < b.readStart ; 
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

	void Reverse( char *r, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
			r[i] = seq[len - 1 - i] ;  
		r[i] = '\0' ;
	}
	
	bool IsPosWeightCompatible( const struct _posWeight &a, const struct _posWeight &b )
	{
		int sumA = a.count[0] + a.count[1] + a.count[2] + a.count[3] ;
		int sumB = b.count[0] + b.count[1] + b.count[2] + b.count[3] ;
	
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
		for ( i = 0 ; i < 4 ; ++i )
			if ( cnt[i] <= 2 )
				++lowCnt ;
		if ( lowCnt >= 2 )
			return true ;
		return false ;
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
	int GetOverlapsFromHits( SimpleVector<struct _hit> &hits, int hitLenRequired, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int hitSize = hits.Size() ;
		
		SimpleVector<struct _pair> hitCoordDiff ;
		hitCoordDiff.Reserve( hitSize ) ;
		SimpleVector<struct _pair> concordantHitCoord ;
		SimpleVector<struct _pair> hitCoordLIS ;
		SimpleVector<struct _hit> finalHits ;

		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
					break ;
			//[i,j) holds the hits onto the same seq on the same strand.	
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
				nh.b = hits[k].readOffset - hits[k].indexHit.offset ;
				hitCoordDiff.PushBack( nh ) ;
			}
			std::sort( hitCoordDiff.BeginAddress(), hitCoordDiff.EndAddress(), CompSortPairBInc ) ;

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
					int diff = hitCoordDiff[e].b - hitCoordDiff[e - 1].b ;
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

				// [s, e) holds the candidate in the array of hitCoordDiff 
				concordantHitCoord.Clear() ;
				for ( k = s ; k < e ; ++k )
				{
					struct _pair nh ;
					int hitId = hitCoordDiff[k].a ; 
					nh.a = hits[ hitId ].readOffset ;
					nh.b = hits[ hitId ].indexHit.offset ;
					concordantHitCoord.PushBack( nh ) ;
				}

				std::sort( concordantHitCoord.BeginAddress(), concordantHitCoord.EndAddress(), CompSortPairBInc ) ;
				//for ( k = 0 ; k < e - s ; ++k )	
				//	printf( "%d (%d-%d): %d %d %d\n", i, s, e, hits[i].indexHit.idx, concordantHitCoord[k].a, concordantHitCoord[k].b ) ;


				// Compute the longest increasing subsequence.
				hitCoordLIS.Clear() ;
				int lisSize = LongestIncreasingSubsequence( concordantHitCoord, hitCoordLIS ) ; 
				if ( lisSize * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}

				// Rebuild the hits.
				finalHits.Clear() ;
				for ( k = 0 ; k < lisSize ; ++k )
				{
					struct _hit nh = hits[i];
					nh.readOffset = hitCoordLIS[k].a ;
					nh.indexHit.offset = hitCoordLIS[k].b ;
					//printf( "%d: %d %d %d %d\n", i, nh.readOffset, nh.indexHit.idx, nh.indexHit.offset, nh.strand ) ;
					finalHits.PushBack( nh ) ;
				}

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
				no.matchCnt = hitLen ;
				no.similarity = 0 ;

				if ( hitLen * 2 < no.seqEnd - no.seqStart + 1 )
				{
					s = e ; 
					continue ;
				}
				no.hitCoords.Reserve( lisSize ) ;
				for ( k = 0 ; k < lisSize ; ++k )
				{
					struct _pair nh ;
					nh.a = finalHits[k].readOffset ;
					nh.b = finalHits[k].indexHit.offset ;
					no.hitCoords.PushBack( nh ) ;
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

		GetOverlapsFromHits( VJhits, 17, overlaps ) ;
		
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
			overlaps.clear() ;
			return 0 ;
		}

		std::vector<struct _overlap> ret ;
		ret.push_back( overlaps[ tagi ] ) ;
		ret.push_back( overlaps[ tagj ] ) ;
		
		overlaps = ret ;

		return 2 ;
	}

	// Extend the overlap to include the overhang parts and filter the overlaps if the overhang does not match well.
	// return: whether this is a valid extension or not
	int ExtendOverlap( char *r, int len, struct _seqWrapper &seq, char *align, struct _overlap &overlap, struct _overlap &extendedOverlap )
	{
		// Check whether the overhang part is compatible with each other or not.
		// Extension to 5'-end ( left end )
		int matchCnt, mismatchCnt, indelCnt ;
		int leftOverhangSize = MIN( overlap.readStart, overlap.seqStart ) ;
		int ret = 1 ;

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
		// Extension to 3'-end ( right end )
		int rightOverhangSize = MIN( len - 1 - overlap.readEnd, seq.consensusLen - 1 - overlap.seqEnd ) ;

		//AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqEnd + 1, 
		AlignAlgo::GlobalAlignment_PosWeight( seq.posWeight.BeginAddress() + overlap.seqEnd + 1, 
				rightOverhangSize,
				r + overlap.readEnd + 1, rightOverhangSize, align ) ;
		GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > 0 )
		{
			rightOverhangSize = 0 ;
			ret = 0 ;
		}

		int mismatchThreshold = 2 ;
		if ( leftOverhangSize >= 2 )
			++mismatchThreshold ;
		if ( rightOverhangSize >= 2 )
			++mismatchThreshold ;
		
		if ( mismatchCnt > mismatchThreshold && (double)mismatchCnt / ( leftOverhangSize + rightOverhangSize ) > 1.5 / kmerLength ) 
			ret = 0 ;

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

		/*if ( !seqs[ extendedOverlap.seqIdx ].isRef && 
			extendedOverlap.readEnd - extendedOverlap.readStart + 1 - extendedOverlap.matchCnt / 2 >= 5 )
		{
			// Exceed the number of mismatch allowed
			extendedOverlap = overlap ;
			ret = 0 ;
		}*/
		return ret ;
	}

	// Obtain the overlaps, each overlap further contains the hits induce the overlap. 
	// Return: the number of overlaps.
	int GetOverlapsFromRead( char *read, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return -1 ;

		SimpleVector<struct _hit> hits ;		    	

		KmerCode kmerCode( kmerLength ) ;
		KmerCode prevKmerCode( kmerLength ) ;

		// Locate the hits from the same-strand case.
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;
		
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			if ( i == kmerLength || !prevKmerCode.IsEqual( kmerCode ) )
			{
				SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

				int size = indexHit.Size() ;

				//if ( size >= 40 )
				//	continue ;

				for ( j = 0 ; j < size ; ++j )
				{
					struct _hit nh ;
					nh.indexHit = indexHit[j] ;
					nh.readOffset = i - kmerLength + 1 ;
					nh.strand = 1 ;
					hits.PushBack( nh ) ;
				}
			}

			prevKmerCode = kmerCode ;
		}
		// Locate the hits from the opposite-strand case.
		char *rcRead =  new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( rcRead[i] ) ;
		
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( rcRead[i] ) ;
			if ( i == kmerLength || !prevKmerCode.IsEqual( kmerCode ) )
			{
				SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

				int size = indexHit.Size() ;

				//if ( size >= 40 )
				//	continue ;

				for ( j = 0 ; j < size ; ++j )
				{
					struct _hit nh ;
					nh.indexHit = indexHit[j] ;
					nh.readOffset = i - kmerLength + 1 ;
					nh.strand = -1 ;
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
		delete[] rcRead ;

		// Find the overlaps.
		std::sort( hits.BeginAddress(), hits.EndAddress() ) ;
		//for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
		//	printf( "- %d %d %d %d\n", it->readOffset, it->indexHit.idx, it->indexHit.offset, it->strand ) ;

		int overlapCnt = GetOverlapsFromHits( hits, 31, overlaps ) ;

		/*for ( i = 0 ; i < overlapCnt ; ++i )
			printf( "%d: %d %s %d. %d %d %d %d\n", i, overlaps[i].seqIdx,seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ;*/ 
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
		k = 0 ;
		for ( i = 1 ; i < overlapCnt ; ++i )
		{
			if ( overlaps[i].strand != overlaps[i - 1].strand )
			{
				overlaps.clear() ;
				return 0 ;
			}
		}

		rcRead = new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		

		// Compute similarity overlaps
		//std::sort( overlaps.begin(), overlaps.end() ) ;
		int firstRef = -1 ;

		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *r ;
			if ( overlaps[i].strand == 1 )
				r = read ;
			else
				r = rcRead ;

			SimpleVector<struct _pair> &hitCoords = overlaps[i].hitCoords ; 	
			
			if ( seqs[ overlaps[i].seqIdx ].isRef )
			{
				if ( firstRef == -1 )
					firstRef = i ;
				/*else
				{
					if ( overlaps[i].matchCnt < 0.95 * overlaps[ firstRef ].matchCnt )
					{
						overlaps[i].similarity = 0 ;  // No need to look into this.
						continue ;
					}
				}*/
			}

			int hitCnt = hitCoords.Size() ;
			int matchCnt = 0, mismatchCnt = 0, indelCnt = 0  ;
			double similarity = 1 ;
			
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
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;
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
					}
				}
				else
				{
					if ( radius == 0 || seqs[ overlaps[i].seqIdx ].isRef )
					{
						similarity = 0 ;
						break ;
					}

					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 < hitCoords[j].b )
					{
						matchCnt += ( hitCoords[j].a - hitCoords[j - 1].a ) + kmerLength ;
						// Make the two kmer hit match on coordinate.
						indelCnt += ( hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) + 
							( hitCoords[j].a + kmerLength - hitCoords[j - 1].a )  ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 < hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						matchCnt += kmerLength + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ( hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) +
							( hitCoords[j].b + kmerLength - hitCoords[j - 1].b ) ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a &&
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						matchCnt += ( hitCoords[j].a - hitCoords[j - 1].a ) + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ABS( ( hitCoords[j].a - hitCoords[j].b ) - 
							( hitCoords[j - 1].a - hitCoords[j - 1].b ) ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ;
						 
						if ( seqs[ overlaps[i].seqIdx ].isRef )
						{
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
			
			//printf( "%d %d %d %lf\n", matchCnt, overlaps[i].seqEnd - overlaps[i].seqStart + 1, overlaps[i].readEnd - overlaps[i].readStart + 1, similarity ) ;
			overlaps[i].matchCnt = matchCnt ;
			if ( similarity == 1 )
				overlaps[i].similarity = (double)matchCnt / ( overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 
								overlaps[i].readEnd - overlaps[i].readStart + 1 ) ;
			else
				overlaps[i].similarity = 0 ;
			
			if ( IsOverlapLowComplex( r, overlaps[i]) )
				overlaps[i].similarity = 0 ;
			
			/*if ( overlaps[i].similarity > 1 )
			{
				printf( "%d: %d %d %d %d\n", matchCnt, overlaps[i].readStart, overlaps[i].readEnd, 
							overlaps[i].seqStart, overlaps[i].seqEnd ) ;
			}
			assert( overlaps[i].similarity <= 1 ) ;*/
		} // for i
		delete[] rcRead ;

		// Release the memory for hitCoords.
		for ( i = 0 ; i < overlapCnt ; ++i )
			overlaps[i].hitCoords.Release() ;

		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < refSeqSimilarity )
				continue ;
			else if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < novelSeqSimilarity )
				continue ;

			//printf( "%lf\n", overlaps[i].similarity ) ;
			overlaps[k] = overlaps[i] ;
			++k ;
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		

		return overlapCnt ;
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

			overlapCnt = GetOverlapsFromRead( seqs[i].consensus, overlaps ) ;
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( overlaps[j].strand == -1 ) // Note that all the seqs are from 5'->3' on its strand.
					continue ;

				if ( i == overlaps[j].seqIdx  )
					continue ;
				struct _overlap extendedOverlap ;

				if ( ExtendOverlap( seqs[i].consensus, seqs[i].consensusLen, seqs[ overlaps[j].seqIdx ], 
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

	int BuildBranchGraph( std::vector<struct _overlap> *adj )
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

			overlapCnt = GetOverlapsFromRead( seqs[i].consensus, overlaps ) ;
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( overlaps[j].strand == -1 ) // Note that all the seqs are from 5'->3' on its strand.
					continue ;

				if ( i == overlaps[j].seqIdx  )
					continue ;
				struct _overlap extendedOverlap ;
				int addMatchCnt = 0 ;				
				// Locally extend the overlap. In this case, the overlap actually just means share.
				// Allow "radius" overhang.
				int seqIdx = overlaps[j].seqIdx ;
				int a, b, k ;
				int matchCnt = 0 ;
				int rightExtend = 0;
				int rightExtendMatchCnt = 0 ; 
				for ( k = 1, a = overlaps[j].readEnd + 1, b = overlaps[j].seqEnd + 1 ; 
					a < seqs[i].consensusLen && b < seqs[ seqIdx ].consensusLen ; ++a, ++b, ++k )
				{
					if ( IsPosWeightCompatible( seqs[i].posWeight[a], seqs[ seqIdx ].posWeight[b] ) )
						++matchCnt ;

					if ( matchCnt > k * 0.75 )
					{
						rightExtendMatchCnt = 2 * matchCnt ;
						rightExtend = k ;
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
						++matchCnt ;

					if ( matchCnt > k * 0.75 )
					{
						leftExtendMatchCnt = 2 * matchCnt ;
						leftExtend = k ;
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
				printf( "branch %d %d: %d %d %d %d %d %lf\n", i, j, extendedOverlap.seqIdx, 
						extendedOverlap.readStart, extendedOverlap.readEnd,
						extendedOverlap.seqStart, extendedOverlap.seqEnd,
						extendedOverlap.similarity ) ;
				
				if ( extendedOverlap.similarity >= repeatSimilarity )
					adj[i].push_back( extendedOverlap ) ;
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
			std::sort( adj[i].begin(), adj[i].end() ) ;
		delete[] align ;
		return 1 ;
	}
public:
	SeqSet( int kl ) 
	{
		kmerLength = kl ;
		minHitRequired = 3 ;
		radius = 10 ;
		
		novelSeqSimilarity = 0.9 ;
		refSeqSimilarity = 0.75 ; 
		repeatSimilarity = 0.95 ; 

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
	
	// Input some baseline sequence to match against.
	void InputRefFa( char *filename ) 
	{
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		
		KmerCode kmerCode( kmerLength ) ;
		while ( fa.Next() )
		{
			// Insert the kmers 
			struct _seqWrapper ns ;
			ns.name = strdup( fa.id ) ;
			ns.isRef = true ;

			int id = seqs.size() ;
			seqs.push_back( ns ) ;

			struct _seqWrapper &sw = seqs[id] ;
			int seqLen = strlen( fa.seq ) ;
			sw.consensus = strdup( fa.seq ) ;	
			sw.consensusLen = strlen( fa.seq );	
			seqIndex.BuildIndexFromRead( kmerCode, fa.seq, seqLen, id ) ;
		}
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

			int id = seqs.size() ;
			ns.consensus = strdup( in.seqs[i].consensus ) ;
			ns.consensusLen = in.seqs[i].consensusLen ;
			ns.posWeight = in.seqs[i].posWeight ;
			seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, ns.consensusLen, id ) ;
			seqs.push_back( ns ) ;
		}
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


	

	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set. 
	//	   -1: not add. -2: only overlapped with novel seq and could not be extended.
	int AddRead( char *read )
	{
		//printf( "%s\n", seq ) ;
		int i, j, k ;
		int len = strlen( read ) ;

		SetPrevAddInfo( -1, -1, -1, -1, -1, 0 ) ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;
		
		overlapCnt = GetOverlapsFromRead( read, overlaps ) ;
				
		if ( overlapCnt <= 0 )
			return -1 ;

		std::sort( overlaps.begin(), overlaps.end() ) ;

#ifdef DEBUG
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			printf( "%d: %d %d %s. %d. %d %d %d %d. %lf.\n", i, overlaps[i].seqIdx, seqs[ overlaps[i].seqIdx ].consensusLen, 
					seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, 
					overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd, 
					overlaps[i].similarity ) ; 
			//if ( !seqs[ overlaps[i].seqIdx ].isRef )
			//	printf( " %s\n",seqs[  overlaps[i].seqIdx ].consensus ) ;
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

			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				for ( j = 0 ; j < k ; ++j )
				{
					// If is contained in extended sequence.
					if ( overlaps[i].readStart >= extendedOverlaps[j].readStart - radius  
						&& overlaps[i].readEnd <= extendedOverlaps[j].readEnd + radius )
						break ;
					// Some extended is a subset of this one.
					if ( extendedOverlaps[j].readStart >= overlaps[i].readStart - radius 
						&& extendedOverlaps[j].readEnd <= overlaps[i].readEnd + radius )
						break ;
				}
				
				struct _seqWrapper &seq = seqs[ overlaps[i].seqIdx ] ; 
				if ( j < k || seq.isRef )
					continue ;

				// Only extend the novel seqs.
				if ( ExtendOverlap( r, len, seq, align, overlaps[i], extendedOverlaps[k] ) == 1 )
				{
					// Double check whether there is subset relationship.
					for ( j = 0 ; j < k ; ++j )
					{
						if ( extendedOverlaps[k].readStart >= extendedOverlaps[j].readStart - radius  
								&& extendedOverlaps[k].readEnd <= extendedOverlaps[j].readEnd + radius )
							break ;
						// Some extended is a subset of this one.
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

						if ( extendedOverlaps[k].readStart >= overlaps[j].readStart &&
							extendedOverlaps[k].readEnd <= overlaps[j].readEnd )
						{
							if ( extendedOverlaps[k].readStart > 0 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							if ( extendedOverlaps[k].readEnd < len - 1 )
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
						if ( extendedOverlaps[k].readStart >= failedExtendedOverlaps[j].readStart &&
							extendedOverlaps[k].readEnd <= failedExtendedOverlaps[j].readEnd )
						{
							if ( extendedOverlaps[k].readStart > 0 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							if ( extendedOverlaps[k].readEnd < len - 1 )
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
				
					if ( ExtendOverlap( r, len, seq, align, overlaps[i], extendedOverlaps[k] ) == 1 )
						++k ;
				}

				if ( k > 2 )
				{
					k = 1 ;
				}
				else if ( k == 2 )
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

#ifdef DEBUG
			for ( i = 0 ; i < k ; ++i )
				printf( "extended %d: %d %s. %d. %d %d %d %d %lf\n", i, extendedOverlaps[i].seqIdx, 
						seqs[ extendedOverlaps[i].seqIdx ].name, extendedOverlaps[i].strand, 
						extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, extendedOverlaps[i].seqStart, 
						extendedOverlaps[i].seqEnd, extendedOverlaps[i].similarity ) ; 
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

				// Copy the original consensus in.
				// The earlier seq has higher weight.
				for ( i = eOverlapCnt - 1 ; i >= 0 ; --i )
				{
					memcpy( newConsensus + seqOffset[i], seqs[ extendedOverlaps[i].seqIdx ].consensus,
						seqs[ extendedOverlaps[i].seqIdx ].consensusLen ) ;	
				}

				// Fill in the gaps
				if ( extendedOverlaps[0].readStart > 0 )
					memcpy( newConsensus, r, extendedOverlaps[0].readStart ) ;		
				
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( extendedOverlaps[i].readStart > extendedOverlaps[i - 1].readEnd + 1 )
						memcpy( newConsensus + seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen,
							r + extendedOverlaps[i - 1].readEnd + 1, 
							extendedOverlaps[i].readStart - extendedOverlaps[i - 1].readEnd - 1 ) ;
				}
				

				int newConsensusLen = 0 ;
				if ( extendedOverlaps[i - 1].readEnd < len )
				{
					memcpy( newConsensus + seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen,
							r + extendedOverlaps[i - 1].readEnd + 1,
							len - extendedOverlaps[i - 1].readEnd -1 ) ;
					newConsensusLen = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen 
					                     + ( len - extendedOverlaps[i - 1].readEnd - 1 ) ;
					newConsensus[ newConsensusLen ] = '\0' ;
				}
				else
				{
					newConsensusLen = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen ;
					newConsensus[ newConsensusLen ] = '\0' ;
				}
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
				seqs[ newSeqIdx ].minRightExtAnchor = seqs[ extendedOverlaps[ eOverlapCnt - 1 ].seqIdx ].minRightExtAnchor ;
				
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
						
						seq.posWeight.SetZero( 0, shift ) ;
					}
					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						int start = extendedOverlaps[0].readStart + seq.consensusLen ;
						seq.posWeight.SetZero( start, len - extendedOverlaps[0].readEnd - 1 ) ;
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
		}

		k = 0 ;
		// See whether there is a reference seq match is sequence if it does not match any novel seq.
		int refSeqIdx = -1 ;
		if ( addNew )
		{
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				if ( seqs[ overlaps[i].seqIdx ].isRef )
				{
					refSeqIdx = overlaps[i].seqIdx ;
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
		}

		if ( addNew )
		{
			// A novel sequence
			// Go through the reference to annotate this read.
			for ( i = 0 ; i < overlapCnt ; ++i )					
			{
				// Check whether this overlap is used.
				for ( j = 0 ; j < k ; ++j )
				{
					;
				}
			}
			// TODO: If overlaps on C-gene, it should not be truncated.  
			
			// Add the sequence to SeqSet
			int idx = seqs.size() ;
			struct _seqWrapper ns ;
			
			if ( refSeqIdx >=0 )
				ns.name = strdup( seqs[ refSeqIdx ].name ) ;
			else
				ns.name = strdup( "unknown" ) ;
			ns.consensus = strdup( read ) ;
			ns.consensusLen = strlen( read ) ;
			if ( overlaps[0].strand == -1 )
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

			if ( nucToNum[ seq.consensus[i] - 'A' ] != maxTag )
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
		seqs.resize( k ) ;
	}
	
	// Return:the number of connections made.
	int ExtendSeqFromSeqOverlap()
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int ret = 0 ;
		novelSeqSimilarity = 0.97 ;

		std::vector<struct _overlap> *adj = new std::vector<struct _overlap>[ seqCnt ] ;
		struct _pair *next = new struct _pair[ seqCnt ] ; // a-index, b-the index in adj[i] 
		struct _pair *prev = new struct _pair[ seqCnt ] ;
		int *containedIn = new int[seqCnt] ;
	
		//BuildSeqOverlapGraph( 100, adj ) ;
		BuildBranchGraph( adj ) ;
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
				for ( j = 0 ; j < size ; ++j )
					if ( adj[i][j].seqIdx == containedIn[ adj[i][0].seqIdx ] )
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
				ns.posWeight.SetZero( 0, newConsensusLen - 1 ) ;
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
	
	// Figure out the gene composition for the read. 
	// Return successful or not.
	int AnnotateRead( char *read, int detailLevel, struct _overlap geneOverlap[4], char *buffer )
	{
		int i, j ;
		
		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;
	
		char BT = '\0' ;
		char chain = '\0' ;
		int len = strlen( read ) ;

		geneOverlap[0].seqIdx = geneOverlap[1].seqIdx = geneOverlap[2].seqIdx = geneOverlap[3].seqIdx = -1 ;
		
		sprintf( buffer, "%d", len ) ;
		overlapCnt = GetOverlapsFromRead( read, overlaps ) ;		
		if ( overlapCnt == 0 )
			return 0 ;
		std::sort( overlaps.begin(), overlaps.end() ) ;
		// Get the coverage of the genes.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *name = seqs[ overlaps[i].seqIdx ].name ;
			if ( BT && name[0] != BT )
				continue ;
			BT = name[0] ;
			
			if ( chain && name[2] != chain )
				continue ;
			chain = name[2] ;

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
			//printf( "%d %s %lf: %d %d\n", i, seqs[ overlaps[i].seqIdx].name, overlaps[i].similarity, overlaps[i].readStart, overlaps[i].readEnd ) ;
			if ( geneType < 0 || geneOverlap[ geneType ].seqIdx != -1 )
				continue ;
			
			geneOverlap[ geneType ] = overlaps[i] ;
		}
		
		// Extend overlap
		if ( detailLevel >= 1 )
		{
			char *align = new char[ 2 * len + 2 ] ;
			char *rvr = new char[len + 1] ; 
			for ( i = 0 ; i < 4 ; ++i )
			{
				// Extend right.
				if ( geneOverlap[i].seqIdx == -1 )
					continue ;
				int seqIdx = geneOverlap[i].seqIdx ;				
				AlignAlgo::GlobalAlignment_OneEnd( seqs[ seqIdx ].consensus + geneOverlap[i].seqEnd + 1, seqs[ seqIdx ].consensusLen - geneOverlap[i].seqEnd, read + geneOverlap[i].readEnd + 1, len - geneOverlap[i].readEnd, 0, align ) ;
				
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( align[j] == EDIT_MATCH || align[j] == EDIT_MISMATCH )
					{
						++geneOverlap[i].readEnd ;
						++geneOverlap[i].seqEnd ;
						
						geneOverlap[i].matchCnt += 2 ;
					}
					else if ( radius > 0 )
					{
						if ( align[j] == EDIT_INSERT )
							++geneOverlap[i].readEnd ;
						else if ( align[j] == EDIT_DELETE )
							++geneOverlap[i].seqEnd ;
					}
					else
						break ;
				}

				// Extend left.
				char *rvs = new char[seqs[ seqIdx ].consensusLen ] ;
				Reverse( rvr, read, geneOverlap[i].readStart ) ;
				Reverse( rvs, seqs[seqIdx].consensus, geneOverlap[i].seqStart ) ;
				//rvr[geneOverlap[i].readStart] = '\0' ;
				//rvs[geneOverlap[i].seqStart] = '\0' ;
				
				AlignAlgo::GlobalAlignment_OneEnd( rvs, geneOverlap[i].seqStart, rvr, geneOverlap[i].readStart, 0, align ) ;
				//AlignAlgo::VisualizeAlignment( rvs, geneOverlap[i].readStart, rvr, rvr[geneOverlap[i].readStart], align ) ;
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( align[j] == EDIT_MATCH || align[j] == EDIT_MISMATCH )
					{
						--geneOverlap[i].readStart ;
						--geneOverlap[i].seqStart ;
						geneOverlap[i].matchCnt += 2 ;
					}
					else if ( radius > 0 )
					{
						if ( align[j] == EDIT_INSERT )
							--geneOverlap[i].readStart ;
						else if ( align[j] == EDIT_DELETE )
							--geneOverlap[i].seqStart ;
					}
					else
						break ;
				}
				delete[] rvs ;

				geneOverlap[i].similarity = (double)( geneOverlap[i].matchCnt ) / 
						( geneOverlap[i].seqEnd - geneOverlap[i].seqStart + 1 + 
							geneOverlap[i].readEnd - geneOverlap[i].readStart + 1 ) ;
			}
			delete[] align ;
			delete[] rvr ;
		}

		// Infer CDR1,2,3.
		char *cdr3 = NULL ;
		if ( detailLevel >= 2 )
		{
			// Infer CDR3.
			if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 )
			{
				// The case that we have anchor.
				// Find the motif for anchor.
				if ( geneOverlap[0].readEnd < geneOverlap[2].readStart )
				{
					int s = geneOverlap[0].readEnd + 1 ;
					int e = geneOverlap[2].readStart - 1 ;

					for ( i = s ; i >= 0 ; --i )
					{
						if ( read[i] == 'T' && read[i + 1] == 'G' && read[i + 2] == 'T' )
							break ;
					}
					if ( i >= 0 )
						s = i ;

					for ( i = e ; i < len - 2 ; ++i )
					{
						if ( read[i] == 'T' && read[i + 1] == 'G' && read[i + 2] == 'G' )
							break ;
					}
					if ( i < len - 2 )
						e = i + 2 ;

					cdr3 = new char[e - s + 2  + 1 ] ;
					memcpy( cdr3, read + s, e - s + 1 ) ;
					cdr3[e - s + 1] = '\0' ;
				}
			}
		}

		// Compute the name
		for ( i = 0 ; i < 4 ; ++i )
		{
			if ( geneOverlap[i].seqIdx == -1 )
				continue ;

			int offset = strlen( buffer ) ;
			int seqIdx = geneOverlap[i].seqIdx ;
			sprintf( buffer + offset, " %s(%d):(%d-%d):(%d-%d):%.2lf",
				seqs[ seqIdx ].name, seqs[ seqIdx ].consensusLen,
				geneOverlap[i].readStart, geneOverlap[i].readEnd, 
				geneOverlap[i].seqStart, geneOverlap[i].seqEnd, geneOverlap[i].similarity * 100 ) ;	
		}
		
		sprintf( buffer + strlen( buffer ), " CDR3=%s", cdr3 == NULL ? "null" : cdr3 ) ;

		if ( cdr3 != NULL )
			delete[] cdr3 ;
		return 1 ;
	}
	
	// Use the refSet to annotate current set.
	void Annotate( SeqSet &refSet )
	{
		int i ;
		char *buffer = new char[1024] ;
		int seqCnt = seqs.size() ;
		struct _overlap geneOverlap[4];
		
		for ( i = 0 ; i < seqCnt  ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;
		
			free( seqs[i].name ) ;
			refSet.AnnotateRead( seqs[i].consensus, 2, geneOverlap, buffer ) ;
			seqs[i].name = strdup( buffer ) ;
		}

		delete[] buffer ;
	}
	
	
	void BreakFalseAssembly( std::vector<struct _Read> reads )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		std::vector<struct _overlap> *adj = new std::vector<struct _overlap>[ seqCnt ] ;
		BuildBranchGraph( adj ) ;
		
		int readCnt = reads.size() ;
		
		// Mark whether to trust the branch overlap portion.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			std::vector<struct _overlap> overlaps ;
			int overlapCnt = GetOverlapsFromRead( reads[i].seq, overlaps ) ;
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
				if ( ExtendOverlap( r, len, seqs[ overlaps[j].seqIdx ], align, 
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

				printf( "%d branch %d %d: %d %d %lf\n", i, j, adj[i][j].seqIdx, 
					adj[i][j].readStart, adj[i][j].readEnd, adj[i][j].similarity ) ;
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
			for ( j = 0 ; j < k ; ++j )
				printf( "%d %s breakCoords %d: %d %d\n", i, seqs[i].name, j, breakCoords[j].a, breakCoords[j].b ) ;
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

	void UpdateMateAdjGraph( int from, int fromStart, int fromEnd, int to, int toStart, int toEnd, std::vector<struct _overlap> *mateAdj )
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
			mateAdj[from].push_back( na ) ;
		}
	}
	
	// return: 0: extension. 1: end-to-inside extension. 2: end-to-end extension
	int GetExtendSeqCoord( int from, struct _overlap mateInfo, int direction, 
			std::vector<struct _overlap> *branchAdj, struct _overlap &coord )
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
		
		/*for ( i = 0 ; i < adjCnt ; ++i )
		{
		}*/
		
		coord.seqIdx = to ;
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
			coord.seqEnd = seqs[to].consensusLen - 1 ;

			if ( branchAdj[from][k].seqStart <= overhangSize - 1 )
				ret = 2 ;
		}
		else
		{
			coord.seqStart = 0 ;
			coord.seqEnd = branchAdj[from][k].seqStart - 1 ;

			if ( branchAdj[from][k].seqEnd >= seqs[ branchAdj[from][k].seqIdx ].consensusLen - overhangSize )
				ret = 2 ;
		}

		return ret ;
	}

	// Use this set of reads to extend,rearrange the seq 
	void ExtendSeqFromReads( std::vector<struct _Read> reads )
	{
		int i, j, k ;
		int readCnt = 0 ;
		int seqCnt = seqs.size() ;
		novelSeqSimilarity = 1.0 ;

		int ret = 0 ;

		std::vector<struct _overlap> *branchAdj = new std::vector<struct _overlap>[ seqCnt ] ;

		// Mate-pair support graph. In here, start,end represnts the rough range with read support
		//	and matchCnt represent how many mates support this connection.
		std::vector<struct _overlap> *nextAdj = new std::vector<struct _overlap>[ seqCnt ] ; 
		std::vector<struct _overlap> *prevAdj = new std::vector<struct _overlap>[ seqCnt ] ; 
		readCnt = reads.size() ;
		char *rc ;
		char *mrc ;
		char *r ;
		char *mr ;
		char *align ;
		
		k = 0 ;
		for ( i = 0 ; i < seqCnt ; ++i )
			if ( seqs[i].consensusLen > k )
				k = seqs[i].consensusLen ;
		align = new char[2 * k + 2] ;

		BuildBranchGraph( branchAdj ) ;
		
		// Build the mate adj graph.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			bool paired = false ;
			if ( i < readCnt - 1 && !strcmp( reads[i].id, reads[i + 1].id ) )
				paired = true ; 
			std::vector<struct _overlap> overlaps ;
			int overlapCnt = GetOverlapsFromRead( reads[i].seq, overlaps ) ;
			if ( paired )
			{
				std::vector<struct _overlap> mateOverlaps ;
				int mateOverlapCnt = GetOverlapsFromRead( reads[i + 1].seq, mateOverlaps ) ;
				//printf( "%d %d\n", overlapCnt, mateOverlapCnt ) ;
				//printf( "%d %s\n%d %s\n", overlaps[0].strand, reads[i].seq, mateOverlaps[0].strand, reads[i + 1].seq ) ;

				if ( overlapCnt == 0 || mateOverlapCnt == 0 
					|| overlaps[0].strand == mateOverlaps[0].strand )
				{
					++i ;
					continue ;
				}
				
				std::sort( overlaps.begin(), overlaps.end() ) ;
				std::sort( mateOverlaps.begin(), mateOverlaps.end() ) ;	
				
				int len = strlen( reads[i].seq ) ;
				int mlen = strlen( reads[i + 1].seq ) ;
				rc = new char[len + 1] ;
				mrc = new char[mlen + 1] ;

				ReverseComplement( rc, reads[i].seq, len ) ;
				ReverseComplement( mrc, reads[i + 1].seq, mlen ) ;
				
				r = reads[i].seq ;
				mr = reads[i + 1].seq ;
				if ( overlaps[0].strand == 1 )
					mr = mrc ;
				else
					r = rc ;

				struct _overlap extendedOverlap, mateExtendedOverlap ;
				struct _overlap tmpExtendOverlap ;
				/*for ( j = 0 ; j < overlapCnt ; ++j )
				{
					printf( "+ %d %d: %d %d %d %lf\n", i, j, overlaps[j].seqIdx, overlaps[j].seqStart, overlaps[j].seqEnd, overlaps[j].similarity) ;
				}
				for ( j = 0 ; j < mateOverlapCnt ; ++j )
				{
					printf( "- %d %d: %d %d %d %lf\n", i + 1, j, mateOverlaps[j].seqIdx, mateOverlaps[j].seqStart, mateOverlaps[j].seqEnd, mateOverlaps[j].similarity) ;
				}*/
				int extendCnt = 0 ;
				bool valid = true ;
				for ( j = 0 ; j < overlapCnt ; ++j )
				{
					if ( ExtendOverlap( r, len, seqs[ overlaps[j].seqIdx ], align, 
						overlaps[j], extendedOverlap ) == 1 )
					{
						/*if ( extendCnt == 0 )
						{
							extendedOverlap = tmpExtendedOverlap ;
							++extendCnt ;
						}
						else if ( tmpExtendedOverlap.similarity == extend) */

						break ;
					}
				}
				if ( j >= overlapCnt )
				{
					delete[] rc ;
					delete[] mrc ;
					++i ;
					continue ;
				}
				
				for ( j = 0 ; j < mateOverlapCnt ; ++j )
				{
					if ( ExtendOverlap( mr, mlen, seqs[ mateOverlaps[j].seqIdx ], align, 
						mateOverlaps[j], mateExtendedOverlap ) == 1 )
						break ;
				}
				if ( j >= mateOverlapCnt )
				{
					delete[] rc ;
					delete[] mrc ;
					++i ;
					continue ;
				}


				int from, to ;
				int fromStart, fromEnd, toStart, toEnd ;
				if ( overlaps[0].strand == 1 )
				{
					from = extendedOverlap.seqIdx ;
					fromStart = extendedOverlap.seqStart ;
					fromEnd = extendedOverlap.seqEnd ;

					to = mateExtendedOverlap.seqIdx ;
					toStart = mateExtendedOverlap.seqStart ;
					toEnd = mateExtendedOverlap.seqEnd ;
				}
				else
				{
					to = extendedOverlap.seqIdx ;
					toStart = extendedOverlap.seqStart ;
					toEnd = extendedOverlap.seqEnd ;

					from = mateExtendedOverlap.seqIdx ;
					fromStart = mateExtendedOverlap.seqStart ;
					fromEnd = mateExtendedOverlap.seqEnd ;
				}
				
				UpdateMateAdjGraph( from, fromStart, fromEnd, to, toStart, toEnd, nextAdj ) ;
				UpdateMateAdjGraph( to, toStart, toEnd, from, fromStart, fromEnd, prevAdj ) ;
				
				printf( "%d(%d %d) %d(%d %d)\n", 
					from, fromStart, fromEnd,
					to, toStart, toEnd ) ;
			
				delete[] rc ;
				delete[] mrc ;
			}
			else
			{
				// Single-end reads, not sure how to use them for now.
				// 	may be look for fragmented secondary alignment?
			}
			
			if ( paired )
				++i ;
		}
		

		// Each Seq will be extend to left and right one step.
		SimpleVector<int> toRemoveSeqIdx ;
		toRemoveSeqIdx.Reserve( seqCnt ) ;
		
		// Find which seq id extend to.
		SimpleVector<struct _pair> matePrevNext ; 
		matePrevNext.ExpandTo( seqCnt ) ;
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
				else if ( prevAdj[i][j].matchCnt == max )
				{
					if ( prevAdj[i][j].seqEnd - prevAdj[i][j].seqStart > 
							prevAdj[i][ prevTag ].seqEnd - prevAdj[i][prevTag].seqStart )
						prevTag = j ;
				}
				printf( "<= %d: %d %d. %d %d %d %d\n", i, prevAdj[i][j].seqIdx, prevAdj[i][j].matchCnt, 
						prevAdj[i][j].readStart, prevAdj[i][j].readEnd,
						prevAdj[i][j].seqStart, prevAdj[i][j].seqEnd ) ;
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
				else if ( nextAdj[i][j].matchCnt == max )
				{
					if ( nextAdj[i][j].seqEnd - nextAdj[i][j].seqStart > 
							nextAdj[i][ nextTag ].seqEnd - nextAdj[i][nextTag].seqStart )
						nextTag = j ;
				}


				printf( "=> %d: %d %d. %d %d %d %d\n", i, nextAdj[i][j].seqIdx, nextAdj[i][j].matchCnt, 
						nextAdj[i][j].readStart, nextAdj[i][j].readEnd,
						nextAdj[i][j].seqStart, nextAdj[i][j].seqEnd ) ;
			}
			
			matePrevNext[i].a = prevTag ;
			matePrevNext[i].b = nextTag ;
		}

		// Figure out the seqs whose extension will create repeat sequence
		// 	e.g.: end-to-end extension from two anchor seq.
		SimpleVector<struct _pair> extensionType ;
		SimpleVector<bool> directlyFilter ;
		extensionType.ExpandTo( seqCnt ) ;
		directlyFilter.ExpandTo( seqCnt ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;

			struct _overlap leftExtend, rightExtend ;
			leftExtend.seqIdx = -1 ;
			rightExtend.seqIdx = -1 ;
			extensionType[i].a = extensionType[i].b = 0 ;
			if ( prevTag >= 0 )
				extensionType[i].a = GetExtendSeqCoord( i, prevAdj[i][ prevTag ], -1, branchAdj, leftExtend ) ;
			if ( nextTag >= 0 )
				extensionType[i].b = GetExtendSeqCoord( i, nextAdj[i][ nextTag ], 1, branchAdj, rightExtend ) ;
		}

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			directlyFilter[i] = false ;
			
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;
			if ( prevTag >= 0 && nextTag >= 0 )
				continue ;
			
			// Keep the second part.
			if ( nextTag >= 0 )
			{
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
			}
		}


		// Do the extesion
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			//if ( inCnt[i] > 0 )
			//	continue ;
			
			// Each seq will try to extend to each direction once.
			if ( directlyFilter[i] )
			{
				toRemoveSeqIdx.PushBack( i ) ;
				continue ;
			}

			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;

			struct _overlap leftExtend, rightExtend ;
			leftExtend.seqIdx = -1 ;
			rightExtend.seqIdx = -1 ;
			if ( prevTag >= 0 )
				GetExtendSeqCoord( i, prevAdj[i][ prevTag ], -1, branchAdj, leftExtend ) ;
			if ( nextTag >= 0 )
				GetExtendSeqCoord( i, nextAdj[i][ nextTag ], 1, branchAdj, rightExtend ) ;
			
			// Compute the sequence for new 
			int newConsensusLen = 0 ; 
			int shift = 0, offset = 0 ;
			
			int start, end ; // we might chop the endings of the current seq.
			start = 0 ;
			end = seqs[i].consensusLen - 1 ;
			if ( leftExtend.seqIdx != -1 )
			{
				shift = leftExtend.seqEnd - leftExtend.seqStart + 1 ;
				start = leftExtend.readStart ;
				newConsensusLen += leftExtend.seqEnd - leftExtend.seqStart + 1 ;
			}
			if ( rightExtend.seqIdx != -1 )
			{
				newConsensusLen += rightExtend.seqEnd - rightExtend.seqStart + 1 ;
				end = rightExtend.readEnd ;
			}
			newConsensusLen += end - start + 1 ;
			if ( newConsensusLen == seqs[i].consensusLen )
				continue ;

			char *newConsensus = (char *)malloc( sizeof( char ) * ( newConsensusLen + 1 ) ) ;
			
			// The ending of current seq might be some random extensions, so
			//    it has lower priority.
			memcpy( newConsensus + shift, seqs[i].consensus + start, end - start + 1 ) ;
			
			if ( leftExtend.seqIdx != -1 )
			{
				memcpy( newConsensus, seqs[ prevAdj[i][prevTag].seqIdx ].consensus + leftExtend.seqStart, 
					leftExtend.seqEnd - leftExtend.seqStart + 1 ) ;
				offset += leftExtend.seqEnd - leftExtend.seqStart + 1 ;
			}
			offset += end - start + 1 ;
			if ( rightExtend.seqIdx != -1 )
			{
				memcpy( newConsensus + offset, seqs[ nextAdj[i][nextTag].seqIdx ].consensus + rightExtend.seqStart, 
					rightExtend.seqEnd - rightExtend.seqStart + 1 ) ;
				offset += rightExtend.seqEnd - rightExtend.seqStart + 1 ;
			}
			newConsensus[offset] = '\0' ;

			struct _seqWrapper ns ;
			ns.isRef = false ;
			ns.consensus = newConsensus ;
			ns.consensusLen = newConsensusLen ;
			
			ns.name = strdup( seqs[i].name ) ;
			ns.posWeight.ExpandTo( newConsensusLen ) ;
			ns.posWeight.SetZero( 0, newConsensusLen ) ;

			for ( j = 0 ; j < newConsensusLen ; ++j )
				if ( newConsensus[j] != 'N' )
					ns.posWeight[j].count[ nucToNum[ newConsensus[j] - 'A' ] ] = 1 ;
			for ( j = 0 ; j < end - start + 1 ; ++j )
				ns.posWeight[j + shift] = seqs[i].posWeight[j + start] ;
		
	
			if ( leftExtend.seqIdx != -1 )
				printf( "prev: %d %s\n", prevAdj[i][prevTag].seqIdx, seqs[ prevAdj[i][prevTag].seqIdx ].consensus ) ;
			printf( "my: %d %s\n", i, seqs[i].consensus ) ;
			if ( rightExtend.seqIdx != -1 )
				printf( "next: %d %s\n", nextAdj[i][nextTag].seqIdx, seqs[ nextAdj[i][nextTag].seqIdx ].consensus ) ;
			printf( "%d new %s\n", i, newConsensus) ;	
			seqs.push_back( ns ) ;
			toRemoveSeqIdx.PushBack( i ) ;
		}
		int size = toRemoveSeqIdx.Size() ;
		for ( i = 0 ; i < size ; ++i )
			ReleaseSeq( toRemoveSeqIdx[i] ) ;	
		Clean( true ) ;

		delete[] branchAdj ;
		delete[] nextAdj ;
		delete[] prevAdj ;
		delete[] align ;
	}

	void Output( FILE *fp )
	{
		int i, j, k ;
		int size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;

			fprintf( fp, ">assemble%d %s\n%s\n", i, seqs[i].name, seqs[i].consensus ) ;
			
			for ( k = 0 ; k < 4 ; ++k )
			{
				for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
					fprintf( fp, "%d ", seqs[i].posWeight[j].count[k] ) ;
				fprintf( fp, "\n" ) ;
			}
		}
	}

} ;


#endif