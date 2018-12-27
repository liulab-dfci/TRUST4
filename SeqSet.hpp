// The data structure holds the set of sequences (can be "assembled" from several reads)
#ifndef _MOURISL_SEQSET_HEADER
#define _MOURISL_SEQSET_HEADER

#include <string.h>
#include <algorithm>
#include <vector>

#include "SimpleVector.hpp"
#include "KmerIndex.hpp"
#include "ReadFiles.hpp"
#include "AlignAlgo.hpp"

struct _posWeight
{
	int count[4] ;
} ;

struct _seqWrapper
{
	char *name ;
	char *consensus ; // This should be handled by malloc/free.
	int consensusLen ;
	SimpleVector<struct _posWeight> posWeight ;
	bool isRef ; // this is from reference.
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
	
	double similarity ;

	SimpleVector<struct _pair> hitCoords ;

	bool operator<( const struct _overlap &b ) const
	{
		if ( similarity != b.similarity )
			return similarity > b.similarity ; // The overlap with higher similarity should come first.
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

	static bool CompSortPairBInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.b < p2.b ;
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
			rcSeq[i] = numToNuc[ 3 - nucToNum[seq[len - 1 - i] - 'A'] ];
		}
		rcSeq[i] = '\0' ;
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

public:
	SeqSet( int kl ) 
	{
		kmerLength = kl ;
		minHitRequired = 3 ;
		radius = 10 ;
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
	
	// Input some baseline sequence to match against.
	void InputRefFa( char *filename ) 
	{
		ReadFiles fa ;
		fa.AddReadFile( filename ) ;
		
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

	// Use the hits to extract overlaps from SeqSet
	int GetOverlapsFromHits( SimpleVector<struct _hit> &hits, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int hitSize = hits.Size() ;
		
		SimpleVector<struct _pair> hitCoordDiff ;
		hitCoordDiff.Reserve( hitSize ) ;
		SimpleVector<struct _pair> concordantHitCoord ;
		SimpleVector<struct _pair> hitCoordLIS ;
		SimpleVector<struct _hit> finalHits ;

		int hitLenRequired = 31  ;

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
				//	printf( "%d: %d %d %d\n", i, hits[i].indexHit.idx, concordantHitCoord[k].a, concordantHitCoord[k].b ) ;


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

		// Locate the hits from the same-strand case.
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;

		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
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
		// Locate the hits from the opposite-strand case.
		char *rcRead =  new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( rcRead[i] ) ;
		
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( rcRead[i] ) ;
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
			}
		}
		delete[] rcRead ;

		// Find the overlaps.
		std::sort( hits.BeginAddress(), hits.EndAddress() ) ;
		//for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
		//	printf( "- %d %d %d %d\n", it->readOffset, it->indexHit.idx, it->indexHit.offset, it->strand ) ;

		int overlapCnt = GetOverlapsFromHits( hits,  overlaps ) ;

		//for ( i = 0 ; i < overlapCnt ; ++i )
		//	printf( "%d: %d %s. %d %d %d %d\n", i, overlaps[i].poaIdx, poas[ overlaps[i].poaIdx ].name, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].poaStart, overlaps[i].poaEnd ) ;
		// Determine whether we want to add this reads by looking at the quality of overlap
		if ( overlapCnt == 0 )
			return 0 ;

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
		
		// Compute wehterh the extension matched right.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			;
		}

		// Compute similarity overlaps
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *r ;
			if ( overlaps[i].strand == 1 )
				r = read ;
			else
				r = rcRead ;

			SimpleVector<struct _pair> &hitCoords = overlaps[i].hitCoords ; 	
			
			int hitCnt = hitCoords.Size() ;
			int matchCnt = 0, mismatchCnt = 0, indelCnt = 0  ;
			double similarity = 0 ;
			
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
								
						AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - hitCoords[j - 1].b - 1,
								r + hitCoords[j - 1].a + kmerLength, hitCoords[j].a - hitCoords[j - 1].a - 1, 
								align ) ;

						GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;		

						if ( !seqs[ overlaps[i].seqIdx ].isRef && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}
					}
				}
				else
				{
					if ( !seqs[ overlaps[i].seqIdx ].isRef )
					{
						similarity = 0 ;
						break ;
					}

					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 < hitCoords[j].b )
					{
						matchCnt += 2 * kmerLength ;
						// Make the two kmer hit match on coordinate.
						indelCnt += ( hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) + 
							( hitCoords[j].a + kmerLength - hitCoords[j - 1].a )  ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 < hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						matchCnt += 2 * kmerLength ;
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
						
						AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - hitCoords[j - 1].b - 1,
								r + hitCoords[j - 1].a + kmerLength, hitCoords[j].a - hitCoords[j - 1].a - 1, 
								align ) ;	
						
						GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;

						if ( !seqs[ overlaps[i].seqIdx ].isRef && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}

					}
				}
			} // for j
			delete[] align ;

			overlaps[i].similarity = (double)matchCnt / ( overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 
								overlaps[i].readEnd - overlaps[i].readStart + 1 ) ;
		} // for i
		delete[] rcRead ;

		k = 0 ; 
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < 0.75 )
				continue ;
			else if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < 0.9 )
				continue ;

			//printf( "%lf\n", overlaps[i].similarity ) ;
			overlaps[k] = overlaps[i] ;
			++k ;
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		

		return overlapCnt ;
	}
	
	// Extend the overlap to include the overhang parts and filter the overlaps if the overhang does not match well.
	// return: whether this is a valid extension or not
	int ExtendOverlap( char *r, int len, struct _seqWrapper &seq, char *align, struct _overlap &overlap, struct _overlap &extendedOverlap )
	{
		// Check whether the overhang part is compatible with each other or not.
		// Extension to 5'-end ( left end )
		int matchCnt, mismatchCnt, indelCnt ;
		int leftOverhangSize = MIN( overlap.readStart, overlap.seqStart ) ;

		AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqStart - leftOverhangSize, leftOverhangSize, 
				r + overlap.readStart - leftOverhangSize, leftOverhangSize, align ) ;

		GetAlignStats( align, false, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > 0 )
			return 0 ;

		// Extension to 3'-end ( right end )
		int rightOverhangSize = MIN( len - 1 - overlap.readEnd, seq.consensusLen - 1 - overlap.seqEnd ) ;

		AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqEnd + 1, rightOverhangSize,
				r + overlap.readEnd + 1, rightOverhangSize, align ) ;
		GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > 0 )
			return 0 ;

		if ( (double)mismatchCnt / ( leftOverhangSize + rightOverhangSize ) > 1.5 / kmerLength ) 
			return 0 ;

		extendedOverlap.seqIdx = overlap.seqIdx ;
		extendedOverlap.readStart = overlap.readStart - leftOverhangSize ;
		extendedOverlap.readEnd = overlap.readEnd + rightOverhangSize ;
		extendedOverlap.seqStart = overlap.seqStart - leftOverhangSize ;
		extendedOverlap.seqEnd = overlap.seqEnd + rightOverhangSize ;
		
		return 1 ;
	}

	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set.
	int AddRead( char *read )
	{
		//printf( "%s\n", seq ) ;
		int i, j, k ;
		int len = strlen( read ) ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;

		overlapCnt = GetOverlapsFromRead( read, overlaps ) ;
		
		if ( overlapCnt == 0 )
			return -1 ;

		for ( i = 0 ; i < overlapCnt ; ++i )
			printf( "%d: %d %s. %d. %d %d %d %d\n", i, overlaps[i].seqIdx, seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, 
					overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ; 
		fflush( stdout ) ;	
		std::sort( overlaps.begin(), overlaps.end() ) ;
		
		// If the read only overlaps with the reference, we will add that to the seq.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( !seqs[ overlaps[i].seqIdx ].isRef )
				break ;
		}
		
		struct _overlap *extendedOverlaps = new struct _overlap[ overlapCnt ];
		k = 0 ;
		int ret = 0 ;
		bool addNew = true ;
		if ( i < overlapCnt )
		{
			// Incorporate to existing sequence.
			char *align = new char[3 * len] ;
			char *rcRead = strdup( read ) ;
			ReverseComplement( rcRead, read, len ) ;
			
			char *r ;
			int readInConsensusOffset = 0 ;
			int seqIdx ;

			if ( overlaps[0].strand == 1 )
				r = read ;
			else
				r = rcRead ;

			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				for ( j = 0 ; j < k ; ++j )
				{
					if ( overlaps[i].readStart >= extendedOverlaps[j].readStart - radius  
						&& overlaps[i].readEnd <= extendedOverlaps[j].readEnd + radius )
						break ;
				}
				
				struct _seqWrapper &seq = seqs[ overlaps[i].seqIdx ] ; 
				if ( j < k || seq.isRef )
					continue ;

				// Only extend the novel seqs.
				if ( ExtendOverlap( r, len, seq, align, overlaps[i], extendedOverlaps[k] ) == 1 )
					++k ;
			}	

			if ( k > 1 )
			{
				int eOverlapCnt = k ;
				addNew = false ;		
				// Merge sequences.
				// Reorder the overlaps to the order on the read coordinate.
				std::sort( extendedOverlaps, extendedOverlaps + ( k - 1 ), CompSortOverlapsOnReadCoord ) ;
				
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
				/*for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					printf( "%d: %d %d %d %d. %d\n", i, extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, 
						extendedOverlaps[i].seqStart, extendedOverlaps[i].seqEnd,
						seqOffset[i] ) ;
				}*/
				// Copy the original consensus in.
				for ( i = 0 ; i < eOverlapCnt ; ++i )
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
				SimpleVector<struct _posWeight> &posWeight = seqs[newSeqIdx].posWeight ;
				posWeight.ExpandTo( newConsensusLen ) ;
				posWeight.SetZero( 0, newConsensusLen ) ;
				

				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					// Though not the most efficient implementation, it seems very straightforward.
					int seqIdx = extendedOverlaps[i].seqIdx ;
					for ( j = 0 ; j < seqs[ seqIdx ].consensusLen ; ++j )
					{
						int l ;
						for ( l = 0 ; l < 4 ; ++l )
							posWeight[ seqOffset[i] + j ].count[l] += seqs[ seqIdx ].posWeight[j].count[l] ;
					}
				}

				// Update the index.
				KmerCode kmerCode( kmerLength ) ;
				if ( seqOffset[0] != 0 )
				{
					seqIndex.UpdateIndex( kmerCode, seqs[ extendedOverlaps[0].seqIdx ].consensus, 
							seqs[ extendedOverlaps[0].seqIdx].consensusLen, seqOffset[0], 
							extendedOverlaps[0].seqIdx, newSeqIdx ) ; 
				}
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					seqIndex.UpdateIndex( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus, 
							seqs[ extendedOverlaps[i].seqIdx].consensusLen, seqOffset[i], 
							extendedOverlaps[i].seqIdx, newSeqIdx ) ;
				}
				
				// Update the index for the gap.
				if ( extendedOverlaps[0].readStart > 0 )
				{
					seqIndex.BuildIndexFromRead( kmerCode, r, extendedOverlaps[0].readStart + kmerLength - 1, newSeqIdx ) ;
				}
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( extendedOverlaps[i].readStart > extendedOverlaps[i - 1].readEnd + 1 )
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
				// Update the name.
				sum = 0 ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					sum += strlen( seqs[ extendedOverlaps[i].seqIdx ].name ) ;
				char* nameBuffer = new char[sum] ;
				
				strcpy( nameBuffer, seqs[ newSeqIdx ].name ) ;
				sum = strlen( nameBuffer ) ;
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( strcmp( seqs[ extendedOverlaps[i].seqIdx ].name, seqs[ extendedOverlaps[i - 1].seqIdx ].name ) )
					{
						nameBuffer[ sum ] = '+' ;
						strcat( nameBuffer + sum + 1, seqs[ extendedOverlaps[i - 1].seqIdx ].name ) ;
						sum = sum + 1 + strlen( seqs[ extendedOverlaps[i - 1].seqIdx ].name ) ;
					}
				}
				free( seqs[ newSeqIdx ].name ) ;
				seqs[ newSeqIdx ].name = strdup( nameBuffer ) ;
				delete[] nameBuffer ;

				// Relase the memory for merged seqs.
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					int seqIdx = extendedOverlaps[i].seqIdx ;
					free( seqs[ seqIdx ].name ) ;
					free( seqs[ seqIdx ].consensus ) ;
					seqs[ seqIdx ].name = NULL ;
					seqs[ seqIdx ].consensus = NULL ;
					seqs[ seqIdx ].posWeight.Release() ;
				}
					
				// Put the new allocated stuff in.
				free( seqs[ newSeqIdx ].consensus ) ;
				seqs[ newSeqIdx ].consensus = newConsensus ;
				seqs[ newSeqIdx ].consensusLen = newConsensusLen ;

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
					//printf( "new consensus %s %d. %s %d\n", seq.consensus, seq.consensusLen, newConsensus, j ) ;
					
					// Update index 
					int shift = extendedOverlaps[0].readStart ;
					KmerCode kmerCode( kmerLength ) ;
					if ( shift > 0 )
					{
						seqIndex.BuildIndexFromRead( kmerCode, r, extendedOverlaps[0].readStart + kmerLength - 1, seqIdx ) ;
						seqIndex.UpdateIndex( kmerCode, seq.consensus, seq.consensusLen, shift, seqIdx, seqIdx ) ; 
					}
					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						int rstart = extendedOverlaps[0].readEnd - kmerLength + 2 ;
						seqIndex.BuildIndexFromRead( kmerCode, r + rstart , 
							( len - rstart ), seqIdx, 
							extendedOverlaps[0].readStart + seq.consensusLen - kmerLength + 1 ) ;
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

					// either one of the ends of read or seq should be 0.
					readInConsensusOffset = 0 ;
					if ( extendedOverlaps[0].seqStart > 0 )
						readInConsensusOffset = extendedOverlaps[0].seqStart ;

					free( seq.consensus ) ;
					seq.consensus = newConsensus ;
					seq.consensusLen = strlen( newConsensus ) ;	
					//printf( "new consensus len %d\n", seq.consensusLen ) ;
				}
				else
					readInConsensusOffset = extendedOverlaps[0].seqStart ;
			}

			// Update the posweight, assume we already compute the new readStart and shift existing posWeight.
			// seqIdx holds the index that we need to update.
			if ( !addNew )
			{
				struct _seqWrapper &seq = seqs[seqIdx] ;
				for ( i = 0 ; i < len ; ++i )
					++seq.posWeight[i + readInConsensusOffset].count[ nucToNum[ r[i] - 'A' ] ] ;
			}
			
			free( rcRead ) ;
			delete[] align ;
		}

		k = 0 ;
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

			ns.name = strdup( seqs[ overlaps[0].seqIdx ].name ) ;
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
				++ns.posWeight[i].count[ nucToNum[ ns.consensus[i] - 'A' ] ] ;
			}
			seqs.push_back( ns ) ;

			// Don't forget to update index.
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, len, idx ) ;			
		
			ret = idx ;
		}

		delete[] extendedOverlaps ;

		return ret ;
	}

	void Output()
	{
		int i ;
		int size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;

			printf( ">%s\n%s\n", seqs[i].name, seqs[i].consensus ) ;
		}
	}
} ;


#endif
