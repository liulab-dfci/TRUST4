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
	char *consensus ;
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
					//printf( "%d: %d %d %d %d\n", i, nh.readOffset, nh.indexHit.idx, nh.indexHit.offset, nh.indexHit.strand ) ;
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
				if ( no.strand == 1 )
				{
					no.seqStart = finalHits[0].indexHit.offset ;
					no.seqEnd = finalHits[ lisSize - 1 ].indexHit.offset + kmerLength - 1 ;
				}
				else
				{
					// if on other strand, the index starting the kmer offset from the other end.
					no.seqStart = finalHits[0].indexHit.offset - kmerLength + 1 ;
					no.seqEnd = finalHits[ lisSize - 1].indexHit.offset ;
				}

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
		char *rcRead = (char *)malloc( sizeof( char ) * ( len + 1 ) ) ;
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
		free( rcRead ) ;

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

		rcRead = (char *)malloc( sizeof( char ) * ( len + 1 ) ) ;
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
					
						for ( k = 0 ; align[k] != -1 ; ++k )
						{
							if ( align[k] == EDIT_MATCH )
								++matchCnt ;
							else if ( align[k] == EDIT_MISMATCH )
								++mismatchCnt ;
							else 
								++indelCnt ;
						}

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
					
						for ( k = 0 ; align[k] != -1 ; ++k )
						{
							if ( align[k] == EDIT_MATCH )
								++matchCnt ;
							else if ( align[k] == EDIT_MISMATCH )
								++mismatchCnt ;
							else 
								++indelCnt ;
						}

						if ( !seqs[ overlaps[i].seqIdx ].isRef && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}

					}
				}
			} // for j
			free( align ) ;

			overlaps[i].similarity = (double)matchCnt / ( overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 
								overlaps[i].readEnd - overlaps[i].readStart + 1 ) ;
		} // for i
		free( rcRead ) ;

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


	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set.
	int AddRead( char *read )
	{
		//printf( "%s\n", seq ) ;
		int i, k ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;

		overlapCnt = GetOverlapsFromRead( read, overlaps ) ;
		
		for ( i = 0 ; i < overlapCnt ; ++i )
			printf( "%d: %d %s. %d %d %d %d\n", i, overlaps[i].seqIdx, seqs[ overlaps[i].seqIdx ].name, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ;
		
		
		std::sort( overlaps.begin(), overlaps.end() ) ;
		
		// If the read only overlaps with the reference, we will add that to the seq.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( !seqs[ overlaps[i].seqIdx ].isRef )
				break ;
		}
		
		struct _pair *chosenOverlapOnRead = new struct _pair[ overlapCnt ];
		k = 0 ;
		if ( i >= overlapCnt )
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

			int idx = seqs.size() ;
			struct _seqWrapper ns ;

			ns.name = NULL ;
			ns.consensuss = strdup( read ) ;
			ns.isRef = false ;
			
			int len = strlen( read ) ;
			ns.posWeight.Reserve( len ) ;
			for ( i = 0 ; i < len ; ++i )
			{
				memset( ns.posWeight[i].count, 0, sizeof( ns.posWeight[i].count ) ) ;
				++ns.posWeight[ nucToNum[ read[i] - 'A' ] ] ;
			}
			seqs.push_back( ns ) ;
		}
		else
		{
			// Incorporate to added sequences.
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				for ( j = 0 ; j < k ; ++j )
				{
					if ( overlaps[])
				}
			}
		}

		delete[] chosenOverlapOnRead ;

		return 0 ;
	}
} ;


#endif
