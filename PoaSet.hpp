#ifndef _MOURISL_POASET_HEADER
#define _MOURISL_POASET_HEADER

#include <string.h>
#include <vector>
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
		else if ( indexHit.offset != b.indexHit.offset )
			return indexHit.offset < b.indexHit.offset ;
		else if ( readOffset != b.readOffset )
			return  readOffset < b.readOffset ;
		
		return false ;
	}
} ;

class POASet
{
private:
	std::vector<struct _poaWrapper> poas ;
	KmerIndex poaIndex ;
	int kmerLength ;

public:
	POASet( int kl ) 
	{
		kmerLength = kl ;
	}
	~POASet() 
	{
		int size ;
		int i ;
		size = poas.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			//poas[i].poa.Release() ;
			free( poas[i].name ) ;	
		}
	}

	int Size()
	{
		return poas.size() ;
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

			int id = poas.size() ;
			poas.push_back( np ) ;

			struct _poaWrapper &pw = poas[id] ;
			int seqLen = strlen( fa.seq ) ;
			poas[id].poa.Initialize( fa.seq, seqLen, true  ) ;
			poaIndex.BuildIndexFromRead( kmerCode, fa.seq, seqLen, id ) ;
		}
	}

	// Compute the length of hit from the read, take the overlaps of kmer into account 
	int GetTotalHitLength( std::vector<struct _hit> &hits )
	{
		int hitSize = hits.size() ;
		int i, j ;
		int ret = 0 ;
		for ( i = 0 ; i < hitSize ; ++i )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
			{
				if ( hits[j].readOffset > hits[j - 1].readOffset + kmerLength - 1 )	
					break ;
			}

			ret += hits[j].readOffset - hits[i].readOffset + kmerLength ;

			i = j ;
		}
		return ret ;
	}

	// Use the hits to extract overlaps from 
	int GetPoaOverlaps( std::vector<struct _hit> hits, std::vector<int> &poaIdx, std::vector<struct _pair> &poaRange, 
				std::vector<struct _pair> &readRange )
	{
		int i, j ;
		int hitSize = hits.size() ;
		for ( i = 0 ; i < hitSize ; ++i )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].strand != hits[i].strand || hits[j].indexHit)
		
			i = j ;
		}
	}

	
	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set.
	int AddRead( char *seq )
	{
		int i, j, k ;
		int len = strlen( seq ) ;
		if ( len < kmerLength )
			return -1 ;

		std::vector<struct _hit> hits ;		    	

		KmerCode kmerCode( kmerLength ) ;

		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( seq[i] ) ;

		for ( ; i < len ; ++i )
		{
			kmerCode.Append( seq[i] ) ;
			SimpleVector<struct _indexInfo> &indexHit = *poaIndex.Search( kmerCode ) ; 
			
			int size = indexHit.Size() ;

			if ( size >= 40 )
				continue ;
			
			for ( j = 0 ; j < size ; ++j )
			{
				struct _hit nh ;
				nh.indexHit = indexHit[j] ;
				nh.readOffset = i - kmerLength + 1 ;
				nh.readStrand = 1 ;
				hits.push_back( nh ) ;
			}
		}
		
		std::sort( hits.begin(), hits.end() ) ;
		//for ( std::vector<struct _hit>::iterator it=hits.begin() ; it != hits.end() ; ++it )
		//	printf( "- %d\n", it->indexHit.strand ) ;
		
		std::vector<int> poaIdx ;
		std::vector<struct _pair> rangePoa, rangeRead ; 
		GetPoaOverlaps( hits,  poaIdx, rangePoa, rangeRead ) ;

		return 0 ;
	}
} ;


#endif
