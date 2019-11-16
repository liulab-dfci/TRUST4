#ifndef _LSONG_KMERINDEX_HEADER
#define _LSONG_KMERINDEX_HEADER

#include <stdio.h>
#include <stdint.h>
#include <map>

#include "defs.h"
#include "KmerCode.hpp"
#include "SimpleVector.hpp"

struct _indexInfo
{
	index_t idx ;
	index_t offset ;
	//int strand ;
} ;


#define KINDEX_HASH_MAX 1000003

class KmerIndex
{
private:
	std::map< uint64_t, SimpleVector<struct _indexInfo> > *index ;
	SimpleVector<struct _indexInfo> nullHit ;

	int GetHash( uint64_t k )
	{
		return k % KINDEX_HASH_MAX ;
	}
public:
	KmerIndex() 
	{
		index = new std::map< uint64_t, SimpleVector<struct _indexInfo> >[KINDEX_HASH_MAX] ;
	}

	~KmerIndex()
	{
		int sum = 0 ;
		//for ( int i = 0 ; i < hashSize ; ++i )
		//{
		//	sum += /*sizeof( hash[i] ) +*/ hash[i].Memory() ;
		//}
		//printf( "%s: %d %d %d\n", __func__, hashSize, sizeof( hash[0] ), hash[0].Memory() ) ;
		if ( index != NULL )
			delete[] index ;
	}

	void Clear()
	{
		//index.clear() ;
		int i ;
		for ( i = 0 ; i < KINDEX_HASH_MAX ; ++i )	
			std::map<uint64_t, SimpleVector<struct _indexInfo> >().swap( index[i] ) ;
	}

	void Insert( KmerCode &kmerCode, index_t idx, index_t offset, int strand )
	{
		if ( !kmerCode.IsValid() )
			return ;
		struct _indexInfo newEntry ;
		newEntry.idx = idx ;
		newEntry.offset = offset ;
		//newEntry.strand = strand ;
		uint64_t kcode = kmerCode.GetCode() ;
		int h = GetHash( kcode ) ;
		index[h][ kcode ].PushBack( newEntry ) ;
		
		//printf( "%d\n", hash[key].Memory() ) ;
	}

	void Remove( KmerCode &kmerCode, index_t idx, index_t offset, int strand )
	{
		if ( !kmerCode.IsValid() )
			return ;

		SimpleVector<struct _indexInfo> &list = *Search( kmerCode ) ;
		int size = list.Size() ;
		int i ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( list[i].idx == idx && list[i].offset == offset )
				break ;
		}
		if ( i < size )
		{
			//printf( "remove %d %d\n", idx, offset ) ;
			list.Remove( i ) ;
		}
	}

	SimpleVector<struct _indexInfo> *Search( KmerCode &kmerCode )
	{
		if ( !kmerCode.IsValid() )
			return &nullHit ;
		uint64_t kcode = kmerCode.GetCode() ;
		int h = GetHash( kcode ) ;

		std::map< uint64_t, SimpleVector<struct _indexInfo> >::iterator it = index[h].find( kcode ) ;
		if ( it == index[h].end() )
			return &nullHit ;
		else
			return &(it->second) ;
	}

	void BuildIndexFromRead( KmerCode &kmerCode, char *s, int len, int id, int shift = 0 )
	{
		int i ;
		int kl = kmerCode.GetKmerLength() ;
		if ( len < kl )
			return ;
		kmerCode.Restart() ;
		//KmerCode rcKmerCode( kl ) ;
		KmerCode prevKmerCode( kl ) ;
		for ( i = 0 ; i < kl - 1 ; ++i )
			kmerCode.Append( s[i] ) ;
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( s[i] ) ;
			if ( kmerCode.IsValid() && ( i == kl || !kmerCode.IsEqual( prevKmerCode ) ) )
			{
				Insert( kmerCode, id, i - kl + 1 + shift, 1 ) ;

				//rcKmerCode.SetCode( kmerCode.GetReverseComplementCode() ) ;
				//Insert( rcKmerCode, id, len - 1 - i, -1 ) ;
			}
			prevKmerCode = kmerCode ;
		}
	}
	
	// When merging or extending sequences, there kmer position will shift and change id.
	void UpdateIndexFromRead( KmerCode &kmerCode, char *s, int len, int shift, int oldId, int id ) 
	{
		int i, j ;
		int kl = kmerCode.GetKmerLength() ;
		if ( len < kl )
			return ;
		kmerCode.Restart() ;

		for ( i = 0 ; i < kl - 1 ; ++i )
			kmerCode.Append( s[i] ) ;
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( s[i] ) ;
			if ( kmerCode.IsValid() )
			{
				SimpleVector<struct _indexInfo> &list = *Search( kmerCode ) ;	
				int size = list.Size() ;
				//printf( "kcode: %lld\n", kmerCode.GetCode() ) ;
				//for ( j = 0 ; j < size ; ++j )
				//	printf( "before %d %d\n", list[j].idx, list[j].offset ) ;
				for ( j = 0 ; j < size ; ++j )
				{
					if ( list[j].idx == oldId && list[j].offset == i - kl + 1 )
					{
						//printf( "update %d->%d: %d %d\n", oldId, id, i - kl + 1, shift ) ;
						list[j].idx = id ;
						list[j].offset += shift ;
						break ;
					}
				}
			}
			SimpleVector<struct _indexInfo> &list = *Search( kmerCode ) ;	
			int size = list.Size() ;

			//for ( j = 0 ; j < size ; ++j )
			//	printf( "test %d %d\n", list[j].idx, list[j].offset ) ;
		}

	}

	void RemoveIndexFromRead( KmerCode &kmerCode, char *s, int len, int id, int offset )
	{
		int i, j ;
		int kl = kmerCode.GetKmerLength() ;
		if ( len < kl )
			return ;
		kmerCode.Restart() ;

		for ( i = 0 ; i < kl - 1 ; ++i )
			kmerCode.Append( s[i] ) ;
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( s[i] ) ;
			if ( kmerCode.IsValid() )
			{
				Remove( kmerCode, id, i - kl + 1 + offset, 1 ) ;
			}
		}
	}
} ;

#endif
