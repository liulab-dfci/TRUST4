// The call handles kmer count

#ifndef _LSONG_KMERCOUNT_HEADER
#define _LSONG_KMERCOUNT_HEADER

#include <map>
#include <algorithm>

#include "KmerCode.hpp"

#define KCOUNT_HASH_MAX 1000003

class KmerCount
{
private:
	std::map<uint64_t, int> *count ;
	int kmerLength ;
	KmerCode kmerCode ;
	int maxReadLen ;

	int *c ;

	int GetHash( uint64_t k )
	{
		return k % KCOUNT_HASH_MAX ;
	}
public:
	KmerCount( int k ): kmerCode( k ) 
	{ 
		kmerLength = k ; 
		maxReadLen = -1 ;
		c = NULL ;
		count = new std::map<uint64_t, int>[KCOUNT_HASH_MAX] ;
	}
	~KmerCount() 
	{
		if ( c != NULL )
			delete[] c ;
		if ( count != NULL )
			delete[] count ;
	}

	int AddCount( char *read )
	{
		int i ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return 0 ;

		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;

		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			if ( kmerCode.IsValid() )
			{
				uint64_t kcode = kmerCode.GetCanonicalKmerCode() ;
				++count[ GetHash(kcode) ][ kcode ] ;
				/*if ( count[ GetHash( kcode ) ][ kcode ] >= 500 )
				{
					printf( "%s\n", read + i - kmerLength + 1 ) ;
				}*/
			}
		}
		
		if ( len > maxReadLen )
			maxReadLen = len ;
		return 1 ;
	}

	void AddCountFromFile( char *file )
	{
		FILE *fp = fopen( file, "r" ) ;
		char buffer[100] ;
		int i ;

		while ( fscanf( fp, "%s", buffer ) != EOF )
		{
			int c = atoi( &buffer[ 1 ] ) ;
			fscanf( fp, "%s", buffer ) ;
			if ( c <= 1 )
				continue ;

			kmerCode.Restart() ;
			for ( i = 0 ; buffer[i] ; ++i )
				kmerCode.Append( buffer[i] ) ;
			uint64_t kcode = kmerCode.GetCode() ;
			count[ GetHash( kcode ) ][ kcode ] = c ;
		}
		fclose( fp ) ;
	}

	void SetBuffer( int sz )
	{
		maxReadLen = sz ;
		if ( c == NULL )
			 c = new int[ sz ] ;
	}

	int GetCountStats( char *read, int &minCount, int &medianCount, double &avgCount )
	{
		int i, k ;
		int sum ;
		if ( maxReadLen == -1 )
			return 0 ;

		if ( c == NULL )
			c = new int[ maxReadLen ] ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
		{
			minCount = medianCount = avgCount = -1 ;
			return 0 ;
		}

		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;
		k = 0 ;
		sum = 0 ;
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			if ( kmerCode.IsValid() )
			{
				uint64_t kcode = kmerCode.GetCanonicalKmerCode() ;
				c[k] = count[ GetHash( kcode ) ][ kcode ] ;
				if ( c[k] <= 0 )
					c[k] = 1 ;
				sum += c[k] ;
				++k ;
			}
		}
		if ( k == 0 )
			return 0 ;

		std::sort( c, c + k ) ; 
		minCount = c[0] ;
		medianCount = c[k / 2] ;
		avgCount = sum / (double)k ;
		for ( i = 0 ; i < len ; ++i )
			if ( read[i] == 'N' )
			{
				if ( minCount >= 0 )
					minCount = 0 ;
				else if ( minCount <= 0 )
					--minCount ;
			}

		return 1 ;
	}

	void Release()
	{
		if ( c != NULL )
			delete[] c ;
		if ( count != NULL )
			delete[] count ;

		c = NULL ;
		count = NULL ;
	}
} ;

#endif
