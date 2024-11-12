// The call handles kmer count

#ifndef _LSONG_KMERCOUNT_HEADER
#define _LSONG_KMERCOUNT_HEADER

#include <map>
#include <algorithm>
#include <pthread.h>

#include "KmerCode.hpp"

class KmerCount
{
private:
	std::map<uint64_t, int> *count ;
	int kmerLength ;
	int maxReadLen ;
	int khashMax ;
	pthread_mutex_t **locks ;

	int *c ;

	int GetHash( uint64_t k )
	{
		return k % khashMax ;
	}
public:
	KmerCount( int k, int hmax = 1000003 )
	{ 
		kmerLength = k ;
		khashMax = hmax ;
		maxReadLen = -1 ;
		c = NULL ;
		count = new std::map<uint64_t, int>[khashMax] ;
		locks = NULL ;
	}

	KmerCount(const KmerCount &b)
	{
		kmerLength = b.kmerLength ;
		khashMax = b.khashMax ;
		maxReadLen = -1 ;
		c = NULL ;
		count = new std::map<uint64_t, int>[khashMax] ;
		locks = NULL ;
	}

	~KmerCount() 
	{
		Release() ;
	}

	void SetPthreadLocks() 
	{
		locks = new pthread_mutex_t *[khashMax] ;
		int i ;
		for (i = 0 ; i < khashMax ; ++i)
		{
			locks[i] = new pthread_mutex_t ;
			pthread_mutex_init(locks[i], NULL) ;
		}
	}

	int AddCount( char *read )
	{
		int i ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return 0 ;

		KmerCode kmerCode(kmerLength) ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;

		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			if ( kmerCode.IsValid() )
			{
				uint64_t kcode = kmerCode.GetCanonicalKmerCode() ;
				int h = GetHash(kcode) ;
				if (locks)
					pthread_mutex_lock(locks[h]) ;
				++count[h][ kcode ] ;
				if (locks)
					pthread_mutex_unlock(locks[h]) ;
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

		KmerCode kmerCode(kmerLength) ;
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

	void Output(FILE *fp)
	{
		int i, j ;
		char *buffer = new char[kmerLength + 1] ;
		for ( i = 0 ; i < khashMax ; ++i )	
		{
			for ( std::map<uint64_t, int>::iterator it = count[i].begin() ; it != count[i].end() ; ++it )
			{
				if ( it->second <= 1 )
					continue ;
				
				for ( j = 0 ; j < kmerLength ; ++j )
				{
					buffer[j] = numToNuc[ ( it->first >> ( 2 * j ) ) & 3 ] ;  
				}
				buffer[j] = '\0' ;
				fprintf( fp, ">%d\n%s\n", it->second, buffer ) ;
			}
		}
		delete[] buffer ;
	}

	void SetBuffer( int sz )
	{
		maxReadLen = sz ;
		if ( c == NULL )
			 c = new int[ sz ] ;
	}

	void SetBuffer()
	{
		if (c == NULL && maxReadLen > 0)
			c = new int[maxReadLen] ;
	}
	
	int GetCount( char *kmer )
	{
		int i ;
		KmerCode kmerCode(kmerLength) ;
		for ( i = 0 ; i < kmerLength ; ++i )
			kmerCode.Append( kmer[i] ) ;
		if ( kmerCode.IsValid() )
		{
			uint64_t kcode = kmerCode.GetCanonicalKmerCode() ;
			int key = GetHash( kcode ) ;
			if ( count[ key ].find( kcode ) == count[key].end() )
				return 0 ;
			else
				return count[ GetHash( kcode ) ][ kcode ] ;
		}
		else
			return 0 ;
	}

	// This is thread-safe when providedBuffer is given
	int GetCountStatsAndTrim( char *read, char *qual, int &minCount, int &medianCount, float &avgCount,
			int *providedBuffer = NULL)
	{
		int i, k ;
		int sum ;
		//minCount = medianCount = avgCount = 1 ;
		//return 0 ;
		if ( maxReadLen == -1 )
			return 0 ;

		int *c = providedBuffer ;
		if ( c == NULL )
		{
			if (this->c == NULL)
				this->c = new int[ maxReadLen ] ;
			c = this->c ;
		}
		int len = strlen( read ) ;
		if ( len < kmerLength )
		{
			minCount = medianCount = avgCount = -1 ;
			return 0 ;
		}

		KmerCode kmerCode(kmerLength) ; 
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
				int key = GetHash( kcode ) ;
				if ( count[ key ].find( kcode ) == count[key].end() )
					c[k] = 0 ;
				else
					c[k] = count[ GetHash( kcode ) ][ kcode ] ;
				if ( c[k] <= 0 )
					c[k] = 1 ;
				sum += c[k] ;
				++k ;
			}
		}
		if ( k == 0 )
		{
			minCount = -len ;
			medianCount = -len ;
			avgCount = -len ;
			// Don't trim the read (kind of) when quality score is not give
			//   this help keep the function consistent that it only trims when
			//   qual score is given.
			if (qual != NULL) 
				read[0] = '\0' ;
			return 0 ;
		}
		
		// Do the trimming.
		if ( qual != NULL )
		{
			int j ;
			for ( i = k - 1 ; i >= 0 ; --i )
			{
				if ( c[i] > 1 )
					break ;
			}
			
			// Locate the low quality region.
			++i ;
			//for ( j = i + kmerLength - 1 ; j < len ; ++j )
			int badCnt = 0 ;
			int trimStart = -1 ;
			for ( j = len - 1 ; j >= i + kmerLength - 1 ; --j )
			{
				if ( qual[j] - 32 <= 15 )	
				{
					++badCnt ;
					if ( badCnt >= 0.1 * ( len - j ) )
						trimStart = j ;
				}
			}
			if ( trimStart > 0 )
			{
				k = trimStart - kmerLength + 1 ;
				read[ trimStart ] = '\0' ;
				qual[ trimStart ] = '\0' ;
			}
			if ( trimStart > 0 && trimStart < kmerLength )
			{
				k = 0 ;
				read[0] = '\0' ;
				qual[0] = '\0' ;
			}
		}

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

  void Clear()
  {
    int i ;
    for (i = 0 ; i < khashMax ; ++i)
      count[i].clear() ;
  }

	void Release()
	{
		if ( c != NULL )
			delete[] c ;
		if ( count != NULL )
			delete[] count ;

		c = NULL ;
		count = NULL ;

		if (locks != NULL)
		{
			int i ;
			for (i = 0 ; i < khashMax ; ++i)
			{
				pthread_mutex_destroy(locks[i]) ;
				delete locks[i] ;
			}
			delete[] locks ;
			locks = NULL ;
		}
	}
} ;

#endif
