#ifndef _LSONG_SIMPLE_VECTOR_HEADER
#define _LSONG_SIMPLE_VECTOR_HEADER

// A light version of vector, which increase the size of the array by 
// a value no more than specified if it got overflow.
// And the type of elements is basic.

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//const int maxInc = -1 ;

template <class T>
class SimpleVector
{
private:
	int size ;
	int capacity ;
	int maxInc ; // The maximal value we can use to increase the capacity.
	int inc ;
	T *s ;
public:
	SimpleVector() : maxInc( -1 )
	{ 
		s = NULL ;
		size = capacity = 0 ;
		inc = 1 ;
	}
	
	SimpleVector( int mi ): maxInc( mi ) 
	{ 
		s = NULL ;
		size = capacity = 0 ;
		inc = 1 ;
	}

	SimpleVector( const SimpleVector &in )
	{
		size = in.size ;
		capacity = in.capacity ;
		if ( capacity > 0 )
		{
			//s = in.s ;
			//if ( in.s == NULL )
			//	printf( "null s. %d %d\n", in.size, in.capacity ) ;
			s = (T *)malloc( sizeof( T ) * capacity ) ;
			memcpy( s, in.s, sizeof( T ) * capacity ) ;
		}
		else 
			s = NULL ;
		inc = in.inc ;
		maxInc = in.maxInc ;
	}
	
	SimpleVector& operator=( const SimpleVector &in )
	{
		if ( this != &in )
		{
			if ( s != NULL )
				free( s ) ;
			size = in.size ;
			capacity = in.capacity ;

			if ( capacity > 0 )
			{
				//s = in.s ;
				s = (T *)malloc( sizeof( T ) * capacity ) ;
				memcpy( s, in.s, sizeof( T ) * capacity ) ;
			}
			else 
				s = NULL ;

			inc = in.inc ;
			maxInc = in.maxInc ;
		}
		return *this ;
	}

	~SimpleVector()
	{
		if ( s != NULL )
			free( s ) ;
		capacity = 0 ;
		size = 0 ;
	}
	
	void Release()
	{
		if ( s != NULL )
			free( s ) ;
		s = NULL ;
		size = capacity = 0 ;
	}

	void Reserve( int sz )
	{
		if ( s != NULL )
			free( s ) ;
		s = (T *)malloc( sizeof( T ) * sz ) ;
		size = 0 ;
		capacity = sz ;
		inc = sz ;

		if ( maxInc > 0 && inc > maxInc )
			inc = maxInc ;
	}


	int PushBack( const T &in )	
	{
		if ( size == capacity )
		{
			//int tmp = capacity ;
			capacity += inc ;
			inc *= 2 ;
			if ( maxInc > 0 && inc > maxInc )
				inc = maxInc ;
			if ( size == 0 )
				s = (T *)malloc( sizeof( T ) * capacity ) ;
			else
				s = (T *)realloc( s, sizeof( T ) * capacity ) ;
			if ( s == NULL ) 
			{
				fprintf( stderr, "%s: Failed to allocate memory.\n", __func__ ) ;
				exit( 1 ) ;
			}
		}
		s[ size ] = in ;
		++size ;
		return size ;
	}

	int PushBack( const SimpleVector<T> &in )
	{
		int newsize = size + in.size ;
		if ( newsize > capacity )
		{
			//int tmp = capacity ;
			capacity = newsize + inc ;
			inc *= 2 ;
			if ( maxInc > 0 && inc > maxInc )
				inc = maxInc ;
			if ( size == 0 )
				s = (T *)malloc( sizeof( T ) * capacity ) ;
			else
				s = (T *)realloc( s, sizeof( T ) * capacity ) ;
			if ( s == NULL ) 
			{
				fprintf( stderr, "%s: Failed to allocate memory.\n", __func__ ) ;
				exit( 1 ) ;
			}
		}
		memcpy( s + size, in.s, sizeof( T ) * in.size ) ;
		size = newsize ;
		return size ;
	}

	T PopBack()
	{
		if ( size == 0 )
		{
			fprintf( stderr, "%s: empty array.\n", __func__ ) ;
			exit( 1 ) ;
		}
		--size ;
		return s[size] ;
	}
	
	int GetInc()
	{
		return inc ;
	}

	void SetInc( int in )
	{
		inc = in ;
	}

	void SetMaxInc( int in )
	{
		maxInc = in ;
	}
	int GetMaxInc()
	{
		return maxInc ;
	}
	int Size()
	{
		return size ;
	}

	int Resize( int s ) 
	{
		size = s ;
		return size ;
	}

	int Capacity()
	{
		return capacity ;
	}

	T &Get( int i )
	{
		if ( i >= size )
		{
			fprintf( stderr, "%s: Access out of the vector.\n", __func__ ) ;
			exit( 1 ) ;
		}
		return s[i] ;
	}

	T &operator[]( int i ) const
	{
		/*if ( i >= size )
		{
			printf( "ERROR\n" ) ;
		}*/
		//assert( i < size ) ;
		/*if ( i >= size )
		{
			fprintf( stderr, "%s: Access out of the vector.\n", __func__ ) ;
			exit( 1 ) ;
		}*/
		return s[i] ;
	}
	
	// Return how many element left.
	int Remove( int ind )
	{
		int i ;
		if ( ind >= size )
		{
			fprintf( stderr, "%s: Access out of the vector.\n", __func__ ) ;
			exit( 1 ) ;
		}

		//if ( size == 1 )
		//	return 0 ;
		for ( i = ind ; i < size - 1 ; ++i )
			s[i] = s[i + 1] ;
		--size ;
		return size ;
	}

	// Allocate less memory. 
	int Shrink()
	{
		if ( size < capacity / 4 )
		{
			capacity /= 2 ;
			inc = capacity ;
			if ( inc > maxInc )
				inc = maxInc ;
			s = (T *)realloc( s, sizeof( T ) * capacity ) ;				
		}
		return capacity ;
	}

	void Clear()
	{
		size = 0 ;
	}

	void QSort( int (*compare)(const void*,const void*) )
	{
		qsort( s, size, sizeof( T ), compare ) ;
	}

	int BinarySearch( const T &v )
	{
		int l, r, m ;
		l = 0 ; 
		r = size - 1 ;

		while ( l <= r )
		{
			m = ( l + r ) / 2 ;
			if ( s[m] == v )
				return m ;
			else if ( s[m] < v )
				l = m + 1 ;
			else
				r = m - 1 ;

		}
		return l - 1 ; // Should be between  l - 1 and l
	}

	void Destroy()
	{
		if ( s != NULL )
			free( s ) ;
		s = NULL ;
		size = capacity = 0 ;
		inc = 1 ;
	}

	void Overwrite( const SimpleVector<T> &in )
	{
		if ( s != NULL )
			free( s ) ;
		s = NULL ;
		if ( in.s != NULL )
			s = (T *)malloc( sizeof( T ) * in.capacity ) ;
		size = in.size ;
		capacity = in.capacity ;
		inc = in.inc ;
		int i ;
		for ( i = 0 ; i < size ; ++i )
			s[i] = in.s[i] ;
	}

	void Reverse()
	{
		int i, j ;
		T tmp ;
		for ( i = 0, j = size - 1 ; i < j ; ++i, --j )
		{
			tmp = s[j] ;
			s[j] = s[i] ;
			s[i] = tmp ;
		}
	}
	
	// Expand the array by given size.
	// Does not care about the value in the new allocated space.
	int ExpandBy( int expandSize )
	{
		int newSize = size + expandSize ;
		if ( newSize <= capacity )
		{
			size = newSize ;
		}
		else
		{
			//int tmp = capacity ;
			capacity = newSize + inc ;
			inc *= 2 ;
			if ( maxInc > 0 && inc > maxInc )
				inc = maxInc ;
			if ( size == 0 )
				s = (T *)malloc( sizeof( T ) * capacity ) ;
			else
				s = (T *)realloc( s, sizeof( T ) * capacity ) ;
			if ( s == NULL ) 
			{
				fprintf( stderr, "%s: Failed to allocate memory.\n", __func__ ) ;
				exit( 1 ) ;
			}
			size = newSize ;
		}
		return size ;
	}

	int ExpandTo( int newSize )
	{
		return ExpandBy( newSize - size ) ;
	}

	void ShiftRight( int shift )
	{
		size = ExpandBy( shift ) ;
		int i ;

		for ( i = size - 1 ; i >= shift ; --i )
			s[i] = s[i - shift] ;
		return ;
	}
	
	// Set the content to zero in the range
	void SetZero( int start, int len )
	{
		memset( s + start, 0, sizeof( T ) * len ) ;
	}

	T *BeginAddress()
	{
		return s ;
	}
	T *EndAddress() 
	{
		return s + size ;
	}
} ;

#endif
