// Implementation of statistical tests
#ifndef _LSONG_STATSTESTS_HEADER
#define _LSONG_STATSTESTS_HEADER

#include "SimpleVector.h"

// Return true if they are different
bool SignTest( SimpleVector<int> &sign )
{
	int size = sign.Size() ;
	int plusCnt = 0 ;
	int i ;
	for ( i = 0 ; i < size ; ++i )
		if ( sign[i] > 0 )
			++plusCnt ;
	printf( "%d %d\n", plusCnt, size ) ;
	// TODO: really implement it.
	if ( plusCnt <= 0.7 * size && plusCnt >= 0.3 * size )
		return true ;
	return false ;
}

#endif
