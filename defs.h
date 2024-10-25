#ifndef _LSONG_DEFS_HEADER
#define _LSONG_DEFS_HEADER

#include <stdint.h>

//#define DEBUG

extern int nucToNum[26] ; 
extern char numToNuc[26] ;

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define ABS(x) (((x)>(0))?(x):(-(x)))

typedef int index_t ; 

#define MAX_SEG_COUNT 127
struct _pair
{
	int a, b ;
} ;

struct _pair64
{
	int64_t a, b ;
} ;

struct _triple
{
	int a, b, c ;
	bool operator<(const struct _triple &other)
	{
		if (a != other.a)
			return a < other.a ;
		else if (b != other.b)
			return b < other.b ;
		else
			return c < other.c ;
	}
} ;
/*struct _pair_b64
{
	int a ;
	uint64_t b ;
} ;*/

#endif
