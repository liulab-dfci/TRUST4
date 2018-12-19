#ifndef _LSONG_DEFS_HEADER
#define _LSONG_DEFS_HEADER

#include <stdint.h>

extern char nucToNum[26] ; 
extern char numToNuc[26] ;

#define EDIT_MATCH 0
#define EDIT_MISMATCH 1
#define EDIT_INSERT 2
#define EDIT_DELETE 3

#define SCORE_MATCH 5
#define SCORE_MISMATCH (-2)
#define SCORE_GAPOPEN (-4)
#define SCORE_GAPEXTEND (-1)
#define SCORE_INDEL (-4)

#define MAX(x,y) (((x)>(y))?(x):(y))
#define ABS(x) (((x)>(0))?(x):(-(x)))

typedef uint32_t index_t ; 

struct _pair
{
	int a, b ;
} ;

/*struct _pair_b64
{
	int a ;
	uint64_t b ;
} ;*/

#endif
