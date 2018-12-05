#ifndef _LSONG_DEFS_HEADER
#define _LSONG_DEFS_HEADER

#include <stdint.h>

extern char nucToNum[26] ; 
extern char numToNuc[26] ;

#define EDIT_MATCH 0
#define EDIT_MISMATCH 1
#define EDIT_INSERT 2
#define EDIT_DELETE 3


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
