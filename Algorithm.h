// Methods related to alignment and string operations
#ifndef _LSONG_ALIGNMENT_HEADER
#define _LSONG_ALIGNMENT_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "KmerCode.h"
#include "SimpleVector.h"
#include "defs.h"

void ReverseComplement( char *s, int len ) ;

// We had a kmer hit from pa and pb on the two sequences.
// Then we can decide whether these two sequences are from the strand or not.
bool IsSameStrand( char *sa, index_t pa, char *sb, index_t pb, int kmerLength ) ;

// Is sequence sa and sb similar and overlapped?
// The overlapped region are stored in rangeA and rangeB for sequence A, B
// hits hold the kmer hits offset, a for seqA, b for seqB.
bool IsOverlap( SimpleVector< struct _pair > &hits, bool sameStrand, int minHitRequired, struct _pair &rangeA, struct _pair &rangeB, SimpleVector< struct _pair > *retLIS ) ;

int GlobalAlignment( char *t, int lent, char *p, int lenp, char *align ) ;
int GlobalAlignment_classic( char *t, int lent, char *p, int lenp, char *align ) ;
void VisualizeAlignment( char *t, int lent, char *p, int lenp, char *align ) ;
#endif
