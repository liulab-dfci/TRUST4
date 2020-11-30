#ifndef _MOURISL_BARCODECORRECTOR_HEADER
#define _MOURISL_BARCODECORRECTOR_HEADER

#include <stdio.h>
#include "ReadFiles.hpp"
#include "SimpleVector.hpp"
#include "defs.h"

struct _trie
{
	struct _trie *next[4] ;

	bool end ;
	int count ;
} ;

class Trie
{
private:
	struct _trie head ;	
	char nucToNum[26] ; 
	
	void ReleaseMemory(struct _trie *node)
	{
		for (int i = 0 ; i < 4 ; ++i)
		{
			if (node->next[i] != NULL)
			{
				ReleaseMemory(node->next[i]) ;
				delete node->next[i] ;
				node->next[i] = NULL ;
			}
		}
	}

	void InitNode(struct _trie &node)
	{
		for (int i = 0 ; i < 4 ; ++i)
			node.next[i] = NULL ;
		node.end = false ;
		node.count = 0 ;
	}
public:
	Trie()
	{
		InitNode(head) ;

		memset(nucToNum, -1, sizeof(nucToNum)) ;
		nucToNum['A' - 'A'] = 0 ;
		nucToNum['C' - 'A'] = 1 ;
		nucToNum['G' - 'A'] = 2 ;
		nucToNum['T' - 'A'] = 3 ;
	}

	~Trie()
	{
		ReleaseMemory(&head) ;
	}

	void Insert(char *s, int weight)
	{
		int i ;
		struct _trie *p = &head ;
		for (i = 0 ; s[i] ; ++i)
			if (nucToNum[s[i] - 'A'] == -1)
				return ;
	
		for (i = 0 ; s[i] ; ++i)
		{
			int tag = nucToNum[s[i] - 'A'] ;
			if (p->next[tag] == NULL) 
			{
				p->next[tag] = new struct _trie ;
				InitNode(*(p->next[tag])) ;
			}
			p = p->next[tag] ;
		}
		p->end = true ;
		p->count += weight ;
	}

	int SearchAndUpdate(char *s, int weight) // Return the count after the update: -1 not found.
	{
		int i ;
		struct _trie *p = &head ;
		for (i = 0 ; s[i] ; ++i)
			if (nucToNum[s[i] - 'A'] == -1)
				return -1 ;
		
		for (i = 0 ; s[i] ; ++i)
		{
			int tag = nucToNum[s[i] - 'A'] ;
			if (p->next[tag] == NULL)
				return -1 ;
			p = p->next[tag] ;
		}
		p->count += weight ;
		return p->count ;
	}
} ;

class BarcodeCorrector
{
private:
	Trie barcodeFreq ;

	void ReverseComplement( char *rcSeq, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
		{
			if ( seq[len - 1 - i] != 'N' )
				rcSeq[i] = numToNuc[ 3 - nucToNum[seq[len - 1 - i] - 'A'] ];
			else
				rcSeq[i] = 'N' ;
		}
		rcSeq[i] = '\0' ;
	}
	
	// Put the reformatted raw barcode into buffer
	void FormatBarcode(char *raw, int start, int end, bool revcomp, char *buffer)
	{
		if ( start == 0 && end == -1 && revcomp == false )
			strcpy(buffer, raw) ;
		else
		{
			int i ;
			int s = start ;
			int e = ( end == -1 ? strlen( raw ) - 1 : end ) ;

			if ( revcomp == false )
			{
				for ( i = s ; i <= e ; ++i )
					buffer[i - s] = raw[i] ;
				buffer[i - s] = '\0' ;
			}
			else
				ReverseComplement( buffer, raw + s, e - s + 1 ) ;
		}
	}
public:
	BarcodeCorrector() {}
	~BarcodeCorrector() {}
	
	void SetWhitelist(char *whitelist)
	{
		FILE *fp = fopen(whitelist, "r") ;
		char buffer[256] ;
		while ( fscanf(fp, "%s", buffer) != EOF)
			barcodeFreq.Insert(buffer, 1) ;
		fclose(fp) ;
	}

	int CollectBackgroundDistribution(ReadFiles &barcodeFile, int start, int end, bool revcomp, int caseCnt = 2000000) 
	{
		int readCnt = 0 ;
		char buffer[256] ;
		while (barcodeFile.Next())
		{	
			FormatBarcode(barcodeFile.seq, start, end, revcomp, buffer) ;
			barcodeFreq.SearchAndUpdate(buffer, 1) ;
			readCnt += 1 ;
			if (readCnt >= caseCnt) 
				break ;
		}
		barcodeFile.Rewind() ;
	}
	
	// return: -1: could not correct. 0: no need for correction. 1: corrected
	int Correct(char *barcode, char *qual)
	{
		char buffer[256] ;
		if (barcodeFreq.SearchAndUpdate(barcode, 0) != -1)
		{
			return 0 ;
		}
		else
		{
			int i, j ;
			char testChr[5] = "ACGT" ;
			SimpleVector<struct _triple> records ;

			strcpy( buffer, barcode) ;
			for ( i = 0 ; barcode[i] ; ++i )
			{
				for (j = 0 ; j < 4; ++j)
				{
					if (testChr[j] == barcode[i])
						continue ;
					buffer[i] = testChr[j] ;
					int cnt = barcodeFreq.SearchAndUpdate(buffer, 0) ;
					buffer[i] = barcode[i] ;

					if (cnt != -1)
					{
						struct _triple nr ;
						nr.a = i ;
						nr.b = j ;
						nr.c = cnt ;

						records.PushBack(nr) ;
					}
				}
			}	

			int bestCnt = -1 ;
			int bestTag = -1 ;
			int bestLowQual = 255 ; // the lowest quality score within the best candidates

			int recordCnt = records.Size() ;
			if (recordCnt == 0)
				return -1 ;

			for ( i = 0 ; i < recordCnt ; ++i )
			{
				if (records[i].c > bestCnt)
				{
					bestCnt = records[i].c ;
					bestTag = i ;
					if (qual != NULL)
						bestLowQual = qual[ records[i].a ] ;
				}
				else if (records[i].c == bestCnt)
				{
					if (qual != NULL && qual[records[i].a] < bestLowQual)
					{
						bestLowQual = qual[records[i].a] ;
						bestTag = i ;
					}
				}
			}
			// If there is a N in the sequence, all the corrections at other bases will be invalid implicitly
			// If there are more than one Ns, all the corrections are also invalid implicitly.
			// So there is no need to worry about that the correction must fix N.
			barcode[ records[bestTag].a ] = testChr[ records[bestTag].b ] ;
			return 1 ;
		}
	}
} ;


#endif
