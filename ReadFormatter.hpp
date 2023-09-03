#ifndef _MOURISL_READ_FORMATTER
#define _MOURISL_READ_FORMATTER

#include <vector>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BufferManager.hpp"

enum {
  FORMAT_READ1,
  FORMAT_READ2,
  FORMAT_BARCODE,
  FORMAT_UMI,
  FORMAT_CATEGORY_COUNT
} ;

struct _segInfo
{
  int start ;
  int end ;
  int strand ; // -1: minus, 1:posiive
  
	bool operator<( const struct _segInfo &b )	const
  {
    return start < b.start ; 
  }
} ;

// Parse for read1, read2, barcode and UMI.
// Current implementation is not thread-safe.
class ReadFormatter
{
private:
  BufferManager<char> _buffers ;
  char _compChar[256] ;

  std::vector<struct _segInfo> _segs[FORMAT_CATEGORY_COUNT] ;
  
  // Return false if it fails to parse the format string.
  bool ParseFormatStringAndAppendEffectiveRange(const char *s, int len) {
    int i;
    int j = 0;  // start, end, strand section
    char buffer[20];
    int blen = 0;
    
    struct _segInfo seg ;
    if (s[2] != ':')
      return false ;
    int category = 0 ;
    if (s[0] == 'r' && s[1] == '1') {
      category = FORMAT_READ1 ;
    } else if (s[0] == 'r' && s[1] == '2') {
      category = FORMAT_READ2 ;
    } else if (s[0] == 'b' && s[1] == 'c') {
      category = FORMAT_BARCODE ;
    } else if (s[0] == 'u' && s[1] == 'm') {
      category = FORMAT_UMI ;
    } else {
      return false ;
    }
    
    seg.strand = 1 ;
    for (i = 3; i <= len; ++i) {
      if (i == len || s[i] == ':') {
        buffer[blen] = '\0';
        if (j == 0) {
          seg.start = atoi(buffer) ;
        } else if (j == 1) {
          seg.end = atoi(buffer) ;
        } else {
          seg.strand = (buffer[0] == '+' ? 1 : -1);
        }

        blen = 0;
        if (i < len && s[i] == ':') {
          ++j;
        }
      } else {
        buffer[blen] = s[i];
        ++blen;
      }
    }
    if (j >= 3 || j < 1) {
      return false;
    }
    
    _segs[category].push_back(seg) ;
    return true;
  }

  void ReverseBuffer(char *buffer, int len)
  {
    int i, j ;
    for (i = 0, j = len - 1 ; i < j ; ++i, --j )
    {
      char tmp = buffer[i] ;
      buffer[i] = buffer[j] ;
      buffer[j] = tmp ;
    }
  }

  void ComplementBuffer(char *buffer, int len)
  {
    int i ;
    for (i = 0 ; i < len ; ++i)
      buffer[i] = _compChar[ (int)buffer[i] ] ;
  }

public:
  ReadFormatter() {
    int i ;
    for (i = 0 ; i < 256 ; ++i)
      _compChar[i] = 'N' ;
    _compChar['A'] = 'T' ;
    _compChar['C'] = 'G' ;
    _compChar['G'] = 'C' ;
    _compChar['T'] = 'A' ;
  } 

  ~ReadFormatter() {
  }

  void AllocateBuffers(int bufferCnt)
  {
    _buffers.Init(bufferCnt) ;
  }

  void Init(const char *formatStr) {
    int i, j;
    if (_buffers.GetBufferCount() == 0)
      AllocateBuffers(2) ;
    for (i = 0 ; formatStr[i] ; ) {
      for (j = i ; formatStr[j] && formatStr[j] != ';' && formatStr[j] != ',' ; ++j)
        ;
      
      if (!ParseFormatStringAndAppendEffectiveRange(formatStr + i, j - i))
      {
        fprintf(stderr, "Format description error in %s\n", formatStr) ;
        exit(1) ;
      }

      if (formatStr[j])
        i = j + 1 ;
      else
        i = j ;
    }
    
    // Sort the order in each specification
    // It seems there are applications 
    //for (i = 0 ; i < FORMAT_CATEGORY_COUNT ; ++i)
    //  std::sort(_segs[i].begin(), _segs[i].end()) ;
  }

  void AddSegment(int start, int end, int strand, int category)
  {
    struct _segInfo ns ;
    ns.start = start ;
    ns.end = end ;
    ns.strand = strand ;
    _segs[category].push_back(ns) ;
    //std::sort(_segs[ category ].begin(), _segs[ category ].end()) ;
    
    if (_buffers.GetBufferCount() == 0)
      AllocateBuffers(2) ;
  }

  int NeedExtract(int category)
  {
    int size = _segs[category].size() ;
    if (size == 0)
      return 0 ;
    else if (size == 1)
    {
      if (_segs[category][0].start == 0  
          && _segs[category][0].end == -1
          && _segs[category][0].strand == 1)
        return 0 ;
    }
    return 1 ;
  }

  // needComplement=true: reverse complement. Otherwise, just reverse
  // retSeqWhenNoExtraction: when needextract==false, return seq instead of buffer
  // The outside program can modify the buffer.
  char* Extract(char *seq, int category, bool needComplement, bool retSeqWhenNoExtraction, int bufferId = 0)
  {
    int len = strlen(seq) ;
    int i, j, k ;
    const std::vector<_segInfo> &seg = _segs[category] ;
    int segSize = seg.size() ;
    int strand = 1 ;
    char *buffer = _buffers.Get(bufferId, len + 1) ;
    
    if (!NeedExtract(category))
    {
      if (retSeqWhenNoExtraction)
        return seq ;
      else
      {
        strcpy(buffer, seq) ;
        return buffer ;
      }
    }

    i = 0 ;
    for (k = 0 ; k < segSize ; ++k)
    {
      int start = seg[k].start ;
      int end = seg[k].end ;
      
      if (start < 0)
        start = len + start ;

      if (end >= len)
        end = len - 1 ;
      else if (end < 0)
        end = len + end ;

      for (j = start ; j <= end ; ++j)
      {
        buffer[i] = seq[j] ;
        ++i ;
      }
      if (seg[k].strand == -1)
        strand = -1 ;
    }
    buffer[i] = '\0' ;

    if (strand == -1)
    {
      ReverseBuffer(buffer, i) ;
      if (needComplement)
        ComplementBuffer(buffer, i) ;
    }
    return buffer ;
  }
} ;

#endif
