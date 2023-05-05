#ifndef _MOURISL_BUFFER_MANAGER
#define _MOURISL_BUFFER_MANAGER

#include <stdio.h>

// Management of multiple buffers
// The main feature is to easily adjust size
template <class T>
class BufferManager
{
private:
  T **_buffers ;
  size_t *_bufferSize ;
  int _bufferCnt ;
public:
  BufferManager() 
  {
    _buffers = NULL ;
    _bufferSize = NULL ;
    _bufferCnt = 0 ;
  }
  ~BufferManager()
  {
    Free() ;
  }

  int GetBufferCount()
  {
    return _bufferCnt ;
  }
  
  void Free()
  {
    int i ;
    for (i = 0 ; i < _bufferCnt ; ++i)
      delete[] _buffers[i] ;
    delete[] _buffers ;
    delete[] _bufferSize ;
  }
  
  void Init(int cnt)
  {
    Free() ;
    _bufferCnt = cnt ;
    _buffers = new T*[cnt] ;
    _bufferSize = new size_t[cnt] ; 
    for (int i = 0 ; i < cnt ; ++i)
    {
      const int defaultSize = 256 ;
      _buffers[i] = new T[defaultSize] ;
      _bufferSize[i] = defaultSize ;
    }
  }

  // Get the ith buffer, planning use it to hold size elements
  T *Get(int i, size_t size)
  {
    if (size > _bufferSize[i])
    {
      delete[] _buffers[i] ;
      _buffers[i] = new T[size] ;
      _bufferSize[i] = size ;
    }
    return _buffers[i] ;
  }
} ;

#endif
