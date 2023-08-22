// The class that read a line froma file
// Does not store the change line symbol
//
#ifndef _MOURISL_FILELINEREADER
#define _MOURISL_FILELINEREADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class FileLineReader
{
private:
	FILE *fp ;
	char c ;
	char *buffer ;
	int bufferSize ;
	char delimit ;
public:
	FileLineReader() 
	{
		fp = NULL ;
		bufferSize = 1000 ;
		buffer = (char *)malloc(sizeof(char) * bufferSize) ;
		buffer[0] = '\0' ;
		delimit = '\n' ;
	}

	~FileLineReader() 
	{
		if (fp != NULL) 
			fclose(fp) ;
		if (buffer != NULL) 
			free(buffer) ;
	}

	int Open(const char *filename)
	{
		fp = fopen(filename, "r") ;
		if (fp == NULL)
			return 0 ;
		return 1 ;
	}

	int IsOpen()
	{
		if (fp == NULL)
			return 0 ;
		return 1 ;
	}

	void Close()
	{
		if (fp != NULL) 
			fclose(fp) ;
		fp = NULL ;
	}

	const char *ReadLine()
	{
		if (feof(fp))
			return NULL ;

		int len = 0 ;	
		while ((c = fgetc(fp)) != EOF)
		{
			buffer[len] = c ;
			++len ;
			if (len >= bufferSize - 1) // not place to hold \0
			{
				bufferSize *= 2 ;
				buffer = (char *)realloc(buffer, sizeof(char) * bufferSize) ;
			}
			if (c == delimit)
			{
				--len ;
				break ;
			}
		}
		buffer[len] = '\0' ;
		if (len == 0 && feof(fp))
			return NULL ;
		return GetLinePtr() ;
	}

	const char *GetLinePtr()
	{
		return (const char*)buffer ;
	}

	FILE *GetFP()
	{
		return (FILE*)fp ;
	}
} ;

#endif
