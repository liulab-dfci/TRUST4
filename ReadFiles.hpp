// The class for reading reads from a file
#ifndef _MOURISL_READS
#define _MOURISL_READS

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "defs.h"
#include "kseq.h"

#define MAX_READ_FILE 100

KSEQ_INIT( gzFile, gzread ) ;

struct _Read
{
	char *id ;
	char *seq ;
	char *qual ;
} ;

class ReadFiles
{
	private:
		gzFile gzFp[MAX_READ_FILE] ;
		kseq_t *inSeq[MAX_READ_FILE] ;
		bool hasMate[MAX_READ_FILE] ;
		int FILE_TYPE[MAX_READ_FILE] ; // 0-FASTA, 1-FASTQ
		int fpUsed ;
		int currentFpInd ;

		void GetFileName( char *in, char *out ) 
		{
			int i, j ;
			int len = (int)strlen( in ) ;
			for ( i = len ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
				;
			for ( j = len ; j >= 0 && in[j] != '/' ; --j )
				;
			if ( i >= 0 && in[i] == '.' )
			{
				in[i] = '\0' ;
				strcpy( out, in + j + 1 ) ;
				in[i] = '.' ;
			}
			else
			{
				strcpy( out, in + j + 1 ) ;
			}
		}

	public:
		char *id ;
		char *seq ;
		char *qual ;

		ReadFiles(): fpUsed(0), currentFpInd(0)
		{
			id = seq = qual = NULL ;
		}

		~ReadFiles()
		{
			int i ;
			if ( id != NULL )
				free( id ) ;
			if ( seq != NULL )
				free( seq ) ;
			if ( qual != NULL )
				free( qual ) ;

			
			for ( i = 0 ; i < fpUsed ; ++i )
			{
				kseq_destroy( inSeq[i] ) ;
				gzclose( gzFp[i] ) ;
			}

		}


		void AddReadFile( char *file, bool fileHasMate )
		{
			if ( fpUsed >= MAX_READ_FILE )
			{
				fprintf( stderr, "The number of read files exceeds the limit %d.\n", MAX_READ_FILE ) ;
				exit( 1 ) ;
			}
			char buffer[1024], fileName[1024] ;

			gzFp[ fpUsed ] = gzopen( file, "r" ) ;
			inSeq[ fpUsed ] = kseq_init( gzFp[ fpUsed ] ) ;

			kseq_read( inSeq[ fpUsed ] ) ;
			if ( inSeq[ fpUsed ]->qual.l == 0 )
			{
				FILE_TYPE[ fpUsed ] = 0 ;
				//qual[0] = '\0' ;
			}
			else 
			{
				FILE_TYPE[ fpUsed ] = 1 ;
			}
			/*else
			{
				fprintf( stderr, "\"%s\"'s format is wrong.\n", file ) ;
				exit( 1 ) ;
			}*/
			
			//printf( "%s %s\n", inSeq[fpUsed]->name.s, inSeq[fpUsed]->comment.s ) ;
			gzrewind( gzFp[ fpUsed ]) ;
			kseq_rewind( inSeq[ fpUsed] ) ;
			//kseq_read( inSeq[ fpUsed ] ) ;
			//printf( "%s %s\n", inSeq[fpUsed]->name.s, inSeq[fpUsed]->comment.s ) ;
			//gzrewind( gzFp[ fpUsed ]) ;
			//kseq_rewind( inSeq[ fpUsed] ) ;
			hasMate[ fpUsed ] = fileHasMate ;	
			++fpUsed ;
		}


		bool HasQuality()
		{
			return ( FILE_TYPE[ currentFpInd ] != 0  ) ;
		}

		void Rewind() 
		{
			int i ;
			for ( i = 0 ; i < fpUsed ; ++i )
			{
				gzrewind( gzFp[i] ) ;
				kseq_rewind( inSeq[i] ) ;
			}
			
			if ( id != NULL )
				free( id ) ;
			if ( seq != NULL )
				free( seq ) ;
			if ( qual != NULL )
				free( qual ) ;
			id = seq = qual = NULL ;
			currentFpInd = 0 ;
		}

		int Next() 
		{
			//int len ;
			//char buffer[2048] ;
			while ( currentFpInd < fpUsed && ( kseq_read( inSeq[ currentFpInd ] ) < 0 ) )
			{
				++currentFpInd ;
			}
			if ( currentFpInd >= fpUsed )
				return 0 ;
			/*printf( "%s %s\n", id, inSeq[currentFpInd ]->comment.s ) ;
			if ( inSeq[currentFpInd]->comment.l )
			{
				id[ inSeq[currentFpInd ]->name.l ] = ' ' ;
				strcpy( &id[ inSeq[currentFpInd]->name.l + 1], inSeq[ currentFpInd]->comment.s ) ;
			}*/
			if ( id != NULL )	
				free( id ) ;
			if ( seq != NULL )
				free( seq ) ;
			if ( qual != NULL )
				free( qual ) ;

			id = strdup( inSeq[ currentFpInd ]->name.s ) ;
			int len = strlen( id ) ;
			if ( ( id[len - 1] == '1' || id[len - 1] == '2' )
					&& id[len - 2] == '/' )
			{
				id[len - 2] = '\0' ;
			}
			
			seq = strdup( inSeq[ currentFpInd ]->seq.s ) ;
			if ( inSeq[ currentFpInd ]->qual.l )
				qual = strdup( inSeq[ currentFpInd]->qual.s ) ;
			else
				qual = NULL ;

			return 1 ;
		}

		int NextWithBuffer( char **id, char **seq, char **qual, bool removeReturn = true, bool stopWhenFileEnds = false ) 
		{
			//int len ;
			//char buffer[2048] ;
			while ( currentFpInd < fpUsed && ( kseq_read( inSeq[ currentFpInd ] ) < 0 ) )
			{
				++currentFpInd ;
				if ( stopWhenFileEnds )
					return -1 ;
			}
			if ( currentFpInd >= fpUsed )
				return 0 ;
			
			if ( *id != NULL )	
				free( *id ) ;
			if ( *seq != NULL )
				free( *seq ) ;
			if ( *qual != NULL )
				free( *qual ) ;

			*id = strdup( inSeq[ currentFpInd ]->name.s ) ;
			*seq = strdup( inSeq[ currentFpInd ]->seq.s ) ;
			/*if ( removeReturn )
			{
				int i ;
				for ( i = strlen( *seq ) - 1 ; i >= 0 ; --i )
					if ( (*seq)[i] != '\n' )
						break ;
				(*seq)[i + 1] = '\0' ;
			}*/
			if ( inSeq[ currentFpInd ]->qual.l )
				*qual = strdup( inSeq[ currentFpInd]->qual.s ) ;
			else
				*qual = NULL ;
			
			return 1 ;
		}

		// Get a batch of reads, it terminates until the buffer is full or 
		// the file ends.
		int GetBatch( struct _Read *readBatch, int maxBatchSize, int &fileInd, bool trimReturn, bool stopWhenFileEnds )
		{
			int batchSize = 0 ;
			while ( batchSize < maxBatchSize ) 
			{
				int tmp = NextWithBuffer( &readBatch[ batchSize].id, &readBatch[batchSize].seq,
							&readBatch[batchSize].qual, trimReturn, stopWhenFileEnds ) ;
				if ( tmp == -1 && batchSize > 0 )
				{
					--currentFpInd ;
					fileInd = currentFpInd ;
					return batchSize ; // Finished read current file.
				}
				else if ( tmp == -1 && batchSize == 0 )
					continue ; // The current read file is empty
				else if ( tmp == 0 && batchSize == 0 )
				{
					fileInd = currentFpInd ;
					return 0 ; // Finished reading	
				}

				++batchSize ;
			}

			fileInd = currentFpInd ;
			return batchSize ;
		}

		int GetFpUsed()
		{
			return fpUsed ;
		}
} ;

// The class handling read in a batch of reads
/*class ReadBatch
{
} ;*/

#endif
