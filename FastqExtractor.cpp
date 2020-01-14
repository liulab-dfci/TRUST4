#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <pthread.h>

#include "SeqSet.hpp"
#include "ReadFiles.hpp"

char usage[] = "./bam-extractor [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t-u STRING: path to single-end read file\n"
		"\t\tor\n"
		"\t-1 STRING -2 STRING: path to paired-end read files"
		"Optional:\n"
		"\t-o STRING: prefix to the output file\n"
		"\t-t INT: number of threads (default: 1)\n" ;

static const char *short_options = "f:u:1:2:o:t" ;
static struct option long_options[] = {
			{ (char *)0, 0, 0, 0} 
			} ;

char buffer[100001] ;
char seqBuffer[100001] ;
char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

struct _threadArg 
{
	struct _Read *readBatch, *readBatch2 ;

	int batchSize ;
	int batchUsed ;
	char *buffer ;

	SeqSet *refSet ;
	
	int tid ;
} ;

void PrintLog( const char *fmt, ... )
{
	va_list args ;
	va_start( args, fmt ) ;
	vsprintf( buffer, fmt, args ) ;

	time_t mytime = time(NULL) ;
	struct tm *localT = localtime( &mytime ) ;
	char stime[500] ;
	strftime( stime, sizeof( stime ), "%c", localT ) ;
	fprintf( stderr, "[%s] %s\n", stime, buffer ) ;
}

bool IsLowComplexity( char *seq )
{
	int cnt[5] = {0, 0, 0, 0, 0} ;
	int i ;
	for ( i = 0 ; seq[i] ; ++i )
	{
		if ( seq[i] == 'N' )
			++cnt[4] ;
		else
			++cnt[ nucToNum[ seq[i] - 'A' ] ] ;
	}

	if ( cnt[0] >= i / 2 || cnt[1] >= i / 2 || cnt[2] >= i / 2 || cnt[3] >= i / 2 || cnt[4] >= i / 10 )
		return true ;

	int lowCnt = 0 ; 
	for ( i = 0 ; i < 4 ; ++i )
		if ( cnt[i] <= 2 )
			++lowCnt ;
	if ( lowCnt >= 2 )
		return true ;
	return false ;
}

int IsGoodCandidate( char *read, char *buffer, SeqSet *refSet )
{
	if ( !IsLowComplexity( read ) && refSet->HasHitInSet( read, buffer ) )
		return 1 ;
	return 0 ;
}

void OutputSeq( FILE *fp, const char *name, char *seq, char *qual )
{
	if ( qual != NULL )
		fprintf( fp, "@%s\n%s\n+\n%s\n", name, seq, qual ) ;
	else
		fprintf( fp, ">%s\n%s\n\n", name, seq ) ;
}

void *ProcessReads_Thread( void *pArg )
{
	int i ;	
	struct _threadArg &arg = *((struct _threadArg *)pArg ) ;
	for ( i = 0 ; i < arg.batchSize ; ++i )
	{
		if ( i % arg.tid != 0 )
			continue ;

		int goodCandidate = 0 ;
		if ( IsGoodCandidate( arg.readBatch[i].seq, arg.buffer, arg.refSet ) )
			++goodCandidate ;
		
		if ( !goodCandidate && arg.readBatch2 && IsGoodCandidate( arg.readBatch2[i].seq, arg.buffer, arg.refSet ) )
			++goodCandidate ;

		if ( !goodCandidate ) 
			arg.readBatch[i].id[0] = '\0' ;
	}

	pthread_exit( NULL ) ;
}

int main( int argc, char *argv[] )
{
	int i ;
	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ;
	option_index = 0 ;
	FILE *fpRef = NULL ;
	char prefix[1024] = "toassemble" ;
	int kmerLength = 9  ;
	SeqSet refSet( kmerLength ) ;
	int threadCnt = 1 ;
	
	ReadFiles reads ;
	ReadFiles mateReads ;
	bool hasMate = false ;
	
	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			fpRef = fopen( optarg, "r" ) ;
			refSet.InputRefFa( optarg ) ;	
		}
		else if ( c == 'o' )
		{
			strcpy( prefix, optarg ) ;
		}
		else if ( c == '1' )
		{
			reads.AddReadFile( optarg, true ) ;
			hasMate = true ;
		}
		else if ( c == '2' )
		{
			mateReads.AddReadFile( optarg, true ) ;
			hasMate = true ;
		}
		else if ( c == 'u' )
		{
			reads.AddReadFile( optarg, false ) ;
		}
		else if ( c == 't' )
		{
			threadCnt = atoi( optarg ) ;
		}
		else
		{
			fprintf( stderr, "Unknown parameter %s\n", optarg ) ;
			return EXIT_FAILURE ;
		}
	}


	if ( fpRef == NULL )
	{
		fprintf( stderr, "Need to use -f to specify the receptor genome sequence.\n" ) ;
		return EXIT_FAILURE ;
	}
	
	PrintLog( "Start to extract candidate reads from read files." ) ;
	
	int hitLenRequired = 31 ;
	if ( !hasMate )
		hitLenRequired = 27 ;
	int len = 0 ;
	for ( i = 0 ; i < 1000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		len += strlen( reads.seq ) ;
	}
	if ( len / (i * 5) > hitLenRequired )
		hitLenRequired = len / (i * 5) ;
	refSet.SetHitLenRequired( hitLenRequired ) ;
	
	FILE *fp1 = NULL ;
	FILE *fp2 = NULL ;
	if ( !hasMate )
	{
		sprintf( buffer, "%s.fq", prefix ) ;
		fp1 = fopen( buffer, "w" ) ;
	}
	else
	{
		sprintf( buffer, "%s_1.fq", prefix ) ;
		fp1 = fopen( buffer, "w" ) ;
		sprintf( buffer, "%s_2.fq", prefix ) ;
		fp2 = fopen( buffer, "w" ) ;
	}
	
	pthread_t *threads ;
	struct _threadArg *threadArgs ;
	pthread_attr_t attr ;
	void *pthreadStatus ;

	if ( threadCnt == 1 )
	{
		while ( reads.Next() )
		{
			if ( hasMate && !mateReads.Next() )
			{
				fprintf( stderr, "The two mate-pair read files have different number of reads.\n" ) ;
				exit( 1 ) ;
			}

			int goodCandidate = 0 ;
			
			if ( IsGoodCandidate( reads.seq, seqBuffer, &refSet ) )
				++goodCandidate ;
			if ( !goodCandidate && hasMate && IsGoodCandidate( mateReads.seq, seqBuffer, &refSet ) )
				++goodCandidate ;
		
			if ( goodCandidate )
			{
				OutputSeq( fp1, reads.id, reads.seq, reads.qual ) ;
				if ( hasMate )
					OutputSeq( fp2, reads.id, reads.seq, reads.qual ) ;
			}
		}
	}
	else
	{
		int maxBatchSize = 512 * threadCnt ;
		int batchSize ;
		
		struct _Read *readBatch = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		struct _Read *readBatch2 = NULL ;
		if ( hasMate )
			readBatch2 = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		int fileInd1, fileInd2 ;
		
		struct _threadArg *args = (struct _threadArg *)malloc( sizeof( struct _threadArg ) * threadCnt ) ;

		for ( i = 0 ; i < threadCnt ; ++i )
		{
			args[i].tid = i ;
			args[i].readBatch = readBatch ;
			args[i].readBatch2 = readBatch2 ;

			args[i].buffer = (char *)malloc( sizeof( char ) * 10001 ) ;
			args[i].refSet = &refSet ;
		}

		while ( 1 )
		{
			batchSize = reads.GetBatch( readBatch, maxBatchSize, fileInd1, true, true ) ;
			if ( hasMate )
			{
				int tmp = mateReads.GetBatch( readBatch2, maxBatchSize, fileInd2, true, true ) ;
				if ( tmp != batchSize )
				{
					fprintf( stderr, "The two mate-pair read files have different number of reads.\n" ) ;
					exit ( 1 ) ;
				}
			}

			if ( batchSize == 0 )
				break ; 
			//printf( "batchSize=%d\n", batchSize ) ;

			for ( i = 0 ; i < threadCnt ; ++i )
			{
				args[i].batchSize = batchSize ;
				pthread_create( &threads[i], &attr, ProcessReads_Thread, (void *)&args[i] ) ;	
			}

			for ( i = 0 ; i < threadCnt ; ++i )
				pthread_join( threads[i], &pthreadStatus ) ;

			for ( i = 0 ; i < batchSize ; ++i )
			{
				if ( readBatch[i].id[0] == '0' )
					continue ;
				OutputSeq( fp1, readBatch[i].id, readBatch[i].seq, readBatch[i].qual ) ;
				if ( readBatch2 != NULL )
					OutputSeq( fp2, readBatch2[i].id, readBatch2[i].seq, readBatch2[i].qual ) ;
			}
		}
		
		// Release memory.
		for ( i = 0 ; i < batchSize ; ++i )
		{
			int j ;
			free( readBatch[i].id ) ;
			free( readBatch[i].seq ) ;
			if ( readBatch[i].qual )
				free( readBatch[i].qual ) ;
			if ( readBatch2 != NULL )
			{
				free( readBatch2[i].id ) ;
				free( readBatch2[i].seq ) ;
				if ( readBatch2[i].qual )
					free( readBatch2[i].qual ) ;
			}
		}
		for ( i = 0 ; i < threadCnt ; ++i )
			free( args[i].buffer ) ;
		free( readBatch ) ;
		if ( hasMate )
			free( readBatch2 ) ;
	}
	fclose( fp1 ) ;
	if ( hasMate )
		fclose( fp2 ) ;
	return 0 ;
}
