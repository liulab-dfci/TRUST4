#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <pthread.h>

#include "SeqSet.hpp"
#include "ReadFiles.hpp"
#include "BarcodeCorrector.hpp"
#include "ReadFormatter.hpp"
#include "BarcodeTranslator.hpp"

char usage[] = "./fastq-extractor [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t-u STRING: path to single-end read file\n"
		"\t\tor\n"
		"\t-1 STRING -2 STRING: path to paired-end read files\n"
		"Optional:\n"
		"\t-o STRING: prefix to the output file (default: toassemble)\n"
		"\t-t INT: number of threads (default: 1)\n" 
		"\t--barcode STRING: path to the raw barcode file (default: not used)\n"
		//"\t--barcodeStart INT: the start position of barcode in the barcode sequence (default: 0)\n"
		//"\t--barcodeEnd INT: the end position of barcode in the barcode sequence (default: length-1)\n"
		//"\t--barcodeRevComp: whether the barcode need to be reverse complemented (default: not used)\n"
		"\t--barcodeWhitelist STRING: path to the barcode whitelist (default: not used)\n"
		"\t--barcodeTranslate STRING: path to the barcode translate file (default: not used)\n"
		//"\t--read1Start INT: the start position of sequence in read 1 (default: 0)\n"
		//"\t--read1End INT: the end position of sequence in read 1 (default: length-1)\n"
		//"\t--read2Start INT: the start position of sequence in read 2 (default: 0)\n"
		//"\t--read2End INT: the end position of sequence in read 2 (default: length-1)\n"
		"\t--UMI STRING: path to the raw UMI file (default: not used)\n"
		//"\t--umiStart INT: the start position of barcode in the barcode sequence (default: 0)\n"
		//"\t--umiEnd INT: the end position of barcode in the barcode sequence (default: length-1)\n"
		//"\t--umiRevComp: whether the barcode need to be reverse complemented (default: not used)\n"
		"\t--readFormat STRING: format for read, barcode and UMI files (example: r1:0:-1,r2:0:-1,bc:0:15,um:16:-1 for paired-end files with barcode and UMI)\n"
		;

static const char *short_options = "f:u:1:2:o:t:" ;
static struct option long_options[] = {
			{ "barcode", required_argument, 0, 10000},
			{ "barcodeStart", required_argument, 0, 10001},
			{ "barcodeEnd", required_argument, 0, 10002},
			{ "barcodeRevComp", no_argument, 0, 10003},
			{ "barcodeWhitelist", required_argument, 0, 10004},
			{ "read1Start", required_argument, 0, 10005},
			{ "read1End", required_argument, 0, 10006},
			{ "read2Start", required_argument, 0, 10007},
			{ "read2End", required_argument, 0, 10008},
			{ "UMI", required_argument, 0, 10009},
			{ "umiStart", required_argument, 0, 10010},
			{ "umiEnd", required_argument, 0, 10011},
			{ "umiRevComp", no_argument, 0, 10012},
			{ "readFormat", required_argument, 0, 10013},
			{ "barcodeTranslate", required_argument, 0, 10014},
			{ (char *)0, 0, 0, 0} 
			} ;

char buffer[100001] ;
char seqBuffer[100001] ;
int nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

struct _threadArg 
{
	struct _Read *readBatch, *readBatch2, *barcodeBatch ;
	
	int threadCnt ;
	int batchSize ;
	int batchUsed ;
	char *buffer ;

	SeqSet *refSet ;
	
	BarcodeCorrector *barcodeCorrector ;
	int barcodeStart, barcodeEnd ;
	bool barcodeRevComp ;

	int umiStart, umiEnd ;
	bool umiRevComp ;

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
	if ( !IsLowComplexity( read ) && refSet->HasHitInSet( read ) != 0)
		return 1 ;
	return 0 ;
}

void OutputSeq( FILE *fp, const char *name, char *seq, char *qual, ReadFormatter &readFormatter, int readCategory )
{
	if ( qual != NULL )
		fprintf( fp, "@%s\n%s\n+\n%s\n", name, readFormatter.Extract(seq, readCategory, true, true, 0),
				readFormatter.Extract(qual, readCategory, false, 1)) ;
	else
		fprintf( fp, ">%s\n%s\n", name, readFormatter.Extract(seq, readCategory, true, true, 0) ) ;
}

// Maybe barcode read quality could be useful in future.
void OutputBarcode( FILE *fp, const char *name, char *barcode, char *qual, 
	ReadFormatter &readFormatter, int readCategory, BarcodeCorrector *barcodeCorrector,
	BarcodeTranslator *barcodeTranslator)
{
	if (barcode && barcode[0] != '\0')
	{
		char *bcBuffer = readFormatter.Extract(barcode, readCategory, true, true, 0) ;
		int result = 0 ;
		if ( barcodeCorrector != NULL )	
			result = barcodeCorrector->Correct(bcBuffer, qual ) ;
		if (result >= 0)
		{
			if (barcodeTranslator)
			{
				std::string newbc = barcodeTranslator->Translate(bcBuffer, strlen(bcBuffer)) ;
				fprintf( fp, ">%s\n%s\n", name, newbc.c_str()) ;

			}
			else
				fprintf( fp, ">%s\n%s\n", name, bcBuffer ) ;
		}
		else
			fprintf(fp, ">%s\nmissing_barcode\n", name) ;
	}
	else
		fprintf( fp, ">%s\nmissing_barcode\n", name ) ;
}

void *ProcessReads_Thread( void *pArg )
{
	int i ;	
	struct _threadArg &arg = *((struct _threadArg *)pArg ) ;
	for ( i = 0 ; i < arg.batchSize ; ++i )
	{
		if ( i % arg.threadCnt != arg.tid )
			continue ;
		int goodCandidate = 0 ;
		if ( IsGoodCandidate( arg.readBatch[i].seq, arg.buffer, arg.refSet ) )
			++goodCandidate ;
		
		if ( !goodCandidate && arg.readBatch2 && IsGoodCandidate( arg.readBatch2[i].seq, arg.buffer, arg.refSet ) )
			++goodCandidate ;

		if ( !goodCandidate ) 
			arg.readBatch[i].id[0] = '\0' ;
		/*else if (arg.barcodeBatch != NULL)
		{
			// Process the barcode. 
			char *buffer = arg.buffer ;
			char *barcode = arg.barcodeBatch[i].seq ;
			//printf("%d: %s %s. %s %s\n", goodCandidate, arg.readBatch[i].id, arg.readBatch[i].seq, 
			//		arg.barcodeBatch[i].id, arg.barcodeBatch[i].seq) ;
			if ( arg.barcodeStart == 0 && arg.barcodeEnd == -1 && arg.barcodeRevComp == false )
			{
				if ( arg.barcodeCorrector != NULL )
					strcpy( buffer, barcode ) ;
				else
					continue ; // no need to process the barcode
			}
			else
			{
				int s = arg.barcodeStart ;
				int e = ( arg.barcodeEnd == -1 ? strlen( barcode ) - 1 : arg.barcodeEnd ) ;
				if ( arg.barcodeRevComp == false )
				{
					for ( j = s ; j <= e ; ++j )
						buffer[j - s] = barcode[j] ;
					buffer[j - s] = '\0' ;
				}
				else
					arg.refSet->ReverseComplement( buffer, barcode + s, e - s + 1 ) ;
			}
				
			if ( arg.barcodeCorrector != NULL )
			{
				int result = 0 ;
				result = arg.barcodeCorrector->Correct( buffer, arg.barcodeBatch[i].qual ) ;
				if (result < 0)
					buffer[0] = '\0' ;
			}
		
			strcpy( barcode, buffer ) ;	
		}*/
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
	ReadFiles barcodeFile ;
	ReadFiles umiFile ;
	BarcodeCorrector barcodeCorrector ;
	BarcodeTranslator barcodeTranslator ;
	bool hasMate = false ;
	bool hasBarcode = false ;
	bool hasBarcodeWhitelist = false ;
	bool hasUmi = false ;
	int barcodeStart = 0 ;
	int barcodeEnd = -1 ;
	int read1Start = 0 ;
	int read1End = -1 ;
	int read2Start = 0 ;
	int read2End = -1 ;
	bool barcodeRevComp = false ;
	int umiStart = 0 ;
	int umiEnd = -1 ;
	bool umiRevComp = false ;
	ReadFormatter readFormatter ;

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
		else if ( c == 10000 ) // barcode
		{
			hasBarcode = true ;
			barcodeFile.AddReadFile( optarg, false ) ;
		}
		else if ( c == 10001 ) // barcodeStart
		{	
			barcodeStart = atoi( optarg ) ;
		}
		else if ( c == 10002 ) // barcodeEnd
		{
			barcodeEnd = atoi( optarg ) ;
		}
		else if ( c == 10003 ) // barcodeRevComp
		{
			barcodeRevComp = true ;
		}
		else if ( c == 10004 ) // barcodeWhitelist
		{
			hasBarcodeWhitelist = true ;
			barcodeCorrector.SetWhitelist( optarg ) ;
		}
		else if ( c == 10005 ) // read1Start
		{
			read1Start = atoi( optarg ) ;
		}
		else if ( c == 10006 ) // read1Start
		{
			read1End = atoi( optarg ) ;
		}
		else if ( c == 10007 ) // read1Start
		{
			read2Start = atoi( optarg ) ;
		}
		else if ( c == 10008 ) // read1Start
		{
			read2End = atoi( optarg ) ;
		}
		else if ( c == 10009 ) // umi
		{
			hasUmi = true ;
			umiFile.AddReadFile( optarg, false ) ;
		}
		else if ( c == 10010 ) // umiStart
		{	
			umiStart = atoi( optarg ) ;
		}
		else if ( c == 10011 ) // umiEnd
		{
			umiEnd = atoi( optarg ) ;
		}
		else if ( c == 10012 ) // umiRevComp
		{
			umiRevComp = true ;
		}
		else if (c == 10013) // readFormat
		{
			readFormatter.Init(optarg) ;
		}
		else if (c == 10014) // barcodeTranslate
		{
			barcodeTranslator.SetTranslateTable(optarg) ;
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

	if (read1Start != 0 || read1End != -1)
		readFormatter.AddSegment(read1Start, read1End, 1, FORMAT_READ1) ;
	if (read2Start != 0 || read2End != -1)
		readFormatter.AddSegment(read1Start, read1End, 1, FORMAT_READ2) ;
	if (barcodeStart != 0 || barcodeEnd != -1 || barcodeRevComp)
		readFormatter.AddSegment(barcodeStart, barcodeEnd, barcodeRevComp ? -1 : 1, FORMAT_BARCODE) ;
	if (umiStart != 0 || umiEnd != -1 || umiRevComp)
		readFormatter.AddSegment(umiStart, umiEnd, umiRevComp ? -1 : 1, FORMAT_UMI) ;
	
	PrintLog( "Start to extract candidate reads from read files." ) ;
	
	int hitLenRequired = 27 ;
	if ( !hasMate )
		hitLenRequired = 23 ;
	int len = 0 ;
	for ( i = 0 ; i < 1000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		len += strlen( reads.seq ) ;
	}
	if ( i == 0 )
	{
		fprintf( stderr, "Read file is empty.\n" ) ;
		return EXIT_FAILURE ;
	}
	if ( len / (i * 5) > hitLenRequired )
		hitLenRequired = len / (i * 5) ;
	refSet.SetHitLenRequired( hitLenRequired ) ;
	reads.Rewind() ;
	
	if ( hasBarcode && hasBarcodeWhitelist )
	{
		barcodeCorrector.CollectBackgroundDistribution(barcodeFile, readFormatter) ;
	}

	FILE *fp1 = NULL ;
	FILE *fp2 = NULL ;
	FILE *fpBc = NULL ;
	FILE *fpUmi = NULL ;
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

	if ( hasBarcode )
	{
		sprintf( buffer, "%s_bc.fa", prefix ) ;	
		fpBc = fopen( buffer, "w" ) ;
	}

	if ( hasUmi )
	{
		sprintf( buffer, "%s_umi.fa", prefix ) ;	
		fpUmi = fopen( buffer, "w" ) ;
	}
	
	if ( threadCnt == 1 )
	{
		while ( reads.Next() )
		{
			if ( hasMate && !mateReads.Next() )
			{
				fprintf( stderr, "The two mate-pair read files have different number of reads.\n" ) ;
				exit( 1 ) ;
			}

			if ( hasBarcode && !barcodeFile.Next() )
			{
				fprintf( stderr, "Read file and barcode file  have different number of reads.\n" ) ;
				exit( 1 ) ;
			}

			if ( hasUmi && !umiFile.Next() )
			{
				fprintf( stderr, "Read file and UMI file have different number of reads.\n" ) ;
				exit( 1 ) ;
			}

			int goodCandidate = 0 ;
			
			if ( IsGoodCandidate( reads.seq, seqBuffer, &refSet ) )
				++goodCandidate ;
			if ( !goodCandidate && hasMate && IsGoodCandidate( mateReads.seq, seqBuffer, &refSet ) )
				++goodCandidate ;
			if ( goodCandidate )
			{
				OutputSeq( fp1, reads.id, reads.seq, reads.qual, readFormatter, FORMAT_READ1 ) ;
				if ( hasMate )
					OutputSeq( fp2, reads.id, mateReads.seq, mateReads.qual, readFormatter, FORMAT_READ2 ) ;
				if ( hasBarcode )
					OutputBarcode( fpBc, reads.id, barcodeFile.seq, barcodeFile.qual, 
						readFormatter, FORMAT_BARCODE, 
						hasBarcodeWhitelist ? &barcodeCorrector : NULL, 
						barcodeTranslator.isSet() ? &barcodeTranslator : NULL) ;
				if ( hasUmi )
					OutputBarcode( fpUmi, reads.id, umiFile.seq, umiFile.qual,
							readFormatter, FORMAT_UMI, NULL, NULL) ;
			}
		}
	}
	else
	{
		int maxBatchSize = 512 * threadCnt ;
		int batchSize ;
		
		struct _Read *readBatch = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
		struct _Read *readBatch2 = NULL ;
		struct _Read *barcodeBatch = NULL ;
		struct _Read *umiBatch = NULL ;
		if ( hasMate )
			readBatch2 = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
		if ( hasBarcode )
			barcodeBatch = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
		if ( hasUmi )
			umiBatch = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
		int fileInd1, fileInd2, fileIndBc, fileIndUmi ;
		
		pthread_t *threads = (pthread_t *)malloc( sizeof( pthread_t ) * threadCnt ) ;
		struct _threadArg *args = (struct _threadArg *)malloc( sizeof( struct _threadArg ) * threadCnt ) ;
		pthread_attr_t attr ;
		pthread_attr_init( &attr ) ;
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;

		for ( i = 0 ; i < threadCnt ; ++i )
		{
			args[i].threadCnt = threadCnt ;
			args[i].tid = i ;
			args[i].readBatch = readBatch ;
			args[i].readBatch2 = readBatch2 ;
			args[i].barcodeBatch = barcodeBatch ;
			args[i].barcodeStart = barcodeStart ;
			args[i].barcodeEnd = barcodeEnd ;
			args[i].barcodeRevComp = barcodeRevComp ;
			if ( hasBarcodeWhitelist )
				args[i].barcodeCorrector = &barcodeCorrector ;
			else
				args[i].barcodeCorrector = NULL ;
			args[i].buffer = (char *)malloc( sizeof( char ) * 100001 ) ;
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
			
			if ( hasBarcode )
			{
				int tmp = barcodeFile.GetBatch( barcodeBatch, maxBatchSize, fileIndBc, true, true ) ;
				if ( tmp != batchSize )
				{
					fprintf( stderr, "Read file and barcode file have different number of reads.\n" ) ;
					exit ( 1 ) ;
				}
			}
			
			if ( hasUmi )
			{
				int tmp = umiFile.GetBatch( umiBatch, maxBatchSize, fileIndUmi, true, true ) ;
				if ( tmp != batchSize )
				{
					fprintf( stderr, "Read file and umi file have different number of reads.\n" ) ;
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
				pthread_join( threads[i], NULL ) ;

			for ( i = 0 ; i < batchSize ; ++i )
			{
				if ( readBatch[i].id[0] == '\0' )
					continue ;
				OutputSeq( fp1, readBatch[i].id, readBatch[i].seq, readBatch[i].qual, readFormatter, FORMAT_READ1 ) ;
				if ( readBatch2 != NULL )
					OutputSeq( fp2, readBatch[i].id, readBatch2[i].seq, readBatch2[i].qual, readFormatter, FORMAT_READ2) ;
				if ( hasBarcode )
					OutputBarcode( fpBc, readBatch[i].id, barcodeBatch[i].seq, barcodeBatch[i].qual, 
						readFormatter, FORMAT_BARCODE,
						hasBarcodeWhitelist ? &barcodeCorrector : NULL, 
						barcodeTranslator.isSet() ? &barcodeTranslator : NULL) ;
				
				if ( hasUmi )
					OutputBarcode( fpUmi, readBatch[i].id, umiBatch[i].seq, umiBatch[i].qual, 
						readFormatter, FORMAT_UMI, NULL, NULL) ; 
			}
		}
		
		// Release memory.
		for ( i = 0 ; i < maxBatchSize ; ++i )
		{
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

			if ( barcodeBatch != NULL )
			{
				free( barcodeBatch[i].id ) ;
				free( barcodeBatch[i].seq ) ;
				if ( barcodeBatch[i].qual )
					free( barcodeBatch[i].qual ) ;
			}
			
			if ( umiBatch != NULL )
			{
				free( umiBatch[i].id ) ;
				free( umiBatch[i].seq ) ;
				if ( umiBatch[i].qual )
					free( umiBatch[i].qual ) ;
			}
		}
		reads.id = NULL ;
		reads.seq = NULL ;
		reads.qual = NULL ;
		if ( hasMate )
		{
			mateReads.id = NULL ;
			mateReads.seq = NULL ;
			mateReads.qual = NULL ;
		}
		if ( hasBarcode )
		{
			barcodeFile.id = NULL ;
			barcodeFile.seq = NULL ;
			barcodeFile.qual = NULL ;
		}
		if ( hasUmi )
		{
			umiFile.id = NULL ;
			umiFile.seq = NULL ;
			umiFile.qual = NULL ;
		}

		free( threads ) ;
		for ( i = 0 ; i < threadCnt ; ++i )
			free( args[i].buffer ) ;
		free( args ) ;
		free( readBatch ) ;
		if ( hasMate )
			free( readBatch2 ) ;
		if ( hasBarcode )
			free( barcodeBatch ) ;
		pthread_attr_destroy( &attr ) ;
	}
	fclose( fp1 ) ;
	if ( hasMate )
		fclose( fp2 ) ;
	if ( hasBarcode )
		fclose( fpBc ) ;
	fclose( fpRef ) ; 
	PrintLog( "Finish extracting reads." ) ;
	return 0 ;
}
