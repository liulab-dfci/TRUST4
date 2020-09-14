#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <pthread.h>

#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include "alignments.hpp"
#include "SeqSet.hpp"

char usage[] = "./bam-extractor [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t-b STRING: path to BAM file\n"
		"Optional:\n"
		"\t-o STRING: prefix to the output file\n"
		"\t-t INT: number of threads (default: 1)\n"
		"\t-u: the flag or order of unaligned read-pair is not ordinary (default: not used)\n" 
		"\t--barcode STRING: the barcode field in the bam file (default: not used)\n"
		"\t--UMI STRING: the UMI field in the bam file (default: not used)\n"
		"\t--mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)\n";

static const char *short_options = "f:b:o:t:u" ;
static struct option long_options[] = {
			{ "barcode", required_argument, 0, 10000 },
			{ "UMI", required_argument, 0, 10001 },
			{ "mateIdSuffixLen", required_argument, 0, 10002},
			{ (char *)0, 0, 0, 0} 
			} ;

char buffer[100001] ;
char buffer2[100001] ;
char bufferQual[100001] ;
char bufferQual2[100001] ;
char seqBuffer[100001] ;
char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;


struct _interval
{
	int chrId ;
	int start, end ;

	bool operator<( const struct _interval &b ) const
	{
		if ( chrId != b.chrId )
			return chrId < b.chrId ;
		else if ( start != b.start )
			return start < b.start ;
		else 
			return end < b.end ;
	}
} ;

struct _candidate
{
	char *mate1 ;
	char *mate2 ;

	char *qual1 ;
	char *qual2 ;
} ;

struct _unmappedCandidate
{
	char *name ;

	char *mate1 ;
	char *mate2 ;

	char *qual1 ;
	char *qual2 ;
	
	char *barcode ;
	char *UMI ;
} ; 

// Pre-determined information
struct _threadInfo
{
	int tid ;
	
	FILE *fp1, *fp2 ;
	FILE *fpBc, *fpUMI ;
	SeqSet *refSet ;
	char *seqBuffer ;
	
	struct _unmappedCandidate *outputQueue ;
	int *oqCnt ;

	int *freeThreads ;
	int *ftCnt ;
	bool *initThreads ;
	

	pthread_mutex_t *lockOutput ;
	pthread_mutex_t *lockFreeThreads ;
	pthread_cond_t *condFreeThreads ;
} ;

struct _threadArg
{
	std::vector<struct _unmappedCandidate> candidates ;
	struct _threadInfo info ;
} ;


bool ValidAlternativeChrom( char *chrom )
{
	//if ( strstr( chrom, "_random" ) || strstr( chrom, "_alt" ) )
	if ( strstr( chrom, "_" ) || strstr( chrom, "." ) )
	{
		//if ( chrom[0] == 'c' && chrom[3] >= '0' && chrom[3] <= '9' )
		return true ;
	}
	//else if ( ( chrom[0] == 'G' && chrom[1] == 'I' ) || ( chrom[0] == 'G' && chrom[1] == 'L' ) )
	//	return true ;
	return false ;
}

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

void TrimName( std::string &name, int trimLen )
{
	int len = name.length() ;
	if (trimLen == -1)
	{
		if ( ( name[len - 1] == '1' || name[len - 1] == '2' ) 
				&& name[len - 2] == '/' )
		{
			name.erase( len - 2, 2 ) ; // Haven't tested yet.
		}
	}
	else
	{
		name.erase(len - trimLen, trimLen) ;
	}
}

void OutputSeq( FILE *fp, const char *name, char *seq, char *qual )
{
	if ( qual != NULL )
		fprintf( fp, "@%s\n%s\n+\n%s\n", name, seq, qual ) ;
	else
		fprintf( fp, ">%s\n%s\n", name, seq ) ;
}

// Maybe barcode read quality could be useful in future.
void OutputBarcode( FILE *fp, const char *name, char *barcode, char *qual )
{
	if ( barcode )
		fprintf( fp, ">%s\n%s\n", name, barcode ) ;
	else
		fprintf( fp, ">%s\nmissing_barcode\n", name ) ;
}

void *ProcessUnmappedReads_Thread( void *pArg )
{
	int i ;
	struct _threadArg &arg = *( (struct _threadArg *)pArg ) ;
	struct _threadInfo &info = arg.info ;
	int cnt = arg.candidates.size() ;
	SimpleVector<int> pick ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		if ( arg.candidates[i].mate2 == NULL )
		{
			// single-end 
			if ( !IsLowComplexity( arg.candidates[i].mate1 ) 
				&& info.refSet->HasHitInSet( arg.candidates[i].mate1, info.seqBuffer ) )
			{
				pick.PushBack( i ) ;
			}
		}
		else
		{
			// paired-end
			if (  ( !IsLowComplexity( arg.candidates[i].mate1 ) && !IsLowComplexity( arg.candidates[i].mate2 ) ) 
				&& ( info.refSet->HasHitInSet( arg.candidates[i].mate1, info.seqBuffer ) 
					|| info.refSet->HasHitInSet( arg.candidates[i].mate2, info.seqBuffer ) ) )
			{
				pick.PushBack( i ) ;
			}
		}
	}
	
	// Output
	cnt = pick.Size() ;
	pthread_mutex_lock( info.lockOutput ) ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		int k = pick[i] ;
		info.outputQueue[ *info.oqCnt ] = arg.candidates[k] ;
		++*(info.oqCnt) ;
		arg.candidates[k].name = NULL ;
	}

	if ( *info.oqCnt > 2 * arg.candidates.size() )
	{
		for ( i = 0 ; i < *info.oqCnt ; ++i )
		{
			OutputSeq( info.fp1, info.outputQueue[i].name, info.outputQueue[i].mate1, info.outputQueue[i].qual1 ) ;
			free( info.outputQueue[i].mate1 ) ;
			free( info.outputQueue[i].qual1 ) ;
			if ( info.outputQueue[i].mate2 != NULL )
			{
				OutputSeq( info.fp2, info.outputQueue[i].name, info.outputQueue[i].mate2, info.outputQueue[i].qual2 ) ;
				free( info.outputQueue[i].mate2 ) ;
				free( info.outputQueue[i].qual2 ) ;
			}
			if ( info.fpBc )
			{
				OutputBarcode( info.fpBc, info.outputQueue[i].name, info.outputQueue[i].barcode, NULL ) ;
				if ( info.outputQueue[i].barcode )
					free( info.outputQueue[i].barcode ) ;
			}

			if ( info.fpUMI )
			{
				OutputBarcode( info.fpUMI, info.outputQueue[i].name, info.outputQueue[i].UMI, NULL ) ;
				if ( info.outputQueue[i].UMI )
					free( info.outputQueue[i].UMI ) ;
			}
			free( info.outputQueue[i].name ) ;
		}
		*info.oqCnt = 0 ;
	}
	pthread_mutex_unlock( info.lockOutput ) ;
	
	// Release memory
	cnt = arg.candidates.size() ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		if ( arg.candidates[i].name == NULL )
			continue ;

		free( arg.candidates[i].name ) ;
		free( arg.candidates[i].mate1 ) ;
		free( arg.candidates[i].qual1 ) ;
		if ( arg.candidates[i].mate2 != NULL )
		{
			free( arg.candidates[i].mate2 ) ;
			free( arg.candidates[i].qual2 ) ;
		}

		if ( arg.candidates[i].barcode != NULL )
			free( arg.candidates[i].barcode ) ;
		if ( arg.candidates[i].UMI != NULL )
			free( arg.candidates[i].UMI ) ;
	}

	pthread_mutex_lock( info.lockFreeThreads ) ;
	info.freeThreads[ *( info.ftCnt ) ] = info.tid ;
	++*(info.ftCnt) ;
	if ( *info.ftCnt == 1 )
		pthread_cond_signal( info.condFreeThreads ) ;
	pthread_mutex_unlock( info.lockFreeThreads ) ;
	pthread_exit( NULL ) ;
	return NULL ;
}

// Framework of work distribute model
void InitWork( pthread_t **threads, struct _threadArg **threadArgs, pthread_attr_t &attr, int threadCnt )
{
	int i ;
	struct _threadInfo info ;
	
	*threads = new pthread_t[ threadCnt ] ;
	*threadArgs = new struct _threadArg[ threadCnt ] ;
	pthread_attr_init( &attr ) ;
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;

	int *freeThreads = new int[ threadCnt ] ;
	int *ftCnt = new int ;
	bool *initThreads = new bool[threadCnt] ;
	*ftCnt = threadCnt ;
	pthread_mutex_t *lockOutput = new pthread_mutex_t ;
	pthread_mutex_t *lockFreeThreads = new pthread_mutex_t ;
	pthread_cond_t *condFreeThreads = new pthread_cond_t ;
	
	pthread_mutex_init( lockOutput, NULL ) ;
	pthread_mutex_init( lockFreeThreads, NULL ) ;
	pthread_cond_init( condFreeThreads, NULL ) ;
	for ( i = 0 ; i < threadCnt ; ++i )
	{
		struct _threadInfo &info = (*threadArgs)[i].info ;
		info.tid = i ;
		info.freeThreads = freeThreads ;
		info.initThreads = initThreads ;
		info.ftCnt = ftCnt ;
		info.lockOutput = lockOutput ;
		info.lockFreeThreads = lockFreeThreads ;
		info.condFreeThreads = condFreeThreads ;

		freeThreads[i] = i ;
		initThreads[i] = false ;
	}
}

// Customize your own parameters here.
void InitCustomData( FILE *fp1, FILE *fp2, FILE *fpBc, FILE *fpUMI, SeqSet *refSet, struct _threadArg *threadArgs, int threadCnt )
{
	int i ; 
	struct _unmappedCandidate *outputQueue = new struct _unmappedCandidate[2048 * 5] ;
	int *oqCnt = new int ;
	*oqCnt = 0 ;
	for ( i = 0 ; i < threadCnt ; ++i )
	{
		threadArgs[i].info.fp1 = fp1 ;			
		threadArgs[i].info.fp2 = fp2 ;	
		threadArgs[i].info.fpBc = fpBc ;
		threadArgs[i].info.fpUMI = fpUMI ;
		threadArgs[i].info.refSet = refSet ;
		threadArgs[i].info.seqBuffer = new char[100001] ;
		threadArgs[i].info.outputQueue = outputQueue ;
		threadArgs[i].info.oqCnt = oqCnt ;
	}
}

void DistributeWork( std::vector<struct _unmappedCandidate> &work,
	struct _threadArg *threadArgs, pthread_t *threads, pthread_attr_t &attr, int threadCnt )
{
	int tid ;
	struct _threadInfo &info = threadArgs[0].info ;

	// Determine which thread to use
	//printf( "Threads status %d/%d.\n", *(info.ftCnt), threadCnt ) ;
	pthread_mutex_lock( info.lockFreeThreads ) ;
	if ( *(info.ftCnt) == 0 )
		pthread_cond_wait( info.condFreeThreads, info.lockFreeThreads ) ;
	tid = info.freeThreads[ *( info.ftCnt ) - 1 ] ;
	--*( info.ftCnt ) ;
	pthread_mutex_unlock( info.lockFreeThreads ) ;
	
	// Make sure the thread ends.
	if ( info.initThreads[tid] )
		pthread_join( threads[ tid ], NULL ) ;

	// Call the thread
	threadArgs[tid].candidates = work ;
	work.clear() ;
	info.initThreads[tid] = true ;
	//printf( "%d %lld %lld %lld\n", tid, threads, &attr, threadArgs + tid ) ;
	pthread_create( &threads[tid], &attr, ProcessUnmappedReads_Thread, (void *)( threadArgs + tid ) ) ;
}

void AddWorkQueue( struct _unmappedCandidate &c, std::vector<struct _unmappedCandidate> &work, 
	struct _threadArg *threadArgs, pthread_t *threads, pthread_attr_t &attr, int workLoad, int threadCnt )
{
	work.push_back( c ) ;
	if ( work.size() >= workLoad )
		DistributeWork( work, threadArgs, threads, attr, threadCnt ) ; 
}

void ReleaseCustomData( struct _threadArg *threadArgs, int threadCnt )
{
	int i ;
	for ( i = 0 ; i < threadCnt ; ++i )
		delete[] threadArgs[i].info.seqBuffer ;
	delete[] threadArgs[0].info.outputQueue ;
	delete threadArgs[0].info.oqCnt ;
}

void FinishWork( std::vector<struct _unmappedCandidate> work,  
	struct _threadArg *threadArgs, pthread_t *threads, pthread_attr_t &attr, int threadCnt )
{	
	int i ;
	DistributeWork( work, threadArgs, threads, attr, threadCnt ) ; 
	struct _threadInfo &info = threadArgs[0].info ;

	for ( i = 0 ; i < threadCnt ; ++i )
	{
		if ( info.initThreads[i] )
			pthread_join( threads[i], NULL ) ;
	}

	for ( i = 0 ; i < *info.oqCnt ; ++i )
	{
		OutputSeq( info.fp1, info.outputQueue[i].name, info.outputQueue[i].mate1, info.outputQueue[i].qual1 ) ;
		free( info.outputQueue[i].mate1 ) ;
		free( info.outputQueue[i].qual1 ) ;
		if ( info.outputQueue[i].mate2 != NULL )
		{
			OutputSeq( info.fp2, info.outputQueue[i].name, info.outputQueue[i].mate2, info.outputQueue[i].qual2 ) ;
			free( info.outputQueue[i].mate2 ) ;
			free( info.outputQueue[i].qual2 ) ;
		}
		if ( info.fpBc )
		{
			OutputBarcode( info.fpBc, info.outputQueue[i].name, info.outputQueue[i].barcode, NULL ) ;
			if ( info.outputQueue[i].barcode )
				free( info.outputQueue[i].barcode ) ;
		}
		if ( info.fpUMI )
		{
			OutputBarcode( info.fpUMI, info.outputQueue[i].name, info.outputQueue[i].UMI, NULL ) ;
			if ( info.outputQueue[i].UMI )
				free( info.outputQueue[i].UMI ) ;
		}
		free( info.outputQueue[i].name ) ;
	}
	// Release memory 
	delete[] threads ;
	pthread_attr_destroy( &attr ) ;
	delete[] info.freeThreads ;
	delete info.ftCnt ;
	delete[] info.initThreads ;
	pthread_mutex_destroy( info.lockOutput ) ;
	pthread_mutex_destroy( info.lockFreeThreads ) ;
	pthread_cond_destroy( info.condFreeThreads ) ;
	delete info.lockOutput ;
	delete info.lockFreeThreads ;
	delete info.condFreeThreads ;
	ReleaseCustomData( threadArgs, threadCnt ) ;
	delete[] threadArgs ;
}


int main( int argc, char *argv[] )
{
	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ;
	option_index = 0 ;
	FILE *fpRef = NULL ;
	char prefix[1024] = "toassemble" ;
	char bcField[1024] = "" ; // barcode field in bam file.
	char umiField[1024] = "" ; // UMI field
	Alignments alignments ;
	bool abnormalUnalignedFlag = false ;
	int kmerLength = 9  ;
	SeqSet refSet( kmerLength ) ;
	int threadCnt = 1 ;
	int mateIdLen = -1 ;

	std::map<std::string, struct _candidate> candidates ; 
	std::map<std::string, int> usedName ; // For single-end case.

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
		else if ( c == 'b' )
		{
			alignments.Open( optarg ) ;
		}
		else if ( c == 'o' )
		{
			strcpy( prefix, optarg ) ;
		}
		else if ( c == 'u' )
		{
			abnormalUnalignedFlag = true ;
		}
		else if ( c == 't' )
		{
			threadCnt = atoi( optarg ) ;
		}
		else if ( c == 10000 ) // barcode
		{
			strcpy( bcField, optarg ) ;
		}
		else if ( c == 10001 ) // UMI
		{
			strcpy( umiField, optarg ) ;
		}
		else if ( c == 10002 ) // mateIdSuffixLen
		{
			mateIdLen = atoi( optarg ) ;
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

	if ( !alignments.IsOpened() )
	{
		fprintf( stderr, "Need to use -b to specify the BAM file.\n" ) ;
		return EXIT_FAILURE ;
	}

	// Scan the header to store the coordinates of the intervals
	char geneName[101] ;
	char chrom[101] ;
	int start, end ;
	char strand[3] ;

	std::vector<struct _interval> genes ;
	PrintLog( "Start to extract candidate reads from bam file." ) ;

	while ( fscanf( fpRef, "%s %s %d %d %s", geneName, chrom, &start, &end, strand ) != EOF )
	{
		struct _interval ni ;
		ni.chrId = alignments.GetChromIdFromName( chrom ) ;
		ni.start = start ;
		ni.end = end ;

		genes.push_back( ni ) ;

		fscanf( fpRef, "%s", buffer ) ;
	}

	int geneCnt = genes.size() ;
	std::sort( genes.begin(), genes.end() ) ;

	// Go through the BAM file to output the read id.
	alignments.GetGeneralInfo( true ) ;
	alignments.Rewind() ;
	
	int hitLenRequired = 21 ;
	if ( alignments.fragStdev == 0 )
		hitLenRequired = 17 ; // For single end, be more aggressive.
	if ( alignments.readLen / 5 > hitLenRequired )
		hitLenRequired = alignments.readLen / 5 ;
	//refSet.SetRadius( 1 ) ;	
	refSet.SetHitLenRequired( hitLenRequired ) ;

	int tag = 0 ;
	
	FILE *fp1 = NULL ;
	FILE *fp2 = NULL ;
	FILE *fpBc = NULL ; // the file holding barcode.
	FILE *fpUMI = NULL ; // the file holding UMI
	if ( alignments.fragStdev == 0 )
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
	if ( bcField[0] )
	{
		sprintf( buffer, "%s_bc.fa", prefix ) ;
		fpBc = fopen( buffer, "w" ) ;
	}
	if ( umiField[0] )
	{
		sprintf( buffer, "%s_umi.fa", prefix ) ;
		fpUMI = fopen( buffer, "w" ) ;
	}

	pthread_t *threads ;
	struct _threadArg *threadArgs ;
	pthread_attr_t attr ;
	std::vector<struct _unmappedCandidate> threadsWorkQueue ;
	if ( threadCnt > 1 )
	{
		InitWork( &threads, &threadArgs, attr, threadCnt - 1 ) ;
		InitCustomData( fp1, fp2, fpBc, fpUMI, &refSet, threadArgs, threadCnt - 1 ) ;
	}
	// assuming the input is sorted by coordinate.
	while ( alignments.Next() )
	{
		// If not aligned, output it.
		//if ( !alignments.IsPrimary() )
		//	continue ;
		//if ( IsLowComplexity( buffer ) )
		//	continue ;

		if ( !alignments.IsTemplateAligned() 
			|| ( alignments.IsAligned() && ValidAlternativeChrom( alignments.GetChromName( alignments.GetChromId() ) ) ) )
		{
			//if ( !alignments.IsTemplateAligned() && !refSet.HasHitInSet( buffer ) ) 
			//	continue ;
			
			if ( !alignments.IsTemplateAligned() && alignments.fragStdev != 0 
				&& !abnormalUnalignedFlag ) // And the two reads of unaligned template should come together. 
			{
				//printf( "filtered\n" ) ;
				alignments.GetReadSeq( buffer ) ;
				alignments.GetQual( bufferQual ) ;
				
				std::string name( alignments.GetReadId() ) ;
				strcpy( buffer2, buffer ) ;
				alignments.GetQual( bufferQual2 ) ;

				if ( !alignments.Next() )
				{
					fprintf( stderr, "Two reads from the unaligned fragment are not showing up together. Please use -u(--abnormalUnmapFlag from wrapper) option.\n") ;
					return EXIT_FAILURE ;
				}
				std::string mateName( alignments.GetReadId() ) ;
				alignments.GetReadSeq( buffer ) ;
				alignments.GetQual( bufferQual ) ;
				TrimName( name, mateIdLen ) ;
				TrimName( mateName, mateIdLen ) ;
				if ( name.compare( mateName ) != 0 )
				{
					fprintf( stderr, "%s\n%s\n", name.c_str(), mateName.c_str() ) ;
					fprintf( stderr, "Two reads from the unaligned fragment are not showing up together. Please use -u(--abnormalUnmapFlag from wrapper) option.\n") ;
					return EXIT_FAILURE ;
				}


				if ( threadCnt == 1 )
				{
					if ( ( !IsLowComplexity( buffer2 ) && !IsLowComplexity( buffer ) ) && 
						( refSet.HasHitInSet( buffer2, seqBuffer ) || 
						 refSet.HasHitInSet( buffer, seqBuffer ) ) ) 
					{
						if ( !alignments.IsFirstMate() )
						{
							OutputSeq( fp1, name.c_str(), buffer2, bufferQual2 ) ;
							OutputSeq( fp2, name.c_str(), buffer, bufferQual  ) ;
						}
						else
						{
							OutputSeq( fp1, name.c_str(), buffer, bufferQual ) ;
							OutputSeq( fp2, name.c_str(), buffer2, bufferQual2 ) ;
						}
						if ( fpBc != NULL )
							OutputBarcode( fpBc, name.c_str(), alignments.GetFieldZ( bcField ), NULL ) ;
						if ( fpUMI != NULL )
							OutputBarcode( fpUMI, name.c_str(), alignments.GetFieldZ( umiField ), NULL ) ;
					}
				}
				else
				{
					struct _unmappedCandidate nw ;
					nw.name = strdup( name.c_str() ) ;
					if ( !alignments.IsFirstMate() )
					{
						nw.mate1 = strdup( buffer2 ) ;
						nw.qual1 = strdup( bufferQual2 ) ;
						nw.mate2 = strdup( buffer ) ;
						nw.qual2 = strdup( bufferQual ) ;
					}
					else
					{
						nw.mate1 = strdup( buffer ) ;
						nw.qual1 = strdup( bufferQual ) ;
						nw.mate2 = strdup( buffer2 ) ;
						nw.qual2 = strdup( bufferQual2 ) ;
					}

					if ( bcField[0] && alignments.GetFieldZ( bcField ) )
						nw.barcode = strdup( alignments.GetFieldZ( bcField ) ) ;
					else
						nw.barcode = NULL ;
					
					if ( umiField[0] && alignments.GetFieldZ( umiField ) )
						nw.UMI = strdup( alignments.GetFieldZ( umiField ) ) ;
					else
						nw.UMI = NULL ;

					AddWorkQueue( nw, threadsWorkQueue, threadArgs, threads, attr, 2048, threadCnt - 1 ) ;
				}
				continue ;
			}
			
			//printf( "%s %s\n", alignments.GetChromName( alignments.GetChromId() ), alignments.GetReadId() ) ;
			if ( alignments.fragStdev != 0 )
			{
				// reads from alternative chromosomes, or the unmapped flag is not set appropriately
				alignments.GetReadSeq( buffer ) ;
				alignments.GetQual( bufferQual ) ;
				if ( !IsLowComplexity( buffer ) && refSet.HasHitInSet( buffer, seqBuffer ) )
				{
					std::string name( alignments.GetReadId() ) ;
					TrimName( name, mateIdLen ) ;	

					if ( candidates.find( name ) == candidates.end() )
					{
						candidates[name].mate1 = NULL ;
						candidates[name].mate2 = NULL ;
					}
				}
			}
			else 
			{
				// single-end 
				alignments.GetReadSeq( buffer ) ;
				alignments.GetQual( bufferQual ) ;
				if ( threadCnt == 1 ||  alignments.IsAligned() ) 
				{
					// If a read is from alternative chromosome, it could be multiple aligned,
					// so we need to check and mark the read id usage.
					std::string name ;
					if ( alignments.IsAligned() )
					{
						name = std::string( alignments.GetReadId() ) ;
						if ( usedName.find( name ) != usedName.end() )
							continue ;
						
					}
					if ( !IsLowComplexity( buffer ) && refSet.HasHitInSet( buffer, seqBuffer ) )
					{
						//alignments.GetReadSeq( buffer ) ;
						// No need to trim read id for single-end data.
						if ( alignments.IsAligned() )
							usedName[ name ] = 1 ;
						OutputSeq( fp1, alignments.GetReadId(), buffer, bufferQual ) ;
						if ( fpBc != NULL )
							OutputBarcode( fpBc, alignments.GetReadId(), alignments.GetFieldZ( bcField ), NULL ) ;
						if ( fpUMI != NULL )
							OutputBarcode( fpUMI, alignments.GetReadId(), alignments.GetFieldZ( umiField ), NULL ) ;
					}
				}
				else
				{
					struct _unmappedCandidate nw ;
					nw.name = strdup( alignments.GetReadId() ) ;
					nw.mate1 = strdup( buffer ) ;
					nw.qual1 = strdup( bufferQual ) ;
					nw.mate2 = NULL ; 
					nw.qual2 = NULL ; 
					if ( bcField[0] && alignments.GetFieldZ( bcField ) )
						nw.barcode = strdup( alignments.GetFieldZ( bcField ) ) ;
					else
						nw.barcode = NULL ;
					if ( umiField[0] && alignments.GetFieldZ( umiField ) )
						nw.UMI = strdup( alignments.GetFieldZ( umiField ) ) ;
					else
						nw.UMI = NULL ;
					AddWorkQueue( nw, threadsWorkQueue, threadArgs, threads, attr, 2048, threadCnt - 1 ) ;
				}
			}
			continue ;
		}

		if ( !alignments.IsAligned() ) // when reach here, it is parie-end case, and the other mate is aligned.
			continue ;
		
		alignments.GetReadSeq( buffer ) ;
		if ( IsLowComplexity( buffer ) )
			continue ;

		// The aligned reads can reach here.
		int chrId = alignments.GetChromId() ;
		int start = (int)alignments.segments[0].a ;
		int end = (int)alignments.segments[ alignments.segCnt - 1 ].b ;
		//if ( !strcmp( "ERR188021.4488674", alignments.GetReadId() ) )
			//printf( "hi %d %d %d: %d %d\n", chrId, start, end, tag, genes[tag].chrId ) ;
		while ( tag < geneCnt && ( chrId > genes[tag].chrId ||
					( chrId == genes[tag].chrId && start > genes[tag].end ) ) )
			++tag ;
	
		if ( tag >= geneCnt )
			continue ;
		if ( chrId < genes[tag].chrId 
				|| ( chrId == genes[tag].chrId && end <= genes[tag].start ) )
			continue ;
		
		// overlaps with reference gene location.
		if ( alignments.fragStdev != 0 )
		{
			std::string name( alignments.GetReadId() ) ;
			TrimName( name, mateIdLen ) ;

			if ( candidates.find( name ) == candidates.end() )
			{
				candidates[name].mate1 = NULL ;
				candidates[name].mate2 = NULL ;
			}
		}
		else
		{
			std::string name( alignments.GetReadId() ) ;
			if ( usedName.find( name ) != usedName.end() )
				continue ;
			usedName[ name ] = 1 ;
			
			alignments.GetQual( bufferQual ) ;
			//alignments.GetReadSeq( buffer ) ;
			OutputSeq( fp1, alignments.GetReadId(), buffer, bufferQual ) ;
			if ( fpBc != NULL )
				OutputBarcode( fpBc, alignments.GetReadId(), alignments.GetFieldZ( bcField ), NULL ) ;
			if ( fpUMI != NULL )
				OutputBarcode( fpUMI, alignments.GetReadId(), alignments.GetFieldZ( umiField ), NULL ) ;
		}
	}

	if ( threadCnt > 1 )
	{
		FinishWork( threadsWorkQueue, threadArgs, threads, attr, threadCnt - 1 ) ;
	}

	alignments.Rewind() ;
	if ( alignments.fragStdev == 0 ) // Single-end can terminate here.
	{
		fclose( fp1 ) ;
		fclose( fpRef ) ;
		if ( fpBc != NULL )
			fclose( fpBc ) ;
		if ( fpUMI != NULL )
			fclose( fpUMI ) ;
		alignments.Close() ;
		PrintLog( "Finish extracting reads." ) ;
		return 0 ;
	}
	// Case of pair-end data set.
	// Go through the BAM file again to output the candidates 
	PrintLog( "Finish obtaining the candidate read ids." ) ;

	int candidateCnt = candidates.size() ;
	int outputCnt = 0 ;
	//std::map< std::string, int> nameCnt ;
	while ( alignments.Next() )
	{
		if ( !alignments.IsPrimary() )
			continue ;
		if ( !alignments.IsTemplateAligned() && !abnormalUnalignedFlag ) 
			// the sorted bam file should put all unaligned template at last (NO!).
			continue ;

		std::string name( alignments.GetReadId() ) ;
		TrimName( name, mateIdLen ) ;
		
		std::map<std::string, struct _candidate>::iterator it = candidates.find( name ) ;
		if ( it == candidates.end() )
			continue ;

		alignments.GetReadSeq( buffer ) ;
		alignments.GetQual( bufferQual ) ;
		if ( alignments.IsFirstMate() )
		{
			it->second.mate1 = strdup( buffer ) ;
			it->second.qual1 = strdup( bufferQual ) ;
		}
		else
		{
			it->second.mate2 = strdup( buffer ) ;
			it->second.qual2 = strdup( bufferQual ) ;	
		}
		
		/*if ( nameCnt.find( name ) == nameCnt.end() )
		{
			nameCnt[ name ] = 1 ;
		}
		else
		{
			++nameCnt[ name ] ;
			if ( nameCnt[name] >= 3 )
				printf( "Error!! %s\n", name.c_str() ) ;
		}*/

		if ( it->second.mate1 != NULL && it->second.mate2 != NULL )
		{
			OutputSeq( fp1, name.c_str(), it->second.mate1, it->second.qual1 ) ;
			OutputSeq( fp2, name.c_str(), it->second.mate2, it->second.qual2 ) ;
			if ( fpBc != NULL )
				OutputBarcode( fpBc, name.c_str(), alignments.GetFieldZ( bcField ), NULL ) ;
			if ( fpUMI != NULL )
				OutputBarcode( fpUMI, name.c_str(), alignments.GetFieldZ( umiField ), NULL ) ;
			free( it->second.mate1 ) ;
			free( it->second.mate2 ) ;
			free( it->second.qual1 ) ;
			free( it->second.qual2 ) ;

			it->second.mate1 = NULL ;
			it->second.mate2 = NULL ;
			
			++outputCnt ;
			if ( outputCnt == candidateCnt )
				break ;
		}
			
	}
	fclose( fp1 ) ;
	if ( fp2 != NULL )
		fclose( fp2 ) ;
	if ( fpBc != NULL )
		fclose( fpBc ) ;
	if ( fpUMI != NULL )
		fclose( fpUMI ) ;
	fclose( fpRef ) ;
	PrintLog( "Finish extracting reads." ) ;
	return 0 ;
}
