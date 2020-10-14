#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <pthread.h>

#include <vector>

#include "KmerCount.hpp"
#include "SeqSet.hpp"
#include "AlignAlgo.hpp"

char usage[] = "./trust4 [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t[Read file]\n"
		"\t-u STRING: path to single-end read file\n"
		"\t\tor\n"
		"\t-1 STRING -2 STRING: path to paried-end read files\n"
		"\t\tor\n"
		"\t-b STRING: path to BAM alignment file\n"
		"Optional:\n"
		"\t-o STRING: prefix of the output file (default: trust)\n"
		"\t-t INT: number of threads (default: 1)\n"
		"\t-c STRING: the path to the kmer count file\n"
		"\t--skipMateExtension: skip the step of extension assemblies with mate-pair information\n"
		///"\t--noV: do not assemble the full length V gene (default: not used)\n"
		"\t--trimLevel INT: 0: no trim; 1: trim low quality; 2: trim unmatched (default: 1)\n"
		"\t--barcode STRING: the path to the barcode file (default: not used)\n"
		"\t--UMI STRING: the path to the UMI file (default: not used)\n"
		"\t--keepNoBarcode: assemble the reads with missing barcodes. (default: ignore the reads)\n" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

char buffer[10240] = "" ;

static const char *short_options = "f:u:1:2:b:o:c:t:" ;
static struct option long_options[] = {
			{ "debug-ns", required_argument, 0, 10000 },
			{ "trimLevel", required_argument, 0, 10001 },
			{ "barcode", required_argument, 0, 10002 },
			{ "keepNoBarcode", no_argument, 0, 10003 },
			{ "UMI", required_argument, 0, 10004},
			{ "skipMateExtension", no_argument, 0, 10005},
			{ (char *)0, 0, 0, 0} 
			} ;

struct _sortRead
{
	char *id ;
	char *read ;
	char *qual ;
	int minCnt ;
	int medianCnt ;
	double avgCnt ;
	int len ;

	int strand ;

	int mateIdx ; // the index of its mate pair.
	int info ; // some random information, such as it's index in orginal array.

	int barcode ; // this is for cell id
	int umi ; // this is for mrna id.

	struct _overlap geneOverlap[4] ; // Rough annotation of the read

	bool operator<( const struct _sortRead &b ) const 
	{
		if ( minCnt != b.minCnt )
			return minCnt > b.minCnt ;
		else if ( medianCnt != b.medianCnt )
			return medianCnt > b.medianCnt ;
		else if ( avgCnt != b.avgCnt )
			return avgCnt > b.avgCnt ;
		else if ( len != b.len )
			return len > b.len ;
		else if (barcode != b.barcode)
			return barcode < b.barcode ;
		else 
			return strcmp( read, b.read ) < 0 ;
	}
} ;

struct _quickAnnotateReadsThreadArg
{
	int tid ;
	int threadCnt ; 

	SeqSet *refSet ;
	std::vector<struct _sortRead> *pSortedReads ;
	int readCnt ;
} ;

struct _assignReadsThreadArg
{
	int tid ;
	int threadCnt ;

	SeqSet *seqSet ;
	std::vector<struct _assignRead> *pAssembledReads ;
	int assembledReadCnt ;
} ;


bool CompSortReadById( const struct _Read &a, const struct _Read &b )
{
	return strcmp( a.id, b.id ) < 0 ;
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

void *QuickAnnotateReads_Thread( void *pArg )
{
	struct _quickAnnotateReadsThreadArg &arg = *( (struct _quickAnnotateReadsThreadArg *)pArg ) ;
	int start, end ;
	int i, j ;
	std::vector< struct _sortRead> &sortedReads = *arg.pSortedReads ;
	start = arg.readCnt / arg.threadCnt * arg.tid ;
	end = start + arg.readCnt / arg.threadCnt ;
	if ( arg.tid == arg.threadCnt - 1 )
		end = arg.readCnt ;
	struct _overlap geneOverlap[4] ;
	for ( i = start ; i < end ; ++i )
	{
		if ( i == start || strcmp( sortedReads[i].read, sortedReads[i - 1].read ) )
			arg.refSet->AnnotateRead( sortedReads[i].read, 0, geneOverlap, NULL, NULL ) ;
		for ( j = 0 ; j < 4 ; ++j )
			sortedReads[i].geneOverlap[j] = geneOverlap[j] ;
	}
	pthread_exit( NULL ) ;
}

void *AssignReads_Thread( void *pArg )
{
	struct _assignReadsThreadArg &arg = *( (struct _assignReadsThreadArg *)pArg ) ;
	int start, end ;
	int i ;
	std::vector< struct _assignRead> &assembledReads = *arg.pAssembledReads ;
	start = arg.assembledReadCnt / arg.threadCnt * arg.tid ;
	end = start + arg.assembledReadCnt / arg.threadCnt ;
	if ( arg.tid == arg.threadCnt - 1 )
		end = arg.assembledReadCnt ;
	struct _overlap assign ;
	for ( i = start ; i < end ; ++i )
	{
		if ( i == start || strcmp( assembledReads[i].read, assembledReads[i - 1].read ) )
			arg.seqSet->AssignRead( assembledReads[i].read, assembledReads[i].overlap.strand, 
				assembledReads[i].barcode, assign ) ;
		assembledReads[i].overlap = assign ;	
	}
	pthread_exit( NULL ) ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;

	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ;
	option_index = 0 ;
	int indexKmerLength = 9 ;
	int changeKmerLengthThreshold = 4096 ;
	SeqSet seqSet( indexKmerLength ) ; // Only hold the novel seq.
	SeqSet refSet( 9 ) ;
	KmerCount kmerCount( 21 ) ;
	char outputPrefix[1024] = "trust" ;

	ReadFiles reads ;
	ReadFiles mateReads ;
	ReadFiles barcodeFile, umiFile ;
	bool countMyself = true ;
	int maxReadLen = -1 ;
	int firstReadLen = -1 ;
	int trimLevel = 1 ;
	bool hasMate = false ;
	bool hasBarcode = false ;
	bool hasUmi = false ;
	int constantGeneEnd = 200 ;
	bool keepMissingBarcode = false ;
	int threadCnt = 1 ;
	bool skipMateExtension = false ;

	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			//seqSet.InputRefFa( optarg ) ;
			refSet.InputRefFa( optarg ) ;
		}
		else if ( c == 'u' )
		{
			reads.AddReadFile( optarg, false ) ;
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
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;
		}
		else if ( c == 'c' )
		{
			kmerCount.AddCountFromFile( optarg ) ;
			countMyself = false ;
			PrintLog( "Read in the kmer count information from %s", optarg ) ;
		}
		else if ( c == 't' )
		{
			threadCnt = atoi( optarg ) ;
		}
		else if ( c == 10000 ) //debug-ns
		{
			seqSet.InputNovelFa( optarg ) ;
		}
		else if ( c == 10001) //trimLevel
		{
			trimLevel = atoi( optarg ) ;
		}
		else if ( c == 10002 ) // barcode 
		{
			barcodeFile.AddReadFile( optarg, false ) ;
			hasBarcode = true ;
		}
		else if ( c == 10003 ) // keepNoBarcode
		{
			keepMissingBarcode = true ;
		}
		else if ( c == 10004 ) // umi
		{
			umiFile.AddReadFile( optarg, false ) ;
			hasUmi = true ;
		}
		else if ( c == 10005 ) // skipMateExtension
		{
			skipMateExtension = true ;
		}
		else
		{
			fprintf( stderr, "%s", usage ) ;
			return EXIT_FAILURE ;
		}
	}
	/*char align[10000] ;
	char t[] = "ACGGGGTTTTTT" ;
	char p[] = "ACTTTTTTGGGG" ;
	AlignAlgo::GlobalAlignment( t, strlen( t ), p, strlen( p ), align ) ;
	AlignAlgo::VisualizeAlignment( t, strlen( t ), p, strlen( p ), align ) ; */
	
	if ( refSet.Size() == 0 )
	{
		fprintf( stderr, "Need to use -f to specify the receptor genome sequence.\n" ) ;
		return EXIT_FAILURE ;
	}

	if ( trimLevel > 1 ) 
	{
		// We need a more careful rough annotation.
		refSet.ChangeKmerLength( 7 ) ;
	} 
	refSet.SetHitLenRequired( 17 ) ;

	std::vector< struct _sortRead > sortedReads ;
	std::map<std::string, int> barcodeStrToInt ;
	std::map<std::string, int> umiStrToInt ;
	std::vector<std::string> barcodeIntToStr ;

	i = 0 ;
	while ( reads.Next() )
	{
		/*struct _overlap geneOverlap[4] ;
		refSet.AnnotateRead( reads.seq, 1, geneOverlap, buffer ) ;
	
		if ( geneOverlap[0].seqIdx + geneOverlap[1].seqIdx + geneOverlap[2].seqIdx + geneOverlap[3].seqIdx == -4 )
		{
			continue ;
		}*/

		int barcode = -1 ;
		int umi = -1 ;
		if ( hasBarcode )
		{
			barcodeFile.Next() ;
			
			if ( !strcmp( barcodeFile.seq, "missing_barcode" ) && !keepMissingBarcode )
			{
				if ( hasMate )
					mateReads.Next() ;
				continue ;
			}

			std::string s( barcodeFile.seq ) ;
			if ( barcodeStrToInt.find( s ) != barcodeStrToInt.end() )
				barcode = barcodeStrToInt[s] ;
			else
			{
				barcode = barcodeIntToStr.size() ;
				barcodeStrToInt[s] = barcode ;
				barcodeIntToStr.push_back( s ) ;
			}
		}

		if ( hasUmi )
		{
			umiFile.Next() ;
			std::string s(umiFile.seq) ;
			if ( umiStrToInt.find(s) != umiStrToInt.end())
				umi = umiStrToInt[s] ;
			else
			{
				umi = umiStrToInt.size() ;
				umiStrToInt[s] = umi ;
			}
		}

		struct _sortRead nr ;
		int rWeight = 1 ;
		nr.read = strdup( reads.seq ) ;
		nr.id = strdup( reads.id ) ;
		if ( reads.qual != NULL )
			nr.qual = strdup( reads.qual ) ;
		else
			nr.qual = NULL ;
		nr.barcode = barcode ;
		nr.umi = umi ;

		++i ;

		if ( countMyself && i % 100000 == 0 )
			PrintLog( "Read in and count kmers for %d reads.", i ) ;
		else if ( !countMyself && i % 1000000 == 0 )
			PrintLog( "Read in %d reads.", i ) ;
		
		if ( firstReadLen == -1 )
		{
			firstReadLen = strlen( reads.seq ) ;
		}

		struct _sortRead mateR ;
		mateR.read = NULL ;
		if ( mateReads.Next() )
		{
			mateR.read = strdup( mateReads.seq ) ;
			mateR.id = strdup( mateReads.id ) ;
			mateR.barcode = barcode ;
			mateR.umi = umi ;
			if ( mateReads.qual != NULL )
				mateR.qual = strdup( mateReads.qual ) ;
			else
				mateR.qual = NULL ;
			
			++i ;

			if ( countMyself && i % 100000 == 0 )
				PrintLog( "Read in and count kmers for %d reads.", i ) ;
			else if ( !countMyself && i % 1000000 == 0 )
				PrintLog( "Read in %d reads.", i ) ;
			
			// Check chimeric. After reverse-complement, the mate should be after the current read.
			int flen = strlen( mateReads.seq ) ;
			int slen = strlen( nr.read ) ;
			seqSet.ReverseComplement( mateR.read, mateReads.seq, flen ) ;
			int minOverlap = ( flen + slen ) / 10 ;
			int minOverlap2 = ( flen + slen ) / 20 ;
			if ( minOverlap > 31 )
				minOverlap = 31 ;
			if ( minOverlap2 > 31 )
				minOverlap2 = 31 ;
			int offset = -1 ;
			int bestMatchCnt = -1 ;
			
			/*if ( AlignAlgo::IsMateOverlap( mateR.read, flen, nr.read, slen, minOverlap, offset, bestMatchCnt ) >= 0 )
			{
				printf( "Outie\n%s\n%s\n", nr.read, mateR.read ) ;
			}
			else if ( AlignAlgo::IsMateOverlap( nr.read, slen, mateR.read, flen, minOverlap, offset, bestMatchCnt ) >= 0 )
			{
				printf( "Innie\n%s\n%s\n", nr.read, mateR.read ) ;
			}
			else
				printf( "Uncertain\n" ) ;*/
			
			int overlapSize = AlignAlgo::IsMateOverlap( mateR.read, flen, nr.read, slen, minOverlap, offset, bestMatchCnt, false ) ;
			if ( overlapSize >= 0 )
			{
				// Wrong order happened, 
				// Only keep the overlapped portion.
				nr.read[ overlapSize ] = '\0' ;
				if ( nr.qual != NULL )
				{
					nr.qual[ overlapSize ] = '\0' ;

					for ( j = 0 ; j < overlapSize ; ++j )
					{
						if ( mateR.qual[j + offset] > nr.qual[j] || nr.read[j] == 'N' )
						{
							nr.read[j] = mateR.read[j + offset] ;
							nr.qual[j] = mateR.qual[j + offset] ;
						}
					}
				}

				// filter the mate.
				free( mateR.read ) ; free( mateR.id ) ; free( mateR.qual ) ;
				mateR.read = NULL ;
			}
			else if (  ( overlapSize = 
				AlignAlgo::IsMateOverlap( nr.read, slen, mateR.read, flen, minOverlap2, offset, bestMatchCnt ) ) >= 0 )
			{
				if ( bestMatchCnt >= 0.95 * overlapSize )
				{
					char *r = (char *)malloc( sizeof( char ) * ( slen + flen + 1 ) ) ;
					char *q = (char *)malloc( sizeof( char ) * ( slen + flen + 1 ) ) ;
					for ( j = 0 ; j < flen ; ++j )
					{
						r[ offset + j ] = mateR.read[j] ;
						q[ offset + j ] = mateR.qual[j] ;
					}
					int len = offset + j ;
					for ( j = 0 ; j < slen && j < len ; ++j )
					{
						if ( j < offset || nr.qual[j] >= q[j] - 14 || r[j] == 'N' )
						{
							r[j] = nr.read[j] ;
							q[j] = nr.qual[j] ;
						}
					}

					//if ( j > len ) 
					//	len = j ;
					r[len] = q[len] = '\0' ;
					
					if ( 0 )//len > 4 )
					{
						for ( j = 2 ; j <= len ; ++j )
						{
							r[j - 2] = r[j] ;
							q[j - 2] = q[j] ;
						}
						r[len - 4] = '\0' ;
						q[len - 4] = '\0' ;
						len -= 4 ;
					}
					
					free( mateR.read ) ; free( mateR.id ) ; free( mateR.qual ) ;
					mateR.read = NULL ;
					free( nr.read ) ; free( nr.qual ) ;
					nr.read = r ;
					nr.qual = q ;
					++rWeight ;
				}
				else
				{
					// Discard one of the read
					bool useFirst = true ;
					if ( nr.qual != NULL )
					{
						double avgQualR = 0, avgQualMate = 0 ;
						for ( j = offset ; j < slen ; ++j )
							avgQualR += nr.qual[j] - 32 ;
						for ( j = flen - 1 ; j >= flen - overlapSize ; --j )
							avgQualMate += mateR.qual[j] - 32 ;

						avgQualR /= overlapSize ; avgQualMate /= overlapSize ;
						if ( avgQualR + 10 < avgQualMate )
							useFirst = false ;
					}
					if ( useFirst )
					{
						free( mateR.read ) ; free( mateR.id ) ; free( mateR.qual ) ;
						mateR.read = NULL ;
					}
					else
					{
						free( nr.read ) ; free( nr.qual ) ;
						free( mateR.id ) ;

						nr.read = mateR.read ;
						nr.qual = mateR.qual ;
						strcpy( nr.read, mateReads.seq ) ;
						mateR.read = NULL ;
					}
				}
			}
			else
			{
				strcpy( mateR.read, mateReads.seq ) ;
			}
		}
		else if ( hasMate ) 
		{
			fprintf( stderr, "The two mate-pair read files have different number of reads.\n" ) ;
			exit( 1 ) ;
		}
		
		
		if ( !IsLowComplexity( reads.seq ) )
		{
			sortedReads.push_back( nr ) ;
			if ( countMyself )
				kmerCount.AddCount( nr.read ) ;
			
			if ( rWeight == 2 )
			{
				struct _sortRead wr ;
				wr = nr ;
				wr.read = strdup( nr.read ) ;
				if ( nr.qual != NULL )
					wr.qual = strdup( nr.qual ) ;
				else
					wr.qual = NULL ;
				int len = strlen( nr.id ) ;
				wr.id = ( char * )malloc( sizeof( char ) * ( len + 3 ) ) ;
				strcpy( wr.id, nr.id ) ;
				wr.id[len] = '.' ;
				wr.id[len + 1] = '1' ;
				wr.id[len + 2] = '\0' ;

				sortedReads.push_back( wr ) ;
				if ( countMyself )
					kmerCount.AddCount( wr.read ) ;
			}

			int len = strlen( nr.read ) ;
			if ( len > maxReadLen )
				maxReadLen = len ;
		}
		else
		{
			free( nr.read ) ; free( nr.id ) ; free( nr.qual ) ;
		}

		if ( mateR.read != NULL )
		{
			if ( !IsLowComplexity( mateReads.seq ) )
			{
				sortedReads.push_back( mateR ) ;
				if ( countMyself )
					kmerCount.AddCount( mateR.read ) ;

				
				int len = strlen( mateR.read ) ;
				if ( len > maxReadLen )
					maxReadLen = len ;
			}
			else
			{
				free( mateR.read ) ; free( mateR.id ) ; free( mateR.qual ) ;
			}
		}
	
		/*if ( countMyself && i % 100000 == 0 )
			PrintLog( "Read in and count kmers for %d reads.", i ) ;
		else if ( !countMyself && i % 1000000 == 0 )
			PrintLog( "Read in %d reads.", i ) ;*/
	}
	if ( maxReadLen <= 0 )
	{
		return 0 ;
	}
#ifdef DEBUG
	printf( "Finish read in the reads and kmer count.\n") ;
#endif
	int readCnt = sortedReads.size() ;
	kmerCount.SetBuffer( maxReadLen ) ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		if ( trimLevel == 0 )
			kmerCount.GetCountStatsAndTrim( sortedReads[i].read, NULL, 
					sortedReads[i].minCnt, sortedReads[i].medianCnt, sortedReads[i].avgCnt ) ;
		else
			kmerCount.GetCountStatsAndTrim( sortedReads[i].read, sortedReads[i].qual, 
					sortedReads[i].minCnt, sortedReads[i].medianCnt, sortedReads[i].avgCnt ) ;

		if ( sortedReads[i].qual != NULL )
		{
			free( sortedReads[i].qual ) ;
			sortedReads[i].qual = NULL ;
		}
		if ( sortedReads[i].read[0] == '\0' )
		{
			free( sortedReads[i].read ) ;
			free( sortedReads[i].id ) ;
			sortedReads[i].read = NULL ;
		}
	}
	
	k = 0 ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		if ( sortedReads[i].read == NULL )
			continue ;
		sortedReads[k] = sortedReads[i] ;
		sortedReads[k].len = strlen( sortedReads[i].read ) ;
		++k ;
	}
	readCnt = k ;
	sortedReads.resize( k ) ;

	PrintLog( "Found %i reads.", k ) ;
#ifdef DEBUG
	printf( "Finish put in the read kmer count.\n" ) ;
#endif
	
	kmerCount.Release() ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		sortedReads[i].info = i ;
		sortedReads[i].mateIdx = -1 ;
	}
	for ( i = 0 ; i < readCnt - 1 ; ++i )
	{
		if ( !strcmp( sortedReads[i].id, sortedReads[i + 1].id ) )
		{
			sortedReads[i].mateIdx = i + 1 ;
			sortedReads[i + 1].mateIdx = i ;
			++i ;	
		}
	}
	std::sort( sortedReads.begin(), sortedReads.end() ) ;
	

	SimpleVector<int> originToSortedIdx ;
	SimpleVector<bool> goodCandidate ; // Use mate pair alignment information to infer 
					   // whether this read could be a good candidate
	originToSortedIdx.ExpandTo( readCnt ) ;
	goodCandidate.ExpandTo( readCnt ) ;
	for ( i = 0 ; i < readCnt ; ++i )
		originToSortedIdx[ sortedReads[i].info ] = i ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		if ( sortedReads[i].mateIdx != -1 )
			sortedReads[i].mateIdx = originToSortedIdx[ sortedReads[i].mateIdx ] ;
		goodCandidate[i] = false ;
	}
	PrintLog( "Finish sorting the reads." ) ;
	

	// Quickly annoate the reads.
	if ( trimLevel > 1 )
		refSet.SetRadius(0) ;
	if ( threadCnt <= 1 )
	{
		struct _overlap geneOverlap[4] ;
		for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( i == 0 || strcmp( sortedReads[i].read, sortedReads[i - 1].read ) )
				refSet.AnnotateRead( sortedReads[i].read, 0, geneOverlap, NULL, NULL ) ;
			for ( j = 0 ; j < 4 ; ++j )
				sortedReads[i].geneOverlap[j] = geneOverlap[j] ;
		}
	}
	else
	{
		pthread_t *threads = new pthread_t[ threadCnt ] ;
		struct _quickAnnotateReadsThreadArg *args = new struct _quickAnnotateReadsThreadArg[threadCnt] ;
		pthread_attr_t attr ;
		
		pthread_attr_init( &attr ) ;
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;
		
		for ( i = 0 ; i < threadCnt ; ++i )
		{
			args[i].tid = i ;
			args[i].threadCnt = threadCnt ;
			args[i].refSet = &refSet ;
			args[i].pSortedReads = &sortedReads ;
			args[i].readCnt = readCnt ;
			pthread_create( &threads[i], &attr, QuickAnnotateReads_Thread, (void *)( args + i ) ) ;
		}

		for ( i = 0 ; i < threadCnt ; ++i )
			pthread_join( threads[i], NULL ) ;
	

		delete[] threads ;
		delete[] args ;
	}
	PrintLog( "Finish rough annotations." ) ;

	if ( trimLevel > 1 && !hasBarcode )
	{
		// use V gene assignment as barcode to speed up.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( sortedReads[i].geneOverlap[0].seqIdx != -1 && sortedReads[i].geneOverlap[0].similarity > 0.95 )
			{
				sortedReads[i].barcode = sortedReads[i].geneOverlap[0].seqIdx ;
				if ( sortedReads[i].mateIdx != -1 )
					sortedReads[ sortedReads[i].mateIdx ].barcode = sortedReads[i].geneOverlap[0].seqIdx ;
			}
		}
		/*for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( sortedReads[i].barcode == -1 )
			{
				free( sortedReads[i].id ) ;
				free( sortedReads[i].read ) ;
				if ( sortedReads[i].qual != NULL )
					free( sortedReads[i].qual ) ;
				sortedReads[i].read = NULL ;
			}
		}*/

		/*for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( sortedReads[i].geneOverlap[0].seqIdx == -1 && sortedReads[i].geneOverlap[2].seqIdx == -1 )
			{
				free( sortedReads[i].id ) ;
				free( sortedReads[i].read ) ;
				if ( sortedReads[i].qual != NULL )
					free( sortedReads[i].qual ) ;
				sortedReads[i].read = NULL ;
			}
		}*/
	}

	// Remove the redudant sequence before V gene.
	for ( i = 0 ; i < readCnt ; ++i )
	{
		if ( sortedReads[i].read == NULL )
			continue ;
		struct _overlap *geneOverlap = sortedReads[i].geneOverlap ;
		if ( geneOverlap[0].seqIdx == -1 )
			continue ;
		bool mayTrim = false ;
		if ( geneOverlap[0].seqStart < 31 && geneOverlap[0].similarity > 0.9 )
			mayTrim  = true ;

		if ( geneOverlap[0].similarity > 0.95 && 
			geneOverlap[0].seqStart <= refSet.GetSeqConsensusLen( geneOverlap[0].seqIdx ) / 3 )
			mayTrim = true ; 
		
		if ( trimLevel > 1 )
			mayTrim = true ;

		if ( !mayTrim )
			continue ;
		
		int trimBase = geneOverlap[0].readStart ;
		if ( trimLevel > 1 && refSet.GetSeqName(geneOverlap[0].seqIdx)[0] == 'T') // Some bad annotation
		{	
			if ( geneOverlap[0].seqEnd + 25 < refSet.GetSeqConsensusLen( geneOverlap[0].seqIdx ) )
				trimBase = geneOverlap[0].readEnd + 2 ;
			else if ( geneOverlap[0].similarity < 0.97 )
			{
				trimBase = ( geneOverlap[0].readStart + geneOverlap[0].readEnd ) / 2 ;
			}
		}
		//printf("%s\n", sortedReads[i].read) ;	
		//printf("0: %lf %d %d %s: %d %d\n", geneOverlap[0].similarity, mayTrim, trimBase,  
		//	refSet.GetSeqName(geneOverlap[0].seqIdx), geneOverlap[0].readStart, geneOverlap[0].readEnd) ;
		if ( trimBase <= 0 )
			continue ;

		if ( geneOverlap[2].seqIdx != -1 
			&& geneOverlap[2].readStart < trimBase && trimLevel <= 1 )
			continue ;
		if ( geneOverlap[3].seqIdx != -1 
			&& geneOverlap[3].readStart < trimBase && trimLevel <= 1 )
			continue ;

		if ( sortedReads[i].len - trimBase < 31 )
		{
			free( sortedReads[i].id ) ;
			free( sortedReads[i].read ) ;
			if ( sortedReads[i].qual != NULL )
				free( sortedReads[i].qual ) ;
			sortedReads[i].read = NULL ;
			continue ;
		}
		
		if ( geneOverlap[0].strand >= 0 )
		{
			// Remove the first bases
			int len = sortedReads[i].len ;
			for ( j = trimBase ; j <= len ; ++j )
			{
				sortedReads[i].read[j - trimBase] = sortedReads[i].read[j] ;
				if ( sortedReads[i].qual != NULL )
					sortedReads[i].qual[j - trimBase] = sortedReads[i].qual[j] ;

			}
		}
		else
		{
			// Remove the last bases.
			int len = sortedReads[i].len ;
			sortedReads[i].read[len - trimBase] = '\0' ;
			if ( sortedReads[i].qual != NULL )
				sortedReads[i].qual[len - trimBase] = '\0' ;
		}

		for ( j = 0 ; j < 4 ; ++j )
		{
			if ( geneOverlap[j].seqIdx == -1 )
				continue ;
			geneOverlap[j].readStart -= trimBase ;
			geneOverlap[j].readEnd -= trimBase ;
			if ( geneOverlap[j].readStart < 0 )
				geneOverlap[j].readStart = 0 ;
			if ( geneOverlap[j].readEnd < 0 )
			{
				geneOverlap[j].readEnd = 0 ;
				geneOverlap[j].seqIdx = -1 ;
			}
		}
		sortedReads[i].len -= trimBase ;
	}
	
	// Remove redundant sequences after constant gene ; a
	// and also constant gene for non IGH genes for long reads.
	for ( i = 0 ; i < readCnt ; ++i )
	{
		int len = sortedReads[i].len ;
		int gidx = 3 ;
		struct _overlap *geneOverlap = sortedReads[i].geneOverlap ;
		if ( sortedReads[i].read == NULL )
			continue ;
		for ( gidx = 2 ; gidx <= 3 ; ++gidx)
			if ( geneOverlap[gidx].seqIdx != -1 )
				break ;
		if ( gidx > 3 )
			continue ;
		
		if ( gidx == 2 && refSet.GetSeqName(geneOverlap[gidx].seqIdx)[2] == 'H' )
		{
			gidx = 3 ;
			if ( geneOverlap[gidx].seqIdx == -1 )
				continue ;
		}

		bool mayTrim = false ;
		if ( gidx == 3 && geneOverlap[3].seqStart < indexKmerLength && geneOverlap[3].similarity > 0.95 )
				//&& geneOverlap[3].readEnd + indexKmerLength >= sortedReads[i].len 
				//&& refSet.GetSeqName( geneOverlap[3].seqIdx )[2] != 'H' )
			mayTrim  = true ;
		if ( trimLevel > 1 )
			mayTrim = true ;
		if ( !mayTrim )
			continue ;
		
		int trimBase = len - geneOverlap[gidx].readEnd - 1 ;
		if ( trimLevel > 1 && refSet.GetSeqName(geneOverlap[gidx].seqIdx)[0] == 'T') // Some bad annotation
		{
			if ( geneOverlap[gidx].seqStart >= 25 )
				trimBase = len - ( geneOverlap[gidx].readStart - 1 ) ;
			else if ( geneOverlap[gidx].similarity < 0.97 )
			{
				trimBase = len - ( ( geneOverlap[gidx].readStart + geneOverlap[gidx].readEnd ) / 2 ) - 1 ;
			}
		}
		//printf("%s %d %s\n", sortedReads[i].id, gidx, sortedReads[i].read) ;
		//printf("%d: %lf %d %d %s: %d %d\n", gidx, geneOverlap[gidx].similarity, mayTrim, trimBase,  
		//	refSet.GetSeqName(geneOverlap[gidx].seqIdx), geneOverlap[gidx].readStart, geneOverlap[gidx].readEnd) ;
		/*if ( firstReadLen > 200 && gidx == 3 && trimLevel < 2 && refSet.GetSeqName( geneOverlap[3].seqIdx )[2] != 'H' )
		{
			trimBase = ( len - geneOverlap[3].readStart ) + geneOverlap[3].seqStart ;
		}*/
		
		if ( trimBase <= 0 )
			continue ;

		if ( gidx == 3 && geneOverlap[2].seqIdx != -1 
				&& geneOverlap[2].readStart + trimBase >= sortedReads[i].len && trimLevel <= 1 )
			continue ;
		if ( geneOverlap[0].seqIdx != -1 
				&& geneOverlap[0].readStart + trimBase >= sortedReads[i].len && trimLevel <= 1 )
			continue ;

		if ( sortedReads[i].len - trimBase < 31 )
		{
			free( sortedReads[i].id ) ;
			free( sortedReads[i].read ) ;
			if ( sortedReads[i].qual != NULL )
				free( sortedReads[i].qual ) ;
			sortedReads[i].read = NULL ;
			continue ;
		}

		if ( geneOverlap[gidx].strand < 0 )
		{
			// Remove the first bases
			for ( j = trimBase ; j <= len ; ++j )
			{
				sortedReads[i].read[j - trimBase] = sortedReads[i].read[j] ;
				if ( sortedReads[i].qual != NULL )
					sortedReads[i].qual[j - trimBase] = sortedReads[i].qual[j] ;

			}

			geneOverlap[3].seqIdx = -1 ;
		}
		else
		{
			// Remove the last bases.
			sortedReads[i].read[len - trimBase] = '\0' ;
			if ( sortedReads[i].qual != NULL )
				sortedReads[i].qual[len - trimBase] = '\0' ;
			geneOverlap[3].seqIdx = -1 ;
		}
		for ( j = 0 ; j < 4 ; ++j )
		{
			if ( geneOverlap[j].seqIdx == -1 )
				continue ;
			if ( geneOverlap[j].readStart + trimBase >= len )
			{
				geneOverlap[j].readStart = len - 1 ;
				geneOverlap[j].seqIdx = -1 ;
			}
			if ( geneOverlap[j].readEnd + trimBase >= len )
				geneOverlap[j].readEnd = len - 1 ;
		}
		//printf("%s\n", sortedReads[i].read);
		sortedReads[i].len -= trimBase ;
	}

	// Remove the reads that are too short if this is a long read data set.
	if ( firstReadLen > 200 )
	{
		for ( i = 0 ; i < readCnt ; ++i )
			if ( sortedReads[i].read != NULL && sortedReads[i].len < firstReadLen / 3 )
			{
				free( sortedReads[i].id ) ;
				free( sortedReads[i].read ) ;
				if ( sortedReads[i].qual != NULL )
					free( sortedReads[i].qual ) ;
				sortedReads[i].read = NULL ;
				continue ;
			}

		seqSet.SetIsLongSeqSet( true ) ;
	}
	
	// Remove the sequences whose J gene is not associate with C gene.
	/*for ( i = 0 ; i < readCnt ; ++i )
	  {
	  struct _overlap *geneOverlap = sortedReads[i].geneOverlap ;
	  if ( sortedReads[i].read == NULL )
	  continue ;

	  if ( geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 && geneOverlap[3].seqIdx == -1 )
	  {
	  if ( geneOverlap[2].similarity > 0.9 && 
	  sortedReads[i].len - geneOverlap[2].readEnd > 50 ) // With 50 base but could not identify C gene
	  {
	  free( sortedReads[i].id ) ;
	  free( sortedReads[i].read ) ;
	  if ( sortedReads[i].qual != NULL )
	  free( sortedReads[i].qual ) ;
	  sortedReads[i].read = NULL ;
	  }
	  }
	  }*/

	k = 0 ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		if ( sortedReads[i].read != NULL )
		{
			originToSortedIdx[i] = k ;
			sortedReads[k] = sortedReads[i] ;
			++k ;
		}
		else
			originToSortedIdx[i] = -1 ;
	}

	// Don't forget to update the mateidx.
	for ( i = 0 ; i < k ; ++i )
	{
		if ( sortedReads[i].mateIdx != -1 )
			sortedReads[i].mateIdx = originToSortedIdx[ sortedReads[i].mateIdx ] ;
	}

	originToSortedIdx.Release() ;
	sortedReads.resize( k ) ;
	readCnt = k ;

	std::vector<int> rescueReadIdx ;
	std::vector<int> assembledReadIdx ;
	int assembledReadCnt = 0 ;
	int prevAddRet = -1 ;
	
	/*if ( seqSet.GetSeqCnt() > 0 )
	{
		for ( i = 0 ; i < readCnt ; ++i )
			assembledReadIdx.push_back( i ) ;
		readCnt = 0 ;
	}*/
	
	if ( firstReadLen / 2 < 31 )
	{
		int l = firstReadLen / 2 ;
		if ( l < 21 )
			l = 21 ;
		seqSet.SetHitLenRequired( l ) ;
	}

	if ( hasBarcode )
		seqSet.SetHitLenRequired( 17 ) ;

	if ( firstReadLen > 200 || trimLevel > 1 )
		changeKmerLengthThreshold /= 2 ;

	for ( i = 0 ; i < readCnt ; ++i )
	{
		/*++assembledReadCnt ;
		assembledReadIdx.push_back( i ) ;
		continue ;*/
		static struct _overlap geneOverlap[4] ;

#ifdef DEBUG
		printf( "%s %s %d %lf\n", sortedReads[i].id, sortedReads[i].read, sortedReads[i].minCnt, sortedReads[i].avgCnt ) ;
		fflush( stdout ) ;
#endif
		int addRet = -1 ;
		
		if ( i == 0 || strcmp( sortedReads[i].read, sortedReads[i - 1].read ) 
			|| sortedReads[i].barcode != sortedReads[i - 1].barcode )
		{
			//printf( "new stuff\n" ) ;
			//buffer[0] = '\0' ;
			//refSet.AnnotateRead( sortedReads[i].read, 0, geneOverlap, NULL, NULL ) ;
			for ( j = 0 ; j < 4 ; ++j )
				geneOverlap[j] = sortedReads[i].geneOverlap[j] ;
			
			// If the order of V,D,J,C is wrong from this read, then we ignore this.
			//   probably from read through in cyclic fragment.
			// We do not need to worry about the strand here, since if the strand is -1, the coordinate
			//   is relative to the reversed-complementary read.
			bool filter = false ;
			int strand = 0 ;
			/*for ( j = 0 ; j < 4 ; ++j )
			{
				if ( geneOverlap[j].seqIdx == -1 )
					continue ;
				printf( "%d: %s %d %d; %d %d; %d %lf\n", j, refSet.GetSeqName( geneOverlap[j].seqIdx ), 
					geneOverlap[j].readStart, geneOverlap[j].readEnd,
					geneOverlap[j].seqStart, geneOverlap[j].seqEnd,
					geneOverlap[j].strand, geneOverlap[j].similarity ) ;
			}*/
			for ( j = 0 ; j < 4 ; ++j )
			{
				int l ;
				if ( geneOverlap[j].seqIdx == -1 )
					continue ;
				for ( l = j + 1 ; l < 4 ; ++l )
				{
					if ( geneOverlap[l].seqIdx == -1 )
						continue ;
					
					if ( geneOverlap[j].readEnd - 10 > geneOverlap[l].readStart ) 
					{
						filter = true ;	
						break ;
					}
				}
				if ( filter )
					break ;
			}

			if ( geneOverlap[3].seqIdx != -1 && geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx == -1 ) // From constant gene.
			{
				//if ( geneOverlap[3].readEnd - geneOverlap[3].readStart + 1 < sortedReads[i].len )
				//{
				if ( geneOverlap[3].seqStart >= constantGeneEnd )
					filter = true ;
				else if ( geneOverlap[3].seqStart >= 100 && ( geneOverlap[3].strand == 1 
					|| geneOverlap[3].readEnd - geneOverlap[3].readStart + 1 < sortedReads[i].len ) )
					filter = true ;
				//}
				
			}
			/*for ( j = 0 ; j < 4 ; ++j )
				printf( "%d ", geneOverlap[j].seqIdx ) ;
			printf( "\n" ) ;*/
			if ( filter ) 
				addRet = -1 ;
			else
			{
				char name[5] ;
				name[0] = '\0' ;
				strand = 0 ;
				int ambiguousStrand = 0 ;
				for ( j = 0 ; j < 4 ; ++j )
					if ( geneOverlap[j].seqIdx != -1 )
					{
						char *s = refSet.GetSeqName( geneOverlap[j].seqIdx ) ;
						name[0] = s[0] ; name[1] = s[1] ; name[2] = s[2] ; name[3] = s[3] ;
						name[4] = '\0' ;
						if ( strand != 0 && strand != geneOverlap[j].strand )
							ambiguousStrand = 1 ;
						strand = geneOverlap[j].strand ;
					}
				if ( ambiguousStrand )
					strand = 0 ;

				double similarityThreshold = 0.9 ;
				if ( sortedReads[i].minCnt >= 20 )
					similarityThreshold = 0.97 ;
				else if ( sortedReads[i].minCnt >= 2
					|| ( sortedReads[i].minCnt >= 5 && firstReadLen > 200 ) )
				{
					//double tmp = 1.0 - 4.0 / sortedReads[i].len ;
					//similarityThreshold = tmp < 0.95 ? tmp : 0.95 ;
					//if ( similarityThreshold < 0.9 )
					//	similarityThreshold = 0.9 ;
					similarityThreshold = 0.95 ;
				}

				if ( name[0] == 'T' && similarityThreshold < 0.95 ) // TCR
					similarityThreshold = 0.95 ;
				
				// When using barcode, we could be more casual.
				if ( hasBarcode || trimLevel > 1 )
					similarityThreshold = 0.9 ;

				//if ( similarityThreshold > 0.9 && geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 
				//	&& geneOverlap[0].readEnd < geneOverlap[2].readStart )
				//	similarityThreshold = 1.0 ;

				addRet = seqSet.AddRead( sortedReads[i].read, name, strand, sortedReads[i].barcode, 
						sortedReads[i].minCnt, trimLevel > 1, similarityThreshold ) ;
				
				if ( addRet < 0 )
				{
					// Check the overall hit length.
					int matchCnt = 0 ;
					for ( j = 0 ; j < 4 ; ++j )
					{
						if ( geneOverlap[j].seqIdx != -1 )
							matchCnt += geneOverlap[j].matchCnt / 2 ;
					}
					filter = true ;
					if ( matchCnt >= 31 )
						filter = false ;
					else
					{
						// The anchor is very close to the end of V or J gene
						if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 
							&& geneOverlap[0].readEnd < geneOverlap[2].readStart )
						{
							filter = false ;
						}
						else if ( geneOverlap[0].seqIdx != -1 )
						{
							if ( geneOverlap[0].seqEnd >= 
								refSet.GetSeqConsensusLen( geneOverlap[0].seqIdx ) - 17 )
								filter = false ;
						}
						else if ( geneOverlap[2].seqIdx != -1 )
						{
							if ( geneOverlap[2].seqStart <= 17 )
								filter = false ;
						}
					}
					for ( j = 0 ; j < 4 ; ++j )
						if ( geneOverlap[j].seqIdx != -1 )
						{
							break ;
						}
					
					if ( !filter )
					{
						addRet = seqSet.InputNovelRead( refSet.GetSeqName( geneOverlap[j].seqIdx ), 
							sortedReads[i].read, geneOverlap[j].strand, sortedReads[i].barcode ) ;

					}
					else if ( goodCandidate[i] )
					{
						// The mate is matched well, we use motif to check whether this could be
						//   from CDR3.
						if ( seqSet.HasMotif( sortedReads[i].read, -sortedReads[ sortedReads[i].mateIdx ].strand ) )
							addRet = seqSet.InputNovelRead( "Novel", sortedReads[i].read, 
									-sortedReads[ sortedReads[i].mateIdx ].strand,
									sortedReads[i].barcode ) ;
					}
					//printf( "hello %d %d. %d %d %d %d\n", i, addRet, geneOverlap[0].seqIdx,
					//		geneOverlap[1].seqIdx, geneOverlap[2].seqIdx, geneOverlap[3].seqIdx ) ;
				}
				sortedReads[i].strand = strand ;
			}
		}
		else
		{
			//printf( "saved time\n" ) ;
			if ( prevAddRet != -1 && prevAddRet != -3 )
				addRet = seqSet.RepeatAddRead( sortedReads[i].read ) ;
			else if ( prevAddRet == -3 )
				addRet = -3 ;

			sortedReads[i].strand = sortedReads[i - 1].strand ;
		}
		
		if ( addRet == -2 )
			rescueReadIdx.push_back( i ) ;
		else if ( addRet >= 0 )
		{
			++assembledReadCnt ;
			assembledReadIdx.push_back( i ) ;

			if ( sortedReads[i].mateIdx > i ) // This handles the case mateIdx = -1
			{
				bool good = false ;
				bool maySpan = false ;
				if ( geneOverlap[0].seqIdx != -1 && geneOverlap[0].similarity >= 0.9 
						&& sortedReads[i].strand == 1 )
				{
					good = true ;
					if ( geneOverlap[2].seqIdx != -1 && geneOverlap[2].readStart > geneOverlap[0].readEnd )
						maySpan = true ;
					if ( geneOverlap[3].seqIdx != -1 && geneOverlap[3].readStart > geneOverlap[0].readEnd )
						maySpan = true ;
				}

				for ( j = 2 ; j <= 3 ; ++j )
				{
					if ( geneOverlap[j].seqIdx != -1 && geneOverlap[j].similarity >= 0.9 
							&& sortedReads[i].strand == -1 )
					{
						good = true ;
						if ( geneOverlap[0].seqIdx != -1 && geneOverlap[j].readStart > geneOverlap[0].readEnd )
							maySpan = true ;
					}

				}

				if ( maySpan ) // If the read can span V,J, then there is no need to force add its mate.
					good = false ;

				goodCandidate[ sortedReads[i].mateIdx ] = good ;
			}
		}

		if ( assembledReadCnt > 0 && assembledReadCnt % 10000 == 0 )
			seqSet.UpdateAllConsensus() ;

		if ( ( i + 1 ) % 100000 == 0 )
			PrintLog( "Processed %d reads (%d are used for assembly).", i + 1, assembledReadCnt ) ;
		
		prevAddRet = addRet ;
#ifdef DEBUG
		printf( "done\n" ) ;
#endif

		if ( seqSet.Size() > changeKmerLengthThreshold && indexKmerLength < 16 && !hasBarcode )
		{
			changeKmerLengthThreshold *= 4 ;
			indexKmerLength += 2 ;
			seqSet.ChangeKmerLength( indexKmerLength ) ;
		}
	}
	seqSet.UpdateAllConsensus() ;
	PrintLog(  "Assembled %d reads.", assembledReadCnt ) ;

	// Go through the second round.
	// TODO: user-defined number of rounds.
	/*seqSet.ResetPosWeight() ;
	reads.Rewind() ;
	while ( reads.Next() )
	{
		//printf( "%s %s\n", reads.id, reads.seq ) ;
		//fflush( stdout ) ;
		seqSet.AddRead( reads.seq ) ;
		//printf( "done\n" ) ;
	}*/

	int rescueReadCnt = rescueReadIdx.size() ;
	if ( firstReadLen > 200 )
		rescueReadCnt = 0 ;
	PrintLog( "Try to rescue %d reads for assembly.", rescueReadCnt) ;
	assembledReadCnt = 0 ;
	for ( i = 0 ; i < rescueReadCnt ; ++i )
	{
#ifdef DEBUG
		printf( "%s %s %d %lf\n", sortedReads[ rescueReadIdx[i] ].id, 
			sortedReads[ rescueReadIdx[i] ].read, sortedReads[  rescueReadIdx[i]  ].medianCnt, 
			sortedReads[ rescueReadIdx[i] ].avgCnt ) ;
		fflush( stdout ) ;
#endif
		int addRet = -1 ;
		char name[2] = "" ;
		
		double similarityThreshold = 0.9 ;
		if ( sortedReads[ rescueReadIdx[i] ].minCnt >= 20 )
			similarityThreshold = 0.97 ;
		else if ( sortedReads[ rescueReadIdx[i] ].minCnt >= 2 )
		{
			//double tmp = 1.0 - 4.0 / sortedReads[i].len ;
			//similarityThreshold = tmp < 0.95 ? tmp : 0.95 ;
			//if ( similarityThreshold < 0.9 )
			//	similarityThreshold = 0.9 ;
			similarityThreshold = 0.95 ;
		}
		
		int strand = 0 ;
		addRet = seqSet.AddRead( sortedReads[ rescueReadIdx[i] ].read, name, strand, sortedReads[ rescueReadIdx[i] ].barcode, 
			1, trimLevel > 1 ? true : false, similarityThreshold ) ;
		sortedReads[ rescueReadIdx[i] ].strand = strand ;

		if ( addRet >= 0 )
		{
			++assembledReadCnt ;
			assembledReadIdx.push_back( rescueReadIdx[i] ) ;
		}
#ifdef DEBUG
		printf( "done\n" ) ;
#endif
	}
	seqSet.UpdateAllConsensus() ;
	PrintLog( "Rescued %d reads.", assembledReadCnt ) ;
	//exit( 1 ) ;
	/*seqSet.Clean( true ) ;
	fprintf( stderr, "Found %d raw contigs.\nStart to merge raw contigs.\n", seqSet.GetSeqCnt() ) ;
	
	i = seqSet.Assemble() ;
	seqSet.UpdateAllConsensus() ;
	fprintf( stderr, "Reduce %d raw contigs.\n", i ) ;
	
	fprintf( stderr, "Annotate the contigs.\n" ) ;
	seqSet.Annotate( refSet ) ;*/

	// Output the preliminary assembly.
	//seqSet.Clean( true ) ;
	FILE *fp ;
	if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_raw.out", outputPrefix ) ;
		fp = fopen( buffer, "w" ) ;
	}
	else
		fp = stdout ;
	
	if ( hasBarcode )
		seqSet.Output( fp, &barcodeIntToStr ) ;
	else
		seqSet.Output( fp, NULL ) ;
	fflush( fp ) ;
	
	if ( outputPrefix[0] != '-' )
		fclose( fp ) ;

	//sprintf( buffer, "%s_assembled_reads.fa", outputPrefix ) ;

	std::vector<struct _assignRead> assembledReads ;
	assembledReadCnt = assembledReadIdx.size() ;
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		struct _assignRead nr ;
		nr.id = sortedReads[ assembledReadIdx[i] ].id ;
		nr.read = sortedReads[ assembledReadIdx[i] ].read ;
		nr.barcode = sortedReads[ assembledReadIdx[i] ].barcode ;
		nr.umi = sortedReads[ assembledReadIdx[i] ].umi ;
		nr.info = assembledReadIdx[i] ;
		nr.overlap.seqIdx = -1 ;
		nr.overlap.strand = sortedReads[ assembledReadIdx[i] ].strand ;
		assembledReads.push_back( nr ) ;
	}
	//std::sort( assembledReads.begin(), assembledReads.end(), CompSortReadById ) ;	
	
	sprintf( buffer, "%s_assembled_reads.fa", outputPrefix ) ;
	fp = fopen( buffer, "w" ) ;
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		buffer[0] = '\0' ;
		if ( hasBarcode )
		{
			char *p = buffer + strlen( buffer ) ;
			sprintf( p, " barcode:%s", barcodeIntToStr[ assembledReads[i].barcode ].c_str() ) ;
		}
		if ( hasUmi )
		{
			char *p = buffer + strlen( buffer ) ;
			sprintf( p, " umi:%d", assembledReads[i].umi ) ;
		}

		fprintf( fp, ">%s %d %d %d%s\n%s\n", assembledReads[i].id, assembledReads[i].overlap.strand,
			sortedReads[ assembledReadIdx[i] ].minCnt, 
			sortedReads[ assembledReadIdx[i] ].medianCnt, buffer,
			assembledReads[i].read ) ;
	}
	fclose( fp ) ;

	if ( skipMateExtension )
	{
		FILE *fp ;
		if ( outputPrefix[0] != '-' )
		{
			sprintf( buffer, "%s_final.out", outputPrefix ) ;
			fp = fopen( buffer, "w" ) ;
		}
		else
			fp = stdout ;

		if ( hasBarcode )
			seqSet.Output( fp, &barcodeIntToStr ) ;
		else
			seqSet.Output( fp, NULL ) ;
		fflush( fp ) ;

		if ( outputPrefix[0] != '-' )
			fclose( fp ) ;

		return 0 ;	
	}

	SeqSet extendedSeq( indexKmerLength > 17 ? indexKmerLength : 17 ) ;
	extendedSeq.InputSeqSet( seqSet, false ) ;
	struct _overlap assign ;
	memset( &assign, -1, sizeof( assign ) ) ;

	/*fp = fopen( "assign_results.out", "r" ) ;
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		fread( &assembledReads[i].overlap, sizeof( assembledReads[i].overlap ), 1, fp ) ;
		//printf( "%d %d %d\n", assembledReads[i].overlap.seqIdx, assembledReads[i].overlap.readStart, 
		//	assembledReads[i].overlap.readEnd ) ;
	}*/
	//fp = fopen( "assign_results.out", "w" ) ;

	// Add the distant constant genes from the annotation
	/*int refSetSize = refSet.GetSeqCnt() ;
	int constantGeneAddCnt = 0 ;
	for ( i = 0 ; i < refSetSize ; ++i )
	{
		char *name = GetSeqName( i ) ;
		if ( refSet.GetGeneType( name ) != 3 || refSet.GetSeqConsensusLen( i ) <= constantGeneEnd )
			continue ;
		extendedSeq.InputRefSeq( name, refSet.GetSeqConsensus( i ) + constantGeneEnd ) ;
	}*/
	if ( firstReadLen > 200 )
		extendedSeq.SetIsLongSeqSet( true ) ;
	extendedSeq.SetNovelSeqSimilarity( 0.95 ) ;
	if ( threadCnt <= 1 )
	{
		for ( i = 0 ; i < assembledReadCnt ; ++i )
		{
			if ( i == 0 || strcmp( assembledReads[i].read, assembledReads[i - 1].read ) )
				extendedSeq.AssignRead( assembledReads[i].read, assembledReads[i].overlap.strand, 
						assembledReads[i].barcode, assign ) ;
			assembledReads[i].overlap = assign ;	
			if ( ( i + 1 ) % 100000 == 0 )
			{
				PrintLog( "Processed %d reads for extension.", i + 1 ) ;
			}
			//fprintf( fp, "%s %s %d\n", assembledReads[i].id, assembledReads[i].read, assign.seqIdx ) ;
			//fwrite( &assign, sizeof( assign ), 1, fp ) ;
		}
	}
	else
	{
		pthread_t *threads = new pthread_t[ threadCnt ] ;
		struct _assignReadsThreadArg *args = new struct _assignReadsThreadArg[threadCnt] ;
		pthread_attr_t attr ;
		
		pthread_attr_init( &attr ) ;
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;
		
		for ( i = 0 ; i < threadCnt ; ++i )
		{
			args[i].tid = i ;
			args[i].threadCnt = threadCnt ;
			args[i].seqSet = &extendedSeq ;
			args[i].pAssembledReads = &assembledReads ;
			args[i].assembledReadCnt = assembledReadCnt ;
			pthread_create( &threads[i], &attr, AssignReads_Thread, (void *)( args + i ) ) ;
		}

		for ( i = 0 ; i < threadCnt ; ++i )
			pthread_join( threads[i], NULL ) ;
	

		delete[] threads ;
		delete[] args ;
	}
	extendedSeq.SetNovelSeqSimilarity( 0.9 ) ;
	extendedSeq.RecomputePosWeight( assembledReads ) ;
	//fclose( fp ) ;
	//exit( 1 ) ;

#ifdef DEBUG
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		printf( "%s %d %lf %d\n", assembledReads[i].id, assembledReads[i].overlap.seqIdx, assembledReads[i].overlap.similarity,
			assembledReads[i].overlap.strand ) ;
	}
#endif

	/*extendedSeq.BreakFalseAssembly( assembledReads ) ;
	
	if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_contig.out", outputPrefix ) ;
		fp = fopen( buffer, "w" ) ;
	}
	else
		fp = stdout ;

	extendedSeq.Output( fp ) ;
	fflush( fp ) ;
	
	if ( outputPrefix[0] != '-' )
		fclose( fp ) ;*/

	int avgReadLen = 1 ;
	for ( i = 0 ; i < assembledReadCnt && i < 1000 ; ++i )
		avgReadLen += strlen( assembledReads[i].read ) ;
	if ( i > 0 )
		avgReadLen /= i ;
	
	PrintLog( "Extend assemblies by mate pair information." ) ;
	extendedSeq.ExtendSeqFromReads( assembledReads, 17, refSet ) ; //( avgReadLen / 30 < 31 ) ? 31 : ( avgReadLen / 3 )  ) ;
	extendedSeq.UpdateAllConsensus() ;
		
	/*if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_extended.out", outputPrefix ) ;
		fp = fopen( buffer, "w" ) ;
	}
	else
		fp = stdout ;

	extendedSeq.Output( fp ) ;
	fflush( fp ) ;
	
	if ( outputPrefix[0] != '-' )
		fclose( fp ) ;*/
	
		
	
	
	// Now the assembled reads should be sorted by their read id.
	//PrintLog( "Rescue low frequency reads." ) ;
	SeqSet lowFreqSeqSet( indexKmerLength ) ;
	KmerCount lowFreqKmerCount( 31 ) ;
	std::vector<int> lowFreqReadsIdx ;
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		continue ;
		if ( assembledReads[i].overlap.seqIdx != -1 ) 
		//if ( assembledReads[i].overlap.seqIdx != -1 && assembledReads[i + 1].overlap.seqIdx != -1 
		//	&& assembledReads[i].overlap.seqIdx == assembledReads[i + 1].overlap.seqIdx )
		{
			continue ;
		}

		if ( sortedReads[ assembledReads[i].info ].minCnt < 4 )
			continue ;
		
		int nCnt = 0 ;
		for ( j = 0 ; assembledReads[i].read[j] ; ++j )
			if ( assembledReads[i].read[j] == 'N' )
				++nCnt ;
		if ( nCnt >= 1 )
		{
			continue ;
		}

		lowFreqKmerCount.AddCount( assembledReads[i].read ) ;
		lowFreqReadsIdx.push_back(i) ;
	}
	k = 0 ;
	int lowFreqCnt = lowFreqReadsIdx.size() ;
	std::vector<struct _sortRead> lowFreqReads ;
	for ( i = 0 ; i < lowFreqCnt ; ++i )
	{
		// At least one of the mate should have at least 2 kmer count.
		int a = lowFreqReadsIdx[i] ;
		int minCntA, tmp ;
		double tmpd ;
		lowFreqKmerCount.GetCountStatsAndTrim( assembledReads[a].read, NULL, 
			minCntA, tmp, tmpd ) ;

		// Then these two pairs must overlap with each other.
		struct _sortRead nr ;
		nr.id = NULL ;
		nr.read = assembledReads[a].read ;
		nr.qual = NULL ;
		nr.minCnt = minCntA ; 
		nr.medianCnt = tmp ;
		nr.avgCnt = tmpd ;
		nr.strand = assembledReads[a].overlap.strand ;
		lowFreqReads.push_back( nr ) ;
	}
	std::sort( lowFreqReads.begin(), lowFreqReads.end() ) ;
	/*for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		char *r = strdup( assembledReads[i].read ) ;
		lowFreqMergedReads.push_back( r ) ;
		++k ;
		
	}*/
	lowFreqCnt = lowFreqReads.size() ;
	//PrintLog( "Found %d low frequency reads.", lowFreqCnt ) ;
	for ( i = 0 ; i < lowFreqCnt ; ++i )
	{
		char name[10] = "" ;
		name[0] = '\0' ;
		//printf( "%s\n", lowFreqReads[i].read ) ;
		//fflush( stdout ) ;
		//printf( ">r%d\n%s\n", i, lowFreqReads[i].read ) ;
		extendedSeq.AssignRead( lowFreqReads[i].read, lowFreqReads[i].strand, lowFreqReads[i].barcode, assign ) ;
		if ( assign.seqIdx != -1 )
		{
			if ( assign.similarity < 1.0 )
				extendedSeq.AddAssignedRead( lowFreqReads[i].read, assign ) ;
		}
		//else if ( extendedSeq.AddRead( lowFreqReads[i].read, name, 1.0 ) == -1 )
		else
		{
			int strand ;
			if ( lowFreqReads[i].minCnt < 2 )
			{
				lowFreqSeqSet.AssignRead( lowFreqReads[i].read, lowFreqReads[i].strand, 
						lowFreqReads[i].barcode, assign ) ;
				if ( assign.seqIdx != -1 )
				{
					lowFreqSeqSet.AddAssignedRead( lowFreqReads[i].read, assign ) ;
				}
			}
			else if ( lowFreqSeqSet.AddRead( lowFreqReads[i].read, name, strand, lowFreqReads[i].barcode, 1, false, 0.97 ) < 0 )
			{
				struct _overlap geneOverlap[4] ;
				//buffer[0] = '\0' ;
				refSet.AnnotateRead( lowFreqReads[i].read, 0, geneOverlap, NULL, NULL ) ;

				int componentCnt = 0 ;
				int lastJ = -1 ;
				for ( j = 0 ; j < 4 ; ++j )
				{
					if ( geneOverlap[j].seqIdx != -1 )
					{
						++componentCnt ;
						lastJ = j ;
					}
				}
				if ( componentCnt >= 1 && lowFreqReads[i].minCnt >= 2 )
				{
					j = lastJ ;
					lowFreqSeqSet.InputNovelRead( refSet.GetSeqName( geneOverlap[j].seqIdx ), 
							lowFreqReads[i].read, geneOverlap[j].strand, lowFreqReads[i].barcode ) ;
				}

			}
		}
	}
	
	extendedSeq.UpdateAllConsensus() ;
	lowFreqSeqSet.UpdateAllConsensus() ;
	extendedSeq.InputSeqSet( lowFreqSeqSet, false ) ;	
	
	
	PrintLog( "Remove redundant assemblies." ) ;
	extendedSeq.ChangeKmerLength( 31 ) ;
	//i = extendedSeq.ExtendSeqFromSeqOverlap( ( avgReadLen / 2 < 31 ) ? 31 : ( avgReadLen / 2 ) ) ;
	i = extendedSeq.RemoveRedundantSeq() ;
	
	
	if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_final.out", outputPrefix ) ;
		fp = fopen( buffer, "w" ) ;
	}
	else
		fp = stdout ;
	
	if ( hasBarcode )
		extendedSeq.Output( fp, &barcodeIntToStr ) ;
	else
		extendedSeq.Output( fp, NULL ) ;
	fflush( fp ) ;
	
	if ( outputPrefix[0] != '-' )
		fclose( fp ) ;
	
	/*for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		fprintf( fp, ">%s\n%s\n", assembledReads[i].id, assembledReads[i].read ) ;
	}
	fclose( fp ) ;*/
	if ( readCnt == 0 )
		readCnt = assembledReadCnt ;	
	for ( i = 0 ; i < readCnt ; ++i )
	{
		free( sortedReads[i].id ) ;
		free( sortedReads[i].read ) ;
		if ( sortedReads[i].qual != NULL )
			free( sortedReads[i].qual ) ;
	}

	PrintLog( "Finish assembly." ) ;
	return 0 ;
}
