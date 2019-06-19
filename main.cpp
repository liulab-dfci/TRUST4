#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>

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
		"\t-c STRING: the path to the kmer count file\n"
		"\t--noTrim: do not trim the reads for assembly (default: trim)\n" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

char buffer[10240] = "" ;

static const char *short_options = "f:u:1:2:b:o:c:" ;
static struct option long_options[] = {
			{ "debug-ns", required_argument, 0, 10000 },
			{ "noTrim", no_argument, 0, 10001 },
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

	int strand ;

	bool operator<( const struct _sortRead &b )
	{
		if ( minCnt != b.minCnt )
			return minCnt > b.minCnt ;
		else if ( medianCnt != b.medianCnt )
			return medianCnt > b.medianCnt ;
		else if ( avgCnt != b.avgCnt )
			return avgCnt > b.avgCnt ;
		else 
			return strcmp( read, b.read ) < 0 ;
	}
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
	char outputPrefix[200] = "trust" ;

	ReadFiles reads ;
	ReadFiles mateReads ;
	bool countMyself = true ;
	int maxReadLen = -1 ;
	bool flagNoTrim = false ;
	bool hasMate = false ;

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
		else if ( c == 10000 )
		{
			seqSet.InputNovelFa( optarg ) ;
		}
		else if ( c == 10001) 
		{
			flagNoTrim =true ;
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
	std::vector< struct _sortRead > sortedReads ;
	i = 0 ;
	while ( reads.Next() )
	{
		/*struct _overlap geneOverlap[4] ;
		refSet.AnnotateRead( reads.seq, 1, geneOverlap, buffer ) ;
	
		if ( geneOverlap[0].seqIdx + geneOverlap[1].seqIdx + geneOverlap[2].seqIdx + geneOverlap[3].seqIdx == -4 )
		{
			continue ;
		}*/

		struct _sortRead nr ;
		int rWeight = 1 ;
		nr.read = strdup( reads.seq ) ;
		nr.id = strdup( reads.id ) ;
		if ( reads.qual != NULL )
			nr.qual = strdup( reads.qual ) ;
		else
			nr.qual = NULL ; 

		++i ;

		if ( countMyself && i % 100000 == 0 )
			PrintLog( "Read in and count kmers for %d reads.", i ) ;
		else if ( !countMyself && i % 1000000 == 0 )
			PrintLog( "Read in %d reads.", i ) ;
		
		struct _sortRead mateR ;
		mateR.read = NULL ;
		if ( mateReads.Next() )
		{
			mateR.read = strdup( mateReads.seq ) ;
			mateR.id = strdup( mateReads.id ) ;
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
					for ( j = 0 ; j < slen ; ++j )
					{
						if ( j < offset || nr.qual[j] >= q[j] - 14 || r[j] == 'N' )
						{
							r[j] = nr.read[j] ;
							q[j] = nr.qual[j] ;
						}
					}

					if ( j > len ) 
						len = j ;
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
		if ( flagNoTrim )
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
		++k ;
	}
	readCnt = k ;
	sortedReads.resize( k ) ;

	PrintLog( "Found %i reads.", k ) ;
#ifdef DEBUG
	printf( "Finish put in the read kmer count.\n" ) ;
#endif
	
	kmerCount.Release() ;
	std::sort( sortedReads.begin(), sortedReads.end() ) ;
	PrintLog( "Finish sorting the reads." ) ;
	
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

	for ( i = 0 ; i < readCnt ; ++i )
	{
		/*++assembledReadCnt ;
		assembledReadIdx.push_back( i ) ;
		continue ;*/

#ifdef DEBUG
		printf( "%s %s %d %lf\n", sortedReads[i].id, sortedReads[i].read, sortedReads[i].minCnt, sortedReads[i].avgCnt ) ;
		fflush( stdout ) ;
#endif
		int addRet = -1 ;
		
		if ( i == 0 || strcmp( sortedReads[i].read, sortedReads[i - 1].read ) )
		{
			//printf( "new stuff\n" ) ;
			struct _overlap geneOverlap[4] ;
			//buffer[0] = '\0' ;
			refSet.AnnotateRead( sortedReads[i].read, 0, geneOverlap, NULL, buffer ) ;
			
			// If the order of V,D,J,C is wrong from this read, then we ignore this.
			//   probably from read through in cyclic fragment. 
			bool filter = false ;
			int strand = 0 ;
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
			if ( geneOverlap[3].seqIdx != -1 && geneOverlap[3].seqStart >= 100 ) // From constant gene.
				filter = true ;

			if ( filter ) 
				addRet = -1 ;
			else
			{
				char name[5] ;
				name[0] = '\0' ;
				strand = 0 ;
				for ( j = 0 ; j < 4 ; ++j )
					if ( geneOverlap[j].seqIdx != -1 )
					{
						char *s = refSet.GetSeqName( geneOverlap[j].seqIdx ) ;
						name[0] = s[0] ; name[1] = s[1] ; name[2] = s[2] ; name[3] = s[3] ;
						name[4] = '\0' ;
						strand = geneOverlap[j].strand ;
					}

				double similarityThreshold = 0.9 ;
				if ( sortedReads[i].minCnt >= 20 )
					similarityThreshold = 0.97 ;
				else if ( sortedReads[i].minCnt >= 2 )
					similarityThreshold = 0.95 ;
				
				addRet = seqSet.AddRead( sortedReads[i].read, name, strand, similarityThreshold ) ;
				
				if ( addRet < 0 )
				{
					for ( j = 0 ; j < 4 ; ++j )
						if ( geneOverlap[j].seqIdx != -1 )
							break ;
					if ( j < 4 )
					{
						addRet = seqSet.InputNovelRead( refSet.GetSeqName( geneOverlap[j].seqIdx ), 
							sortedReads[i].read, geneOverlap[j].strand ) ;

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
			if ( prevAddRet != -1 )
				addRet = seqSet.RepeatAddRead( sortedReads[i].read ) ;
		
			sortedReads[i].strand = sortedReads[i - 1].strand ;
		}
		
		if ( addRet == -2 )
			rescueReadIdx.push_back( i ) ;
		else if ( addRet >= 0 )
		{
			++assembledReadCnt ;
			assembledReadIdx.push_back( i ) ;
		}

		if ( assembledReadCnt > 0 && assembledReadCnt % 10000 == 0 )
			seqSet.UpdateAllConsensus() ;

		if ( ( i + 1 ) % 100000 == 0 )
			PrintLog( "Processed %d reads (%d are used for assembly).", i + 1, assembledReadCnt ) ;
		
		prevAddRet = addRet ;
#ifdef DEBUG
		printf( "done\n" ) ;
#endif

		if ( seqSet.Size() > changeKmerLengthThreshold && indexKmerLength < 16 )
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
		addRet = seqSet.AddRead( sortedReads[ rescueReadIdx[i] ].read, name, 0, 
			( sortedReads[ rescueReadIdx[i] ].minCnt >= 20 ? 0.97 : 0.9 ) ) ;
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
	
	seqSet.Output( fp ) ;
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
		fprintf( fp, ">%s %d %d\n%s\n", assembledReads[i].id, sortedReads[ assembledReadIdx[i] ].minCnt, 
			sortedReads[ assembledReadIdx[i] ].medianCnt,
			assembledReads[i].read ) ;
	}
	fclose( fp ) ;

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
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		if ( i == 0 || strcmp( assembledReads[i].read, assembledReads[i - 1].read ) )
			extendedSeq.AssignRead( assembledReads[i].read, assembledReads[i].overlap.strand, 0.95, assign ) ;
		assembledReads[i].overlap = assign ;	
		if ( ( i + 1 ) % 100000 == 0 )
		{
			PrintLog( "Processed %d reads for extension.", i + 1 ) ;
		}
		//fprintf( fp, "%s %s %d\n", assembledReads[i].id, assembledReads[i].read, assign.seqIdx ) ;
		//fwrite( &assign, sizeof( assign ), 1, fp ) ;
	}
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
		extendedSeq.AssignRead( lowFreqReads[i].read, lowFreqReads[i].strand, 0.95, assign ) ;
		if ( assign.seqIdx != -1 )
		{
			if ( assign.similarity < 1.0 )
				extendedSeq.AddAssignedRead( lowFreqReads[i].read, assign ) ;
		}
		//else if ( extendedSeq.AddRead( lowFreqReads[i].read, name, 1.0 ) == -1 )
		else
		{
			if ( lowFreqReads[i].minCnt < 2 )
			{
				lowFreqSeqSet.AssignRead( lowFreqReads[i].read, lowFreqReads[i].strand, 0.95, assign ) ;
				if ( assign.seqIdx != -1 )
				{
					lowFreqSeqSet.AddAssignedRead( lowFreqReads[i].read, assign ) ;
				}
			}
			else if ( lowFreqSeqSet.AddRead( lowFreqReads[i].read, name, 0, 0.97 ) < 0 )
			{
				struct _overlap geneOverlap[4] ;
				//buffer[0] = '\0' ;
				refSet.AnnotateRead( lowFreqReads[i].read, 0, geneOverlap, NULL, buffer ) ;

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
							lowFreqReads[i].read, geneOverlap[j].strand ) ;
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
	
	extendedSeq.Output( fp ) ;
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

	PrintLog( "Finished." ) ;
	return 0 ;
}
