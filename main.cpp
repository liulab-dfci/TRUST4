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

char usage[] = "./bcr [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t[Read file]\n"
		"\t-u STRING: path to single-end read file\n"
		"\t\tor\n"
		"\t-1 STRING -2 STRING: path to paried-end read files\n"
		"\t\tor\n"
		"\t-b STRING: path to BAM alignment file\n"
		"Optional:\n"
		"\t-o STRING: prefix of the output file (default: batas)\n"
		"\t-c STRING: the path to the kmer count file\n" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

char buffer[10240] = "" ;

static const char *short_options = "f:u:1:2:b:o:c:" ;
static struct option long_options[] = {
			{ "debug-ns", required_argument, 0, 10000 },
			{ (char *)0, 0, 0, 0} 
			} ;

struct _sortRead
{
	char *id ;
	char *read ;
	int minCnt ;
	int medianCnt ;
	double avgCnt ;

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
	int i, j ;

	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ;
	option_index = 0 ;
	SeqSet seqSet( 9 ) ; // Only hold the novel seq.
	SeqSet refSet( 9 ) ;
	KmerCount kmerCount( 21 ) ;
	char outputPrefix[200] = "batas" ;

	ReadFiles reads ;
	ReadFiles mateReads ;
	bool countMyself = true ;
	int maxReadLen = -1 ;

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
		}
		else if ( c == '2' )
		{
			mateReads.AddReadFile( optarg, true ) ;
		}
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;
		}
		else if ( c == 'c' )
		{
			kmerCount.AddCountFromFile( optarg ) ;
			countMyself = false ;
			fprintf( stderr, "Read in the kmer count information from %s", optarg ) ;
		}
		else if ( c == 10000 )
		{
			seqSet.InputNovelFa( optarg ) ;
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

		if ( IsLowComplexity( reads.seq ) )
		{
			continue ;
		}
		struct _sortRead nr ;
		nr.read = strdup( reads.seq ) ;
		nr.id = strdup( reads.id ) ;

		sortedReads.push_back( nr ) ;
		if ( countMyself )
			kmerCount.AddCount( reads.seq ) ;

		++i ;

		int len = strlen( reads.seq ) ;
		if ( len > maxReadLen )
			maxReadLen = len ;
	
		if ( countMyself && i % 100000 == 0 )
			PrintLog( "Read in and count kmers for %d reads.", i ) ;
		else if ( !countMyself && i % 1000000 == 0 )
			PrintLog( "Read in %d reads.", i ) ;
	}
	
	while ( mateReads.Next() )
	{
		/*struct _overlap geneOverlap[4] ;
		refSet.AnnotateRead( mateReads.seq, 0, geneOverlap, buffer ) ;
		
		if ( geneOverlap[0].seqIdx + geneOverlap[1].seqIdx + geneOverlap[2].seqIdx + geneOverlap[3].seqIdx == -4 )
			continue ;*/
		if ( IsLowComplexity( mateReads.seq ) )
			continue ;

		struct _sortRead nr ;
		nr.read = strdup( mateReads.seq ) ;
		nr.id = strdup( mateReads.id ) ;

		sortedReads.push_back( nr ) ;
		if ( countMyself )
			kmerCount.AddCount( mateReads.seq ) ;

		++i ;
		
		int len = strlen( reads.seq ) ;
		if ( len > maxReadLen )
			maxReadLen = len ;
		
		if ( countMyself && i % 100000 == 0 )
			PrintLog( "Read in and count kmers for %d reads.", i ) ;
		else if ( !countMyself && i % 1000000 == 0 )
			PrintLog( "Read in %d reads.", i ) ;
	}

	PrintLog( "Found %i reads.", i ) ;
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
		kmerCount.GetCountStats( sortedReads[i].read, sortedReads[i].minCnt, sortedReads[i].medianCnt, sortedReads[i].avgCnt ) ;
	}

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
	for ( i = 0 ; i < readCnt ; ++i )
	{
#ifdef DEBUG
		printf( "%s %s %d %lf\n", sortedReads[i].id, sortedReads[i].read, sortedReads[i].minCnt, sortedReads[i].avgCnt ) ;
		fflush( stdout ) ;
#endif
		int addRet = -1 ;
		
		if ( i == 0 || strcmp( sortedReads[i].read, sortedReads[i - 1].read ) )
		{
			//printf( "new stuff\n" ) ;
			struct _overlap geneOverlap[4] ;
			refSet.AnnotateRead( sortedReads[i].read, 0, geneOverlap, buffer ) ;

			if ( geneOverlap[3].seqIdx != -1 && geneOverlap[3].seqStart >= 100 ) // From constant gene.
				addRet = -1 ;
			else
			{
				addRet = seqSet.AddRead( sortedReads[i].read ) ;
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
			}
		}
		else
		{
			//printf( "saved time\n" ) ;
			if ( prevAddRet != -1 )
				addRet = seqSet.RepeatAddRead( sortedReads[i].read ) ;
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

		if ( i > 0 && i % 100000 == 0 )
			PrintLog( "Processed %d reads (%d are used for assembly).", i, assembledReadCnt ) ;
		
		prevAddRet = addRet ;
#ifdef DEBUG
		printf( "done\n" ) ;
#endif
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
		addRet = seqSet.AddRead( sortedReads[ rescueReadIdx[i] ].read ) ;
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
	FILE *fp ;
	if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_raw_assembly.out", outputPrefix ) ;
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
		nr.overlap.seqIdx = -1 ;
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

	SeqSet extendedSeq( 17 ) ;
	extendedSeq.InputSeqSet( seqSet, false ) ;
	struct _overlap assign ;
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		if ( i == 0 || strcmp( assembledReads[i].read, assembledReads[i - 1].read ) )
			extendedSeq.AssignRead( assembledReads[i].read, 1.0, assign ) ;
		assembledReads[i].overlap = assign ;
			
		if ( i > 0 && i % 100000 == 0 )
		{
			PrintLog( "Processed %d reads for extension.", i ) ;
		}
	}
	
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
	extendedSeq.ExtendSeqFromReads( assembledReads, 31 ) ; //( avgReadLen / 30 < 31 ) ? 31 : ( avgReadLen / 3 )  ) ;
	
	if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_extended.out", outputPrefix ) ;
		fp = fopen( buffer, "w" ) ;
	}
	else
		fp = stdout ;

	extendedSeq.Output( fp ) ;
	fflush( fp ) ;
	
	if ( outputPrefix[0] != '-' )
		fclose( fp ) ;
	
	
	PrintLog( "Extend assemblies by their overlap." ) ;
	extendedSeq.ChangeKmerSize( 31 ) ;
	i = extendedSeq.ExtendSeqFromSeqOverlap( ( avgReadLen / 2 < 31 ) ? 31 : ( avgReadLen / 2 ) ) ;
	extendedSeq.UpdateAllConsensus() ;
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
	
	for ( i = 0 ; i < readCnt ; ++i )
	{
		free( sortedReads[i].id ) ;
		free( sortedReads[i].read ) ;
	}

	PrintLog( "Finished." ) ;
	return 0 ;
}
