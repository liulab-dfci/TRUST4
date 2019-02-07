#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

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
	SeqSet seqSet( 9 ) ;
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
			seqSet.InputRefFa( optarg ) ;
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
	
	if ( seqSet.Size() == 0 )
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
			continue ;

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
			fprintf( stderr, "Read in and count kmers for %d reads.\n", i ) ;
	}
	
	while ( mateReads.Next() )
	{
		/*struct _overlap geneOverlap[4] ;
		refSet.AnnotateRead( mateReads.seq, 0, geneOverlap, buffer ) ;
		
		if ( geneOverlap[0].seqIdx + geneOverlap[1].seqIdx + geneOverlap[2].seqIdx + geneOverlap[3].seqIdx == -4 )
			continue ;*/
		if ( IsLowComplexity( reads.seq ) )
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
			fprintf( stderr, "Read in and count kmers for %d reads.\n", i ) ;
	}

	fprintf( stderr, "Found %i reads.\n", i ) ;
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
#ifdef DEBUG
	printf( "Finish sorting\n" ) ;
#endif
	
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
				addRet = seqSet.AddRead( sortedReads[i].read ) ;
		}
		else
		{
			//printf( "saved time\n" ) ;
			if ( prevAddRet != -1 )
				addRet = seqSet.RepeatAddRead( sortedReads[i].read ) ;
		}
		if ( i > 0 && i % 100000 == 0 )
			fprintf( stderr, "Processed %d reads.\n", i ) ;
		
		if ( addRet == -2 )
			rescueReadIdx.push_back( i ) ;
		else if ( addRet >= 0 )
		{
			++assembledReadCnt ;
			assembledReadIdx.push_back( i ) ;
		}

		if ( assembledReadCnt > 0 && assembledReadCnt % 10000 == 0 )
			seqSet.UpdateAllConsensus() ;

		prevAddRet = addRet ;
#ifdef DEBUG
		printf( "done\n" ) ;
#endif
	}
	seqSet.UpdateAllConsensus() ;
	fprintf( stderr, "Assembled %d reads.\n", assembledReadCnt ) ;
	
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
	fprintf( stderr, "Try to rescue %d reads for assembly.\n", rescueReadCnt) ;
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
	fprintf( stderr, "Rescued %d reads.\n", assembledReadCnt ) ;
	
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

	std::vector<struct _Read> assembledReads ;
	assembledReadCnt = assembledReadIdx.size() ;
	for ( i = 0 ; i < assembledReadCnt ; ++i )
	{
		struct _Read nr ;
		nr.id = sortedReads[ assembledReadIdx[i] ].id ;
		nr.seq = sortedReads[ assembledReadIdx[i] ].read ;
		assembledReads.push_back( nr ) ;
	}
	std::sort( assembledReads.begin(), assembledReads.end(), CompSortReadById ) ;	
	SeqSet extendedSeq( 17 ) ;
	extendedSeq.InputSeqSet( seqSet, false ) ;
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
	
	
	fprintf( stderr, "Extend assemblies by mate pair information and their overlaps.\n" ) ;
	extendedSeq.ExtendSeqFromReads( assembledReads ) ;
	
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
	
	
	if ( outputPrefix[0] != '-' )
	{
		sprintf( buffer, "%s_final.out", outputPrefix ) ;
		fp = fopen( buffer, "w" ) ;
	}
	else
		fp = stdout ;
	
	i = extendedSeq.ExtendSeqFromSeqOverlap() ;
	extendedSeq.UpdateAllConsensus() ;
	
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
	return 0 ;
}
