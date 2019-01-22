#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include "alignments.hpp"

char usage[] = "./bam-extractor [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t-b STRING: path to BAM file\n"
		"\t-o STRING: prefix to the output file\n" ;

static const char *short_options = "f:b:o:" ;
static struct option long_options[] = {
			{ (char *)0, 0, 0, 0} 
			} ;

char buffer[100001] ;


struct _interval
{
	int chrId ;
	int start, end ;

	bool operator<( const struct _interval &b )
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
} ;

bool ValidAlternativeChrom( char *chrom )
{
	if ( strstr( chrom, "_random" ) || strstr( chrom, "_alt" ) )
	{
		if ( chrom[0] == 'c' && chrom[3] >= '0' && chrom[3] <= '9' )
			return true ;
	}
	return false ;
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
	char prefix[127] = "toassemble" ;
	Alignments alignments ;

	std::map<std::string, struct _candidate> candidates ; 

	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			fpRef = fopen( optarg, "r" ) ;
		}
		else if ( c == 'b' )
		{
			alignments.Open( optarg ) ;
		}
		else if ( c == 'o' )
		{
			strcpy( prefix, optarg ) ;
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
	
	int tag = 0 ;
	
	FILE *fp1 ;
	if ( alignments.fragStdev == 0 )
	{
		sprintf( buffer, "%s.fa", prefix ) ;
		fp1 = fopen( buffer, "w" ) ;
	}
	// assuming the input is sorted by coordinate.
	while ( alignments.Next())
	{
		// If not aligned, output it.
		//if ( !alignments.IsPrimary() )
		//	continue ;
		
		if ( !alignments.IsTemplateAligned() 
			|| ValidAlternativeChrom( alignments.GetChromName( alignments.GetChromId() ) ) )
		{
			if ( alignments.fragStdev != 0 )
			{
				std::string name( alignments.GetReadId() ) ;
				int len = name.length() ;
				if ( ( name[len - 1] == '1' || name[len - 1] == '2' ) 
						&& name[len - 2] == '/' )
				{
					name.erase( len - 2, 2 ) ; // Haven't tested yet.
				}

				if ( candidates.find( name ) == candidates.end() )
				{
					candidates[name].mate1 = NULL ;
					candidates[name].mate2 = NULL ;
				}
			}
			else
			{
				alignments.GetReadSeq( buffer ) ;
				fprintf( fp1, ">%s\n%s\n", alignments.GetReadId(), buffer ) ;
			}
			continue ;
		}

		if ( !alignments.IsAligned() ) // when reach here, it is parie-end case, and the other mate is aligned.
			continue ;

		// The aligned reads can reach here.
		int chrId = alignments.GetChromId() ;
		int start = (int)alignments.segments[0].a ;
		int end = (int)alignments.segments[ alignments.segCnt - 1 ].b ;
		//if ( !strcmp( "ERR188021.4488674", alignments.GetReadId() ) )
		//	printf( "hi %d %d %d: %d %d\n", chrId, start, end, tag, genes[tag].chrId ) ;
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
			int len = name.length() ;
			if ( ( name[len - 1] == '1' || name[len - 1] == '2' ) 
					&& name[len - 2] == '/' )
			{
				name.erase( len - 2, 2 ) ; // Haven't tested yet.
			}
			if ( candidates.find( name ) == candidates.end() )
			{
				candidates[name].mate1 = NULL ;
				candidates[name].mate2 = NULL ;
			}
		}
		else
		{
			alignments.GetReadSeq( buffer ) ;
			fprintf( fp1, ">%s\n%s\n", alignments.GetReadId(), buffer ) ;
		}
	}
	alignments.Rewind() ;
	if ( alignments.fragStdev == 0 ) // Single-end can terminate here.
	{
		fclose( fp1 ) ;
		alignments.Close() ;
		return 0 ;
	}
	// Case of pair-end data set.
	// Go through the BAM file again to output the candidates 
	FILE *fp2 ;

	sprintf( buffer, "%s_1.fa", prefix ) ;
	fp1 = fopen( buffer, "w" ) ;
	sprintf( buffer, "%s_2.fa", prefix ) ;
	fp2 = fopen( buffer, "w" ) ;

	while ( alignments.Next() )
	{
		if ( !alignments.IsPrimary() )
			continue ;
		std::string name( alignments.GetReadId() ) ;
		int len = name.length() ;
		if ( ( name[len - 1] == '1' || name[len - 1] == '2' ) 
				&& name[len - 2] == '/' )
		{
			name.erase( len - 2, 2 ) ; // Haven't tested yet.
		}
		
		std::map<std::string, struct _candidate>::iterator it = candidates.find( name ) ;
		if ( it == candidates.end() )
			continue ;

		alignments.GetReadSeq( buffer ) ;
		if ( alignments.IsFirstMate() )
			it->second.mate1 = strdup( buffer ) ;
		else
			it->second.mate2 = strdup( buffer ) ;

		if ( it->second.mate1 != NULL && it->second.mate2 != NULL )
		{
			fprintf( fp1, ">%s\n%s\n", name.c_str(), it->second.mate1 ) ;
			fprintf( fp2, ">%s\n%s\n", name.c_str(), it->second.mate2 ) ;
			free( it->second.mate1 ) ;
			free( it->second.mate2 ) ;
		}
			
	}
	fclose( fp1 ) ;
	fclose( fp2 ) ;

	return 0 ;
}
