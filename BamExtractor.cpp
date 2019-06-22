#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

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
		"\t-u: filter unaligned read-pair (default: no filter)\n" ;

static const char *short_options = "f:b:o:u" ;
static struct option long_options[] = {
			{ (char *)0, 0, 0, 0} 
			} ;

char buffer[100001] ;
char buffer2[100001] ;
char bufferQual[100001] ;
char bufferQual2[100001] ;
char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;


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

	char *qual1 ;
	char *qual2 ;
} ;

bool ValidAlternativeChrom( char *chrom )
{
	if ( strstr( chrom, "_random" ) || strstr( chrom, "_alt" ) )
	{
		if ( chrom[0] == 'c' && chrom[3] >= '0' && chrom[3] <= '9' )
			return true ;
	}
	else if ( ( chrom[0] == 'G' && chrom[1] == 'I' ) || ( chrom[0] == 'G' && chrom[1] == 'L' ) )
		return true ;
	return false ;
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

void TrimName( std::string &name )
{
	int len = name.length() ;
	if ( ( name[len - 1] == '1' || name[len - 1] == '2' ) 
			&& name[len - 2] == '/' )
	{
		name.erase( len - 2, 2 ) ; // Haven't tested yet.
	}
}

void OutputSeq( FILE *fp, const char *name, char *seq, char *qual )
{
	fprintf( fp, "@%s\n%s\n+\n%s\n", name, seq, qual ) ;
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
	bool filterUnalignedFragment = false ;
	SeqSet refSet( 21 ) ;

	std::map<std::string, struct _candidate> candidates ; 

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
			filterUnalignedFragment = true ;
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
	FILE *fp2 ;
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

	// assuming the input is sorted by coordinate.
	while ( alignments.Next() )
	{
		// If not aligned, output it.
		//if ( !alignments.IsPrimary() )
		//	continue ;
		alignments.GetReadSeq( buffer ) ;
		alignments.GetQual( bufferQual ) ;
		//if ( IsLowComplexity( buffer ) )
		//	continue ;

		if ( !alignments.IsTemplateAligned() 
			|| ( alignments.IsAligned() && ValidAlternativeChrom( alignments.GetChromName( alignments.GetChromId() ) ) ) )
		{
			if ( filterUnalignedFragment && !alignments.IsTemplateAligned() )
				continue ;
			
			//if ( !alignments.IsTemplateAligned() && !refSet.HasHitInSet( buffer ) ) 
			//	continue ;
			
			if ( !alignments.IsTemplateAligned() && alignments.fragStdev != 0 )
			{
				//printf( "filtered\n" ) ;
				std::string name( alignments.GetReadId() ) ;
				strcpy( buffer2, buffer ) ;
				alignments.GetQual( bufferQual2 ) ;

				if ( !alignments.Next() )
				{
					fprintf( stderr, "Two reads from the unaligned fragment are not showing up together. Please use -u option.\n") ;
					return EXIT_FAILURE ;
				}
				std::string mateName( alignments.GetReadId() ) ;
				alignments.GetReadSeq( buffer ) ;
				alignments.GetQual( bufferQual ) ;

				if ( name.compare( mateName ) != 0 )
				{
					fprintf( stderr, "Two reads from the unaligned fragment are not showing up together. Please use -u option.\n") ;
					return EXIT_FAILURE ;
				}

				if ( ( refSet.HasHitInSet( buffer2 ) || refSet.HasHitInSet( buffer ) ) 
					&& ( !IsLowComplexity( buffer2 ) && !IsLowComplexity( buffer ) ) )
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
				}
				continue ;
			}
			
			//printf( "%s %s\n", alignments.GetChromName( alignments.GetChromId() ), alignments.GetReadId() ) ;
			if ( alignments.fragStdev != 0 )
			{
				std::string name( alignments.GetReadId() ) ;
				TrimName( name ) ;	

				if ( candidates.find( name ) == candidates.end() )
				{
					candidates[name].mate1 = NULL ;
					candidates[name].mate2 = NULL ;
				}
			}
			else if ( !IsLowComplexity( buffer ) && refSet.HasHitInSet( buffer ) )
			{
				//alignments.GetReadSeq( buffer ) ;
				OutputSeq( fp1, alignments.GetReadId(), buffer, bufferQual ) ;
			}
			continue ;
		}

		if ( !alignments.IsAligned() ) // when reach here, it is parie-end case, and the other mate is aligned.
			continue ;
		
		if ( IsLowComplexity( buffer ) )
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
			TrimName( name ) ;

			if ( candidates.find( name ) == candidates.end() )
			{
				candidates[name].mate1 = NULL ;
				candidates[name].mate2 = NULL ;
			}
		}
		else
		{
			//alignments.GetReadSeq( buffer ) ;
			OutputSeq( fp1, alignments.GetReadId(), buffer, bufferQual ) ;
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
	fprintf( stderr, "Finish obtaining the candidate read ids.\n" ) ;


	while ( alignments.Next() )
	{
		if ( !alignments.IsPrimary() )
			continue ;
		if ( !alignments.IsTemplateAligned() ) // the sorted bam file should put all unaligned template at last.
			break ;

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
		
		if ( it->second.mate1 != NULL && it->second.mate2 != NULL )
		{
			OutputSeq( fp1, name.c_str(), it->second.mate1, it->second.qual1 ) ;
			OutputSeq( fp2, name.c_str(), it->second.mate2, it->second.qual2 ) ;
			free( it->second.mate1 ) ;
			free( it->second.mate2 ) ;
			free( it->second.qual1 ) ;
			free( it->second.qual2 ) ;
		}
			
	}
	fclose( fp1 ) ;
	fclose( fp2 ) ;

	fclose( fpRef ) ;
	return 0 ;
}
