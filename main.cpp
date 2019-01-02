#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

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
		"\t-b STRING: path to BAM alignment file\n" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;


static const char *short_options = "f:u:1:2:b:" ;
static struct option long_options[] = {
			{ (char *)0, 0, 0, 0} 
			} ;

int main( int argc, char *argv[] )
{
	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ;
	option_index = 0 ;
	SeqSet seqSet( 9 ) ;

	ReadFiles reads ;
	ReadFiles mateReads ;

	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			seqSet.InputRefFa( optarg ) ;
		}
		else if ( c == 'u' )
		{
			reads.AddReadFile( optarg ) ;
		}
		else if ( c == '1' )
		{
			reads.AddReadFile( optarg ) ;
		}
		else if ( c == '2' )
		{
			mateReads.AddReadFile( optarg ) ;
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
	
	/*while ( reads.Next() )
	{
		printf( "%s %s\n", reads.id, reads.seq ) ;
		fflush( stdout ) ;
		seqSet.AddRead( reads.seq ) ;
		printf( "done\n" ) ;
	}*/
	
	// Go through the second round.
	// TODO: user-defined number of rounds.
	seqSet.ResetPosWeight() ;
	reads.Rewind() ;
	while ( reads.Next() )
	{
		printf( "%s %s\n", reads.id, reads.seq ) ;
		fflush( stdout ) ;
		seqSet.AddRead( reads.seq ) ;
		printf( "done\n" ) ;
	}
	seqSet.Output() ;
	return 0 ;
}
