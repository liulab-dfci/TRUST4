#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

#include <vector>

#include "KmerCount.hpp"
#include "SeqSet.hpp"
#include "AlignAlgo.hpp"

char usage[] = "./annotator [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the receptor genome sequence\n"
		"\t-a STRING: path to the assembly file\n"
		"\t-r STRING: path to the reads used in the assembly\n"
		"Optional:\n"
		"\t--fasta: the assembly file is in fasta format (default: false)\n";

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

char buffer[10241] = "" ;
char buffer2[10241] = "" ;
char weightBuffer[100241] = "" ;
char seq[10241] = "" ;

static const char *short_options = "f:a:r:" ;
static struct option long_options[] = {
			{ "fasta", no_argument, 0, 10000 },
			{ "radius", required_argument, 0, 10001 },
			{ (char *)0, 0, 0, 0} 
			} ;

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int radius = 10 ;

	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	SeqSet refSet( 7 ) ;
	int c, option_index ;
	FILE *fpAssembly = NULL ;
	struct _overlap geneOverlap[4] ;
	struct _overlap cdr[3] ; // the coordinate for cdr1,2,3
	option_index = 0 ;
	bool ignoreWeight = false ;
	ReadFiles reads ;

	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if ( c == -1 )
			break ;
		
		if ( c == 'f' )
		{
			refSet.InputRefFa( optarg ) ;
		}
		else if ( c == 'a' )
		{
			fpAssembly = fopen( optarg, "r" ) ;
		}
		else if ( c == 'r' )
		{
			reads.AddReadFile( optarg, false ) ;
		}
		else if ( c == 10000 )
		{
			ignoreWeight = true ;
		}
		else if ( c == 10001 )
		{
			radius = atoi( optarg ) ;
		}
		else
		{
			fprintf( stderr, "%s", usage ) ;
			return EXIT_FAILURE ;
		}
	}

	if ( refSet.Size() == 0 )
	{
		fprintf( stderr, "Need to use -f to specify the receptor genome sequence.\n" ) ;
		return EXIT_FAILURE ;
	}
	
	// Obtain the kmer count
	int kmerLength = 31 ;
	KmerCount kmerCount( kmerLength ) ;	
	while ( reads.Next() )
		kmerCount.AddCount( reads.seq ) ;


	if ( fpAssembly == NULL )
	{
		fprintf( stderr, "Need to use -a to specify the assembly file.\n" ) ;
		return EXIT_FAILURE ;
	}
	refSet.SetRadius( radius ) ;
	while ( fgets( buffer, sizeof( buffer ), fpAssembly ) != NULL )
	{
		if ( buffer[0] != '>' )
		{
			printf( "%s", buffer ) ;
			continue ;
		}

		fgets( seq, sizeof( seq ), fpAssembly ) ;
		
		// Read in the four line of pos weight
		SimpleVector<struct _posWeight> posWeight ;
		double depthSum = 0 ;
		if ( !ignoreWeight )
		{
			posWeight.ExpandTo( strlen( seq ) ) ;
			for ( k = 0 ; k < 4 ; ++k )
			{
				fgets( weightBuffer, sizeof( weightBuffer ), fpAssembly ) ;

				int num = 0 ;
				i = 0 ;
				for ( j = 0 ; weightBuffer[j] && weightBuffer[j] != '\n' ; ++j )
				{
					if ( weightBuffer[j] == ' ' )
					{
						posWeight[i].count[k] = num ;
						depthSum += num ;
						++i ;
						num = 0 ;
					}
					else
						num = num * 10 + weightBuffer[j] - '0' ; 
				}
			}
		}

		for ( i = 0 ; buffer[i] && buffer[i] != '\n' && buffer[i] != ' ' ; ++i )
			;
		if ( buffer[i] != ' ' )
			buffer[i] = ' ' ;
		buffer[i + 1] = '\0' ;
		int len = strlen( seq ) ;
		if ( seq[len - 1] == '\n' )
		{
			seq[len - 1] = '\0' ;
			--len ;
		}
		//TODO: adjust the average depth 
		sprintf( buffer + i + 1, "%d %.2lf", len, depthSum / len ) ;
		refSet.AnnotateRead( seq, 2, geneOverlap, cdr, buffer + i + 1 ) ;
		printf( "%s\n%s\n", buffer, seq ) ;

		// Extract the alternative sequence of CDR1,2,3.
		for ( i = 0 ; i < 3 && !ignoreWeight ; ++i )
		{
			if ( cdr[i].seqIdx == -1 )
				continue ;
			
			SimpleVector<struct _pair> alterNuc ;
			alterNuc.Reserve( cdr[i].readEnd - cdr[i].readStart + 1 ) ;
			for ( j = cdr[i].readStart ; j <= cdr[i].readEnd ; ++j )
			{
				int sum = 0 ;
				int max = 0 ;
				for ( k = 0 ; k < 4 ; ++k )
				{
					sum += posWeight[j].count[k] ;
					if ( posWeight[j].count[k] > max )
						max = posWeight[j].count[k] ;
				}

				for ( k = 0 ; k < 4 ; ++k )
				{
					if ( numToNuc[k] != seq[j] && 
						( sum <= posWeight[j].count[k] * 3 
							|| ( sum >= 20 && sum <= posWeight[j].count[k] * 4 ) ) )
					{
						struct _pair na ;
						na.a = j ;
						na.b = k ;
						alterNuc.PushBack( na ) ;
					}
				}
			}

			// Naively build the alternative CDR3.
			// So far, only handle one difference case.
			if ( alterNuc.Size() > 0 ) //#&& alterNuc.Size() <= 2 )
			{
				int size = alterNuc.Size() ;
				char *alterCDR = ( char * )malloc( sizeof( char ) * ( cdr[i].readEnd - cdr[i].readStart + 2 ) ) ;
				int cdrLen = cdr[i].readEnd - cdr[i].readStart + 1 ;
				alterCDR[ cdrLen ] = '\0' ;
				int id = 1 ;
				int n, l ;
				for ( n = 1 ; n < (1<<size) ; ++n )
				{
					memcpy( alterCDR, seq + cdr[i].readStart, cdrLen ) ;
					for ( l = 0 ; l < size ; ++l )
					{
						if ( n & ( 1 << l ) )
							alterCDR[ alterNuc[l].a - cdr[i].readStart ] = numToNuc[ alterNuc[l].b ] ; 	
					}
					//sprintf( buffer + strlen( buffer ), " minorCDR%d=%s", i + 1, alterCDR ) ;
					// Check whether this CDR3 is in the read or not. 
					if ( cdrLen < kmerLength )
					{
						if ( cdr[i].readEnd + kmerLength - cdrLen < len )
						{
							if ( kmerCount.GetCount( seq + cdr[i].readStart ) <= 0)
								continue ;
						}
						else
						{
							if ( kmerCount.GetCount( seq + cdr[i].readStart - ( kmerLength - cdrLen ) ) <= 0 )
								continue ;
						}
					}
					else
					{
						bool fail = false ; 
						// Naive implementation
						for ( l = cdr[i].readStart ; l + kmerLength - 1 <= cdr[i].readEnd ; ++l )
							if ( kmerCount.GetCount( seq + l ) <= 0 )
								fail = true ;
						if ( fail )
							continue ;
					}

					// Output the new seq
					k = 0 ;
					for ( j = 0 ; buffer[j] != ' ' ; ++j, ++k )
					{
						buffer2[k] = buffer[j] ;
					}
					buffer2[k] = '.' ;
					sprintf( buffer2 + k + 1, "%d", id ) ;
					++id ;
					while ( buffer2[k] )
						++k ;
					
					for ( ; buffer[j] != '=' ; ++j, ++k )
						buffer2[k] = buffer[j] ;
					buffer2[k] = '=' ;
					++k ;
					strcpy( buffer2 + k, alterCDR ) ;
					printf( "%s\n%s\n", buffer2, seq ) ;
				}
				free( alterCDR ) ;
			}
		}
		
	}
	fclose( fpAssembly ) ;

	return 0 ;
}
