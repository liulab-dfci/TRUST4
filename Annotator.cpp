#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

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
		"\t--fasta: the assembly file is in fasta format (default: false)\n"
		"\t-o STRING: the prefix of the file containing CDR3 information (default: trust)\n" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

char buffer[10241] = "" ;
char buffer2[10241] = "" ;
char weightBuffer[100241] = "" ;
char seq[10241] = "" ;

static const char *short_options = "f:a:r:p:" ;
static struct option long_options[] = {
			{ "fasta", no_argument, 0, 10000 },
			{ "radius", required_argument, 0, 10001 },
			{ (char *)0, 0, 0, 0} 
			} ;

struct _annotate
{
	struct _overlap geneOverlap[4] ;
	struct _overlap cdr[3] ;
} ;

struct _CDR3info
{
	char *seq ;
	int count ;
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
	SeqSet seqSet( 17 ) ;
	int c, option_index ;
	FILE *fpAssembly = NULL ;
	struct _overlap geneOverlap[4] ;
	struct _overlap cdr[3] ; // the coordinate for cdr1,2,3
	option_index = 0 ;
	bool ignoreWeight = false ;
	ReadFiles reads ;
	char outputPrefix[1024] = "trust" ;

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
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;
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
	
	if ( fpAssembly == NULL )
	{
		fprintf( stderr, "Need to use -a to specify the assembly file.\n" ) ;
		return EXIT_FAILURE ;
	}

	refSet.SetRadius( radius ) ;
	PrintLog( "Annotate assemblies." ) ;
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
		buffer[i] = '\0' ;
		
		int len = strlen( seq ) ;
		if ( seq[len - 1] == '\n' )
		{
			seq[len - 1] = '\0' ;
			--len ;
		}

		if ( !ignoreWeight )
			seqSet.InputNovelSeq( buffer, seq, posWeight ) ;
		else 
			seqSet.InputNovelRead( buffer, seq, 1 ) ;
	}
	fclose( fpAssembly ) ;

	int seqCnt = seqSet.Size() ;
	struct _annotate *annotations = new struct _annotate[ seqCnt ] ;
	for ( i = 0 ; i < seqCnt ; ++i )
	{
		int weightSum = seqSet.GetSeqWeightSum( i ) ; 
		int len = seqSet.GetSeqConsensusLen( i ) ;
		sprintf( buffer, ">%s %d %.2lf", seqSet.GetSeqName( i ), len, (double)weightSum / 500.0 ) ;
		refSet.AnnotateRead( seqSet.GetSeqConsensus( i ), 2, annotations[i].geneOverlap, annotations[i].cdr, buffer + strlen( buffer ) ) ;
		printf( "%s\n%s\n", buffer, seqSet.GetSeqConsensus( i ) ) ;
	}
	
	// Output more CDR3 information 
	if ( reads.GetFpUsed() > 0 )
	{
		struct _overlap assign ;
		std::vector< std::vector<struct _CDR3info> > cdr3Infos ;
		cdr3Infos.resize( seqCnt ) ;
		k = 0 ;
		buffer2[0] = '\0' ;
		PrintLog( "Start to realign reads." ) ;
		while ( reads.Next() )	
		{
			if ( strcmp( reads.seq, buffer2 ) ) 
			{
				strcpy( buffer2, reads.seq ) ;
				seqSet.AssignRead( reads.seq, 0, 0.95, assign ) ;	
				if ( assign.seqIdx == -1 )
					continue ;
			}

			if ( assign.seqIdx == -1 )
				continue ;

			if ( annotations[assign.seqIdx].cdr[2].seqIdx != -1 
					&& assign.seqStart <= annotations[ assign.seqIdx ].cdr[2].readStart 
					&& assign.seqEnd >=  annotations[ assign.seqIdx ].cdr[2].readEnd )
			{
				std::vector<struct _CDR3info> &info = cdr3Infos[ assign.seqIdx ] ;
				int size = info.size() ;
				int cdr3Len =  annotations[ assign.seqIdx ].cdr[2].readEnd -
					annotations[ assign.seqIdx ].cdr[2].readStart + 1 ;
				memcpy( buffer, reads.seq + assign.readStart +
						annotations[ assign.seqIdx ].cdr[2].readStart - assign.seqStart, 
						sizeof( char ) * cdr3Len ) ;
				buffer[cdr3Len] = '\0' ;

				for ( i = 0 ; i < size ; ++i )
					if ( !strcmp( info[i].seq, buffer ) )
					{
						++info[i].count ;
						break ;
					}
				if ( i >= size )
				{
					struct _CDR3info nc ;
					nc.seq = strdup( buffer ) ;
					nc.count = 1 ;

					info.push_back( nc ) ;
				}
			}
			++k ;
			if ( k % 100000 == 0 )
			{
				PrintLog( "Realigned %d reads.", k ) ;
			}
		}

		// Output different CDR3s from each main assembly.
		sprintf( buffer, "%s_cdr3.out", outputPrefix ) ;
		FILE *fpOutput = fopen( buffer, "w" ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( annotations[i].cdr[2].seqIdx == -1 )
				continue ;

			std::vector<struct _CDR3info> &info = cdr3Infos[i] ;
			int size = info.size() ;
			if ( size == 0 )
			{
				struct _CDR3info nc ;
				int len =  annotations[i].cdr[2].readEnd - annotations[i].cdr[2].readStart + 1 ;
				nc.seq = (char *)malloc( sizeof( char) * ( len + 1 ) ) ;
				memcpy( nc.seq, seqSet.GetSeqConsensus( i ) + annotations[i].cdr[2].readStart, len ) ;
				nc.seq[ len ] = '\0' ;
				nc.count = 1 ;

				info.push_back( nc ) ;
				size = 1 ;
			}
			// Output each CDR3 information
			for ( j = 0 ; j < size ; ++j )
			{
				fprintf( fpOutput, "%s\t%d\t", seqSet.GetSeqName( i ), j ) ;
				// The gene ids
				for ( k = 0 ; k < 4 ; ++k )
				{
					if ( k == 1 )
						continue ;
					if ( annotations[i].geneOverlap[k].seqIdx == -1 )
						fprintf( fpOutput, "*\t" ) ;	
					else
						fprintf( fpOutput, "%s\t", 
							refSet.GetSeqName( annotations[i].geneOverlap[k].seqIdx ) ) ;
				}
				// Output CDR1,2
				for ( k = 0 ; k < 2 ; ++k )
				{
					if ( annotations[i].cdr[k].seqIdx == -1 )
						fprintf( fpOutput, "*\t" ) ;
					else
					{
						int len = annotations[i].cdr[k].readEnd - annotations[i].cdr[k].readStart + 1 ;
						memcpy( buffer, seqSet.GetSeqConsensus( i ) + annotations[i].cdr[k].readStart, len ) ;
						buffer[len] = '\0' ;
						fprintf( fpOutput, "%s\t", buffer ) ;
					}
				}
				fprintf( fpOutput, "%s\t%d\n", info[j].seq, info[j].count ) ;
			}
			
		}
		// Free up memory
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = cdr3Infos[i].size() ;
			for ( j = 0 ; j < size; ++j )
				free( cdr3Infos[i][j].seq ) ;
		}
		fclose( fpOutput ) ;
	}
	
	delete[] annotations ;
	return 0 ;
}
