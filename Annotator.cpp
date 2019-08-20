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
		"\t-t INT: number of threads (default: 1)\n"
		"\t-o STRING: the prefix of the file containing CDR3 information (default: trust)\n"
		//"\t--partial: including partial CDR3s in the report (default: false)\n"
		"\t--noImpute: do not impute CDR3 sequence for TCR (default: not set (impute))"
		"\t--notIMGT: the receptor genome sequence is not in IMGT format (default: not set(in IMGT format))\n";

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

char buffer[10241] = "" ;
char weightBuffer[100241] = "" ;
char seq[10241] = "" ;

static const char *short_options = "f:a:r:p:o:t:" ;
static struct option long_options[] = {
			{ "fasta", no_argument, 0, 10000 },
			{ "radius", required_argument, 0, 10001 },
			{ "partial",no_argument, 0, 10002 },
			{ "notIMGT", no_argument, 0, 10003 },
			{ "noImpute", no_argument, 0, 10004 },
			{ (char *)0, 0, 0, 0} 
			} ;

struct _annotate
{
	struct _overlap geneOverlap[4] ;
	struct _overlap cdr[3] ;
	std::vector<struct _overlap> secondaryGeneOverlaps ;
} ;

struct _CDR3info
{
	char *seq ;
	double count ;
	double TPM ;
} ;

struct _annotateReadsThreadArg
{
	int tid ;
	int threadCnt ;
	bool impute ;

	SeqSet *seqSet ;
	SeqSet *refSet ;
	struct _annotate *annotations ;
} ;

struct _assignReadsThreadArg
{
	int tid ;
	int threadCnt ;

	SeqSet *seqSet ;
	std::vector<struct _assignRead> *pAssembledReads ;
	int assembledReadCnt ;
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

bool CompSortAssignedReadsByAssign( const struct _assignRead &a, const struct _assignRead &b )
{
	if ( a.overlap.seqIdx != b.overlap.seqIdx )
		return a.overlap.seqIdx < b.overlap.seqIdx ;
	else
		return strcmp( a.id, b.id ) < 0 ;
}

// We already make sure they overlap before calling this function
bool IsCDR3Compatible( const struct _assignRead &r, const struct _CDR3info &cdr3, const struct _overlap cdr3Coord )
{
	int rOffset ;
	int cOffset ;
	if ( r.overlap.seqStart <= cdr3Coord.readStart )
	{
		rOffset = r.overlap.readStart + cdr3Coord.readStart - r.overlap.seqStart ;
		cOffset = 0 ;
	}
	else if ( r.overlap.seqStart > cdr3Coord.readStart ) 
	{
		rOffset = r.overlap.readStart ;
		cOffset = r.overlap.seqStart - cdr3Coord.readStart ;
	}

	int i ;
	for ( i = 0 ; r.read[i + rOffset] && cdr3.seq[i + cOffset] ; ++i )
	{
		if ( r.read[i + rOffset] != cdr3.seq[i + cOffset] )
			return false ;
	}
	//printf( "%s\n%s\n", r.read, cdr3.seq ) ;
	//printf( "%d %d %d\n", i, rOffset, cOffset ) ;
	return true ;
}


// Could be improved by using BitTable and compress the read covering the same set of CDR3s.
void AbundanceEstimation( std::vector< SimpleVector<int> > &compat, std::vector< struct _CDR3info >& info )
{
	int i, j, k, t ;
	double endD = 1e-6 ;
	int cCnt = info.size() ;
	int rCnt = compat.size() ;
	double sum = 0 ;

	SimpleVector<double> abundance ;
	abundance.ExpandTo( cCnt ) ;
	abundance.SetZero( 0, cCnt ) ;

	for ( i = 0 ; i < rCnt ; ++i )
		if ( compat[i].Size() == 1 )
			++abundance[ compat[i][0] ] ;
	for ( i = 0 ; i < cCnt ; ++i )
		sum += abundance[i] ;
	for ( i = 0 ; i < cCnt ; ++i )
		abundance[i] /= sum ;

	for ( t = 0 ; t < 1000 ; ++t )
	{
		double d = 0 ;
		for ( i = 0 ; i < cCnt ; ++i )
			info[i].count = 0 ;
		// E-step
		for ( i = 0 ; i < rCnt ; ++i )
		{
			int size = compat[i].Size() ;
			sum = 0 ;
			for ( j = 0 ; j < size ; ++j )
			{
				//printf( "%d %d %d\n", i, j, compat[i][j] ) ;
				sum += abundance[ compat[i][j] ] ;
			}
			for ( j = 0 ; j < size ; ++j )
				info[ compat[i][j] ].count += abundance[ compat[i][j] ] / sum ;
		}

		// M-step
		sum = 0 ;
		for ( i = 0 ; i < cCnt ; ++i )
			sum += info[i].count ;
		for ( i = 0 ; i < cCnt ; ++i )
		{
			double tmp = abundance[i] ;
			abundance[i] = info[i].count / sum ;
			double diff = ABS( tmp - abundance[i] ) ;
			if ( diff > d )
				d = diff ;
		}
		//printf( "%lf %lf %lf\n", abundance[0], abundance[1], d ) ;
		if ( d < endD )
			return ;
	}
}

void AnnotationTieBreak( struct _annotate *annotations, SeqSet &seqSet, SeqSet &refSet )
{
	int i, j, k ;
	int seqCnt = seqSet.Size() ;
	int refCnt = refSet.Size() ;

	double *abundance = new double[ refCnt ] ; 
	memset( abundance, 0, sizeof( double ) * refCnt ) ;
	// Collect the abundance
	for ( i = 0 ; i < seqCnt ; ++i )
	{
		int weightSum = seqSet.GetSeqWeightSum( i ) ;
		int len = seqSet.GetSeqConsensusLen( i ) ;
		double avgWeight = (double)weightSum / len ;
	
		for ( k = 0 ; k < 4 ; ++k )
		{
			if ( annotations[i].geneOverlap[k].seqIdx == -1 )
				continue ;
			abundance[ annotations[i].geneOverlap[k].seqIdx ] += avgWeight ;
		}		
	}
	
	// Break ties
	for ( i = 0 ; i < seqCnt ; ++i )
	{
		struct _overlap *geneOverlap = annotations[i].geneOverlap ;
		for ( k = 0 ; k < 4 ; ++k )
		{
			if ( geneOverlap[k].seqIdx == -1 )
				continue ;

			std::vector<struct _overlap> &overlaps = annotations[i].secondaryGeneOverlaps ;
			int size = overlaps.size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( refSet.GetGeneType( refSet.GetSeqName( overlaps[j].seqIdx ) ) == j )
					break ;
			}
			if ( j >= size )
				continue ;
			if ( overlaps[j].readStart == geneOverlap[k].readStart 
				&& overlaps[j].readEnd == geneOverlap[k].readEnd 
				&& overlaps[j].similarity == geneOverlap[k].similarity 
				&& abundance[ overlaps[j].seqIdx ] > abundance[ geneOverlap[k].seqIdx ] ) 
			{
				struct _overlap tmp ;
				tmp = geneOverlap[k] ;
				geneOverlap[k] = overlaps[j] ;
				overlaps[j] = tmp ;
			}
		}
	}
	delete[] abundance ;
}

void *AnnotateReads_Thread( void *pArg )
{
	struct _annotateReadsThreadArg &arg = *( (struct _annotateReadsThreadArg *)pArg ) ;
	int start, end ;
	int i ;
	int seqCnt = arg.seqSet->Size() ;
	char *buffer = new char[10001] ;
	start = seqCnt / arg.threadCnt * arg.tid ;
	end = seqCnt / arg.threadCnt * ( arg.tid + 1 ) ;
	if ( arg.tid == arg.threadCnt - 1 )
		end = seqCnt ;
	for ( i = start ; i < end ; ++i )
	{
		arg.refSet->AnnotateRead( arg.seqSet->GetSeqConsensus( i ), 2, arg.annotations[i].geneOverlap, arg.annotations[i].cdr, 
			&( arg.annotations[i].secondaryGeneOverlaps ) ) ;
		if ( arg.impute && arg.refSet->ImputeCDR3( arg.seqSet->GetSeqConsensus( i ), buffer, 
			arg.annotations[i].geneOverlap, arg.annotations[i].cdr, 
			&( arg.annotations[i].secondaryGeneOverlaps ) ) != -1 )
			arg.seqSet->SetSeqConsensus( i, buffer ) ;
	}
	delete[] buffer ;
	pthread_exit( NULL ) ;
}

void *AssignReads_Thread( void *pArg )
{
	struct _assignReadsThreadArg &arg = *( (struct _assignReadsThreadArg *)pArg ) ;
	int start, end ;
	int i ;
	std::vector< struct _assignRead> &assembledReads = *arg.pAssembledReads ;
	start = arg.assembledReadCnt / arg.threadCnt * arg.tid ;
	end = arg.assembledReadCnt / arg.threadCnt * ( arg.tid + 1 ) ;
	if ( arg.tid == arg.threadCnt - 1 )
		end = arg.assembledReadCnt ;
	struct _overlap assign ;
	for ( i = start ; i < end ; ++i )
	{
		if ( i == start || strcmp( assembledReads[i].read, assembledReads[i - 1].read ) )
			arg.seqSet->AssignRead( assembledReads[i].read, assembledReads[i].overlap.strand, assign ) ;
		assembledReads[i].overlap = assign ;	
	}
	pthread_exit( NULL ) ;
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
	char outputPrefix[1024] = "trust" ;
	FILE *fpReads = NULL ;
	int threadCnt = 1 ;
	bool includePartial = true ;
	bool isIMGT = true ;
	bool impute = true ;

	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if ( c == -1 )
			break ;
		
		if ( c == 'f' )
		{
			//refSet.InputRefFa( optarg ) ;
			strcpy( buffer, optarg ) ;
		}
		else if ( c == 'a' )
		{
			fpAssembly = fopen( optarg, "r" ) ;
		}
		else if ( c == 'r' )
		{
			fpReads = fopen( optarg, "r" ) ;
		}
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;
		}
		else if ( c == 't' )
		{
			threadCnt = atoi( optarg ) ;
		}
		else if ( c == 10000 )
		{
			ignoreWeight = true ;
		}
		else if ( c == 10001 )
		{
			radius = atoi( optarg ) ;
		}
		else if ( c == 10002 )
		{
			includePartial = true ;	
		}
		else if ( c == 10003 )
		{
			isIMGT = false ;
		}
		else if ( c == 10004 )
		{
			impute = false ;
		}
		else
		{
			fprintf( stderr, "%s", usage ) ;
			return EXIT_FAILURE ;
		}
	}

	refSet.InputRefFa( buffer, isIMGT ) ;

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
	
	refSet.SetHitLenRequired( 17 ) ;
	refSet.SetRadius( radius ) ;
	PrintLog( "Start to annotate assemblies." ) ;
	while ( fgets( buffer, sizeof( buffer ), fpAssembly ) != NULL )
	{
		if ( buffer[0] != '>' )
		{
			printf( "%s", buffer ) ;
			continue ;
		}

		fgets( seq, sizeof( seq ), fpAssembly ) ;
		int len = strlen( seq ) ;
		if ( seq[len - 1] == '\n' )
		{
			seq[len - 1] = '\0' ;
			--len ;
		}

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
				posWeight[i].count[k] = num ;
				depthSum += num ;
			}
		}
		for ( i = 0 ; buffer[i] && buffer[i] != '\n' && buffer[i] != ' ' ; ++i )
			;
		buffer[i] = '\0' ;
		
		if ( !ignoreWeight )
		{
			seqSet.InputNovelSeq( buffer + 1, seq, posWeight ) ;
		}
		else
			seqSet.InputNovelRead( buffer + 1, seq, 1 ) ;
	}
	fclose( fpAssembly ) ;

	int seqCnt = seqSet.Size() ;
	struct _annotate *annotations = new struct _annotate[ seqCnt ] ;
	if ( threadCnt <= 1 )
	{
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			refSet.AnnotateRead( seqSet.GetSeqConsensus( i ), 2, annotations[i].geneOverlap, annotations[i].cdr, 
					&annotations[i].secondaryGeneOverlaps ) ;
			
			if ( impute && refSet.ImputeCDR3( seqSet.GetSeqConsensus( i ), buffer, annotations[i].geneOverlap, annotations[i].cdr, 
					&annotations[i].secondaryGeneOverlaps ) != -1 )
			{
				seqSet.SetSeqConsensus( i, buffer ) ;
			}
		}
	}
	else
	{
		pthread_t *threads = new pthread_t[ threadCnt ] ;
		struct _annotateReadsThreadArg *args = new struct _annotateReadsThreadArg[threadCnt] ;
		pthread_attr_t attr ;

		pthread_attr_init( &attr ) ;
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;

		for ( i = 0 ; i < threadCnt ; ++i )
		{
			args[i].tid = i ;
			args[i].threadCnt = threadCnt ;
			args[i].seqSet = &seqSet ;
			args[i].refSet = &refSet ;
			args[i].annotations = annotations ;
			args[i].impute = impute ;
			pthread_create( &threads[i], &attr, AnnotateReads_Thread, (void *)( args + i ) ) ;
		}

		for ( i = 0 ; i < threadCnt ; ++i )
			pthread_join( threads[i], NULL ) ;


		delete[] threads ;
		delete[] args ;
	}
	// Use global information to break ties
	AnnotationTieBreak( annotations, seqSet, refSet ) ;

	// Output the annotation of consensus assemblies
	for ( i = 0 ; i < seqCnt ; ++i )
	{
		int weightSum = seqSet.GetSeqWeightSum( i ) ; 
		int len = seqSet.GetSeqConsensusLen( i ) ;
		sprintf( buffer, ">%s %d %.2lf", seqSet.GetSeqName( i ), len, (double)weightSum / 500.0 ) ;
		refSet.AnnotationToString( seqSet.GetSeqConsensus( i ), annotations[i].geneOverlap, 
			annotations[i].cdr, &annotations[i].secondaryGeneOverlaps, buffer + strlen( buffer ) ) ;
		printf( "%s\n%s\n", buffer, seqSet.GetSeqConsensus( i ) ) ;
	}

	
	// Output more CDR3 information 
	if ( fpReads != NULL )
	{
		struct _overlap assign ;
		std::vector< std::vector<struct _CDR3info> > cdr3Infos ;
		cdr3Infos.resize( seqCnt ) ;
		int strand, minCnt, medCnt ;
		k = 0 ;
		PrintLog( "Start to realign reads for CDR3 analysis." ) ;
		seqSet.Clean( false ) ;

		std::vector<struct _assignRead> cdr3Reads ; // Keep the information of the reads aligned to cdr3 region.
		std::vector<struct _assignRead> assembledReads ;
		int assembledReadCnt ;

		while ( fscanf( fpReads, "%s %d %d %d", buffer, &strand, &minCnt, &medCnt ) != EOF )	
		{
			fscanf( fpReads, "%s", seq ) ; 
			
			struct _assignRead nr ;
			nr.id = strdup( buffer + 1 ) ; // skip the > sign
			nr.read = strdup( seq ) ;
			nr.overlap.strand = strand ;
			assembledReads.push_back( nr ) ;
		}
		assembledReadCnt = assembledReads.size() ;
		
		if ( threadCnt <= 1 )
		{
			for ( i = 0 ; i < assembledReadCnt ; ++i )
			{
				if ( i == 0 || strcmp( assembledReads[i].read, assembledReads[i - 1].read ) ) 
					seqSet.AssignRead( assembledReads[i].read, assembledReads[i].overlap.strand, assign ) ;	
				assembledReads[i].overlap = assign ;	
				if ( ( i + 1 ) % 100000 == 0 )
				{
					PrintLog( "Realigned %d reads.", i + 1 ) ;
				}

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
				args[i].seqSet = &seqSet ;
				args[i].pAssembledReads = &assembledReads ;
				args[i].assembledReadCnt = assembledReadCnt ;
				pthread_create( &threads[i], &attr, AssignReads_Thread, (void *)( args + i ) ) ;
			}

			for ( i = 0 ; i < threadCnt ; ++i )
				pthread_join( threads[i], NULL ) ;


			delete[] threads ;
			delete[] args ;
		}
		for ( i = 0 ; i < assembledReadCnt ; ++i )
		{
			assign = assembledReads[i].overlap ;
			if ( assign.seqIdx == -1 )
				continue ;
					
			int cdr3Len =  annotations[ assign.seqIdx ].cdr[2].readEnd -
				annotations[ assign.seqIdx ].cdr[2].readStart + 1 ;
			
			// Store the read overlap with CDR3 region 
			if ( annotations[ assign.seqIdx ].cdr[2].seqIdx != -1 
				&& assign.seqEnd > annotations[ assign.seqIdx ].cdr[2].readStart + 3
				&& assign.seqStart < annotations[ assign.seqIdx ].cdr[2].readEnd - 3 )
			{
				struct  _assignRead nr ;
				nr.id = assembledReads[i].id ; 
				nr.read = strdup( assembledReads[i].read ) ;
				nr.overlap = assign ;
				if ( assign.strand == -1 )
				{
					seqSet.ReverseComplement( nr.read, assembledReads[i].read, strlen( nr.read ) ) ;
					nr.overlap.strand = 1 ;
				}
				cdr3Reads.push_back( nr ) ;
			}
			
			// Process the CDR3 region
			if ( annotations[assign.seqIdx].cdr[2].seqIdx != -1 
					&& assign.seqStart <= annotations[ assign.seqIdx ].cdr[2].readStart 
					&& assign.seqEnd >=  annotations[ assign.seqIdx ].cdr[2].readEnd )
			{
				std::vector<struct _CDR3info> &info = cdr3Infos[ assign.seqIdx ] ;
				int size = info.size() ;
				
				char *seq = assembledReads[i].read ;
				int offset = assign.readStart +
					annotations[ assign.seqIdx ].cdr[2].readStart - assign.seqStart ;
				if ( assign.strand == 1 )
				{
					memcpy( buffer, seq + offset, sizeof( char ) * cdr3Len ) ;
				}
				else if ( assign.strand == -1 )
				{
					int len = strlen( seq ) ;
					seqSet.ReverseComplement( buffer, seq + ( len - 1 - offset ) - cdr3Len + 1, cdr3Len ) ;
				}
				else
					continue ;
				buffer[cdr3Len] = '\0' ;

				for ( j = 0 ; j < size ; ++j )
					if ( !strcmp( info[j].seq, buffer ) )
					{
						++info[j].count ;
						break ;
					}
				if ( j >= size )
				{
					struct _CDR3info nc ;
					nc.seq = strdup( buffer ) ;
					nc.count = 1 ;

					info.push_back( nc ) ;
				}
				//printf( "%s\n%s %d\n", cdr3Reads[ cdr3Reads.size() - 1 ].read, buffer, i ) ;
			}
		}
		
		// Compute the abundance for each CDR3.
		// If there is no read span the whole CDR3 region, use consensus.
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
			}
		}

		// Distribute the reads.
		std::sort( cdr3Reads.begin(), cdr3Reads.end(), CompSortAssignedReadsByAssign  ) ;
		int cdr3ReadCnt = cdr3Reads.size() ;
		for ( i = 0 ; i < cdr3ReadCnt ; )
		{
			for ( j = i + 1 ; j < cdr3ReadCnt ; ++j )
			{
				if ( cdr3Reads[j].overlap.seqIdx != cdr3Reads[i].overlap.seqIdx )
					break ;
			}
			
			std::vector<struct _CDR3info> &info = cdr3Infos[ cdr3Reads[i].overlap.seqIdx ] ;
			int size = info.size() ;
			if ( size == 1 )
			{
				int cnt = 0 ;
				for ( k = i ; k < j ; ++k )
				{
					if ( k < j - 1 && !strcmp( cdr3Reads[k].id, cdr3Reads[k + 1].id ) )
						++k ;
					++cnt ;
				}
				info[0].count = cnt ;
				i = j ;
				continue ;
			}
			// Create the compatility relation
			std::vector< SimpleVector<int> > compat ; // read k is compatible to CDR3 l
			for ( k = i ; k < j ; ++k )
			{
				int l ;
				SimpleVector<int> nc ;
				
				if ( k < j - 1 && !strcmp( cdr3Reads[k].id, cdr3Reads[k + 1].id ) )
				{
					// Fragment overlap with the region.
					for ( l = 0 ; l < size ; ++l )
					{
						if ( IsCDR3Compatible( cdr3Reads[k], info[l], annotations[ cdr3Reads[i].overlap.seqIdx ].cdr[2] ) 
							&& IsCDR3Compatible( cdr3Reads[k + 1], info[l], 
									annotations[ cdr3Reads[i].overlap.seqIdx ].cdr[2] ) )
						{
							nc.PushBack( l ) ;
						}
					}
					++k ;
				}
				else
				{
					for ( l = 0 ; l < size ; ++l )
					{
						if ( IsCDR3Compatible( cdr3Reads[k], info[l], annotations[ cdr3Reads[i].overlap.seqIdx ].cdr[2] ) )
						{
							//printf( "%d=>%d\n", k, l ) ;
							nc.PushBack( l ) ;
						}
					}
				}

				compat.push_back( nc ) ;
			}
			// EM algorithm to estimate the count
			AbundanceEstimation( compat, info ) ;

			i = j ;
		}

		// Output different CDR3s from each main assembly.
		sprintf( buffer, "%s_cdr3.out", outputPrefix ) ;
		FILE *fpOutput = fopen( buffer, "w" ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( annotations[i].cdr[2].seqIdx == -1 )
				continue ;

			if ( !includePartial && annotations[i].cdr[2].similarity == 0 )
				continue ;

			std::vector<struct _CDR3info> &info = cdr3Infos[i] ;
			int size = info.size() ;
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
				fprintf( fpOutput, "%s\t%.2lf\t%.2lf\n", info[j].seq, annotations[i].cdr[2].similarity, info[j].count ) ;
			}
			
		}
		// Free up memory
		for ( i = 0 ; i < cdr3ReadCnt ; ++i )
		{
			free( cdr3Reads[i].read ) ;
		}
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = cdr3Infos[i].size() ;
			for ( j = 0 ; j < size; ++j )
				free( cdr3Infos[i][j].seq ) ;
		}

		for ( i = 0 ; i < assembledReadCnt ; ++i )
		{
			free( assembledReads[i].id ) ;
			free( assembledReads[i].read ) ;
		}

		fclose( fpOutput ) ;
		fclose( fpReads ) ; 
	}
	
	delete[] annotations ;
	PrintLog( "Finish annotation." ) ;
	return 0 ;
}
