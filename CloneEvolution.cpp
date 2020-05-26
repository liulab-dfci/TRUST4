#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include <string>
#include <map>

#include "defs.h"

#define LINE_WIDTH 50000

char usage[] = "Usage: clone-evo xxx_cdr3.out [OPTIONS] > yyy_evo.out\n"
	"Options:\n"
	"\t-a FLOAT: only consider clonotypes more abundant than specified value (default: 2.0)\n"
	"\t-i STRING: the path to the file containing IGH isotype order, one per line (default: human order)\n"
	;

struct _cdr3
{
	char *seq ;
	int isotype ; 
	int clusterId ; // belongs to this cluster
	int subId ; // the id within this cluster.
	double abund ;
	double similarity ;

	bool operator<( const struct _cdr3 &b ) const
	{
		if ( clusterId != b.clusterId )
			return clusterId < b.clusterId ;
		else if ( isotype != b.isotype )
			return isotype < b.isotype ;
		//else if ( similarity != b.similarity )
		//	return similarity > b.similarity ;
		else
		{
			int tmp = strcmp(seq, b.seq) ;
			if ( tmp != 0 )
				return tmp < 0 ;
			else
				return abund > b.abund ;
		}
	}

	bool operator==(const struct _cdr3 &b) const
	{
		if (clusterId == b.clusterId && isotype == b.isotype && !strcmp(seq, b.seq))
			return true ;
		else
			return false ;
	}
} ;

struct _adj
{
	int v ;
	int next ;
} ;


int HammingDistance( char *a, char *b )
{
	int i, ret = 0 ;
	for ( i = 0 ; a[i] ; ++i )
		if ( a[i] != b[i] )
			++ret ;
	return ret ;
}

void InsertAdj( int u, int v, struct _adj *adj, int &adjUsed )
{
	adj[adjUsed].v = v ;
	adj[adjUsed].next = adj[u].next ;
	adj[u].next = adjUsed ;
	++adjUsed ;
}

bool IsCompatibleLevel( int level, int isotype )
{
	if ( isotype == -1 || level == -1 )
		return true ;
	return level == isotype ;
}


// Test whether t is on the path from "from" to the root
bool IsOnPath( int t, int from, struct _pair *minDist )
{
	int p = from ;
	while ( p != -1 )
	{
		if ( p == t )
			return true ;
		p = minDist[p].b ;
	}
	return false ;
}


// Prim algorithm for minimum spanning tree
// The algorithm is adjust to put the earlier isotypes aroudn the same levels.
int Prim( int **dist, int n,  struct _adj *adj, int &adjUsed, int offset, std::vector<struct _cdr3> &cdr3s )
{
	int i, j ;
	struct _pair *minDist = new struct _pair[n] ; // a-dist, b-the node used to connect
	bool *used = new bool[n] ;
	memset( used, false, sizeof(bool) * n ) ;
	
	int maxIsotype = cdr3s[n - 1].isotype ;
	int minIsotype = -1 ;
	int currLevelStart = 0 ;
	for ( i = 0 ; i < n ; ++i )
		if ( cdr3s[offset + i].isotype != -1 )
		{
			minIsotype = cdr3s[offset + i].isotype ;
			currLevelStart = i ;
			break ;
		}

	double max = -2 ;
	int maxTag = -1 ;
	int currLevel = minIsotype ;
	//fprintf( stderr, "%d %d %d %d %s\n", n, minIsotype, currLevelStart, cdr3s[offset].isotype, cdr3s[offset].seq ) ;
	for ( i = 0 ; i < n ; ++i )
	{
		if ( !IsCompatibleLevel( currLevel, cdr3s[i + offset].isotype ) )
			continue ;
		if ( cdr3s[i + offset].similarity > max 
			|| (cdr3s[i + offset].similarity == max 
				&& cdr3s[i + offset].abund > cdr3s[maxTag + offset].abund ) )
		{
			max = cdr3s[i + offset].similarity ;
			maxTag = i ;
		}
	}
	minDist[maxTag].a = 0 ;
	minDist[maxTag].b = -1 ;
	for ( i = 0 ; i < n ; ++i )
	{
		if ( i == maxTag )
			continue ;
		minDist[i].a = dist[maxTag][i] ;
		minDist[i].b = maxTag ;
	}
	used[maxTag] = true ;
	
	for ( i = 1 ; i < n ; ++i )
	{
		int min = LINE_WIDTH ;
		int minTag = -1 ;
		
		// Test whether current level has already been used up.
		bool goNextLevel = true ;

		for ( j = currLevelStart ; j < n ; ++j )
		{
			if ( cdr3s[j + offset].isotype != currLevel )
				break ;
			if ( !used[j] && cdr3s[j + offset].isotype == currLevel )
			{
				goNextLevel = false ;
				break ;
			}
		}
		if ( goNextLevel )
		{
			if ( j < n )
				currLevel = cdr3s[j + offset].isotype ;
			else
				currLevel = -1 ;
			currLevelStart = j ;
		}

		for ( j = 0 ; j < n ; ++j )
		{
			if ( cdr3s[j + offset].isotype > currLevel )
				break ;
			if ( !used[j] ) // The currLevel condition make sure that earlier isotypes are all used.
			{
				if ( minDist[j].a < min )
				{
					min = minDist[j].a ;
					minTag = j ;
				}
				else if ( minDist[j].a == min && cdr3s[j + offset].abund > cdr3s[minTag + offset].abund )
				{
					minTag = j ;
				}
			}
		}

		used[minTag] = true ;
		InsertAdj( minDist[minTag].b, minTag, adj, adjUsed ) ;
		InsertAdj( minTag, minDist[minTag].b , adj, adjUsed ) ;

		for ( j = 0 ; j < n ; ++j )
		{
			if ( !used[j] )
			{
				if ( dist[minTag][j] < minDist[j].a )
				{
					minDist[j].a = dist[minTag][j] ;
					minDist[j].b = minTag ;
				}
				else if ( dist[minTag][j] == minDist[j].a 
					&& cdr3s[minTag + offset].abund > cdr3s[ minDist[j].b + offset ].abund 
					//&& cdr3s[ minDist[j].b + offset].isotype >= cdr3s[minTag + offset].isotype )
					&& !IsOnPath( minDist[j].b, minTag, minDist ) )
				{
					minDist[j].b = minTag ;
				}
			}
		}
	}

	delete[] used ;
	delete[] minDist ;
	
	return maxTag ; 
}

// Return: whether it is a leaf.
int LongestRootLeafDistMST( int tag, struct _adj *adj, int n, bool *visited, int **dist )
{
	visited[tag] = true ;

	int p = adj[tag].next ;
	int ret = 0 ;
	while ( p != -1 )
	{
		if ( visited[ adj[p].v ] == false )
		{
			int tmp = dist[tag][ adj[p].v ] + LongestRootLeafDistMST( adj[p].v, adj, n, visited, dist ) ;
			if ( tmp > ret )
				ret = tmp ;
		}
		p = adj[p].next ;
	}
	return ret ;
}

//Newick tree format
void PrintMST( FILE *fp, int tag, int parent, bool *visited, int **dist, int offset, struct _adj *adj, 
	std::vector<struct _cdr3> &cdr3, std::vector<std::string> clusterIdToName )
{
	visited[tag] = true ;

	int p = adj[tag].next ;
	int childCnt = 0 ;
	while ( p != -1 )
	{
		if ( visited[adj[p].v] == false )
		{
			if ( childCnt == 0 )
				fprintf( fp, "(" ) ;
			else
				fprintf( fp, "," ) ;
			PrintMST( fp, adj[p].v, tag, visited, dist, offset, adj, cdr3, clusterIdToName ) ;
			++childCnt ;
		}
		p = adj[p].next ;
	}
	if ( childCnt > 0 )
		fprintf( fp, ")" ) ;
	if ( parent == tag ) // root
		fprintf( fp, "%s_%d_%.2lf;\n", clusterIdToName[ cdr3[tag + offset].clusterId ].c_str(), cdr3[tag + offset].subId,
			cdr3[tag + offset].abund ) ; 
	else
		fprintf( fp, "%s_%d_%.2lf:%d", clusterIdToName[ cdr3[tag + offset].clusterId ].c_str(), 
			cdr3[tag + offset].subId, cdr3[tag + offset].abund, dist[parent][tag] ) ;
}

// 0-not from ig. 1-from heavy chain. 2-from light chain.
int GetIgType(char *v, char *j, char *c)
{
	char *p ;
	if ( v[0] != '*' )
		p = v ;
	else if ( j[0] != '*' )
		p = j ;
	else if ( c[0] != '*' )
		p = c ;
	else
		return 0 ;

	if ( p[0] != 'I' || p[1] != 'G' )
		return 0 ;

	if ( p[2] == 'H' )
		return 1 ;
	else
		return 2 ;
	
	return 0 ;
}


char buffer[LINE_WIDTH] ;
char buffer2[LINE_WIDTH] ;
char seq[LINE_WIDTH] ;

std::map<std::string, int> clusterNameToId ;
std::vector<std::string> clusterIdToName ;

int main( int argc, char *argv[] )
{
	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}

	int i, j, k ;

	double minAbund = 2 ;

	std::map<std::string, int> isotypeRank ;

	for ( i = 2 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-a" ) )
		{
			minAbund = atof( argv[i + 1] ) ;
			i += 1 ;
		}
		else 
		{
			fprintf( stderr, "Unknown parameter %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	if ( isotypeRank.size() == 0 )
	{
		char humanIsotype[9][6] = {"IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2"} ;
		for ( i = 0 ; i < 9 ; ++i )
			isotypeRank[ std::string(humanIsotype[i]) ] = i ;
	}
	isotypeRank[ std::string("*") ] = -1 ;

	FILE *fp ;
	fp = fopen( argv[1], "r" ) ;

	std::vector<struct _cdr3> cdr3s ;
	while ( fgets( buffer, sizeof(buffer), fp ) != NULL )
	{
		char clusterIdBuffer[1024]  ;
		char vGene[1024], jGene[1024], cGene[1024] ;
		char similarityBuffer[1024] ;
		double score, abund, similarity ;
		sscanf( buffer, "%s %d %s %s %s %s %s %s %s %lf %lf %s", clusterIdBuffer, &k, 
			vGene, buffer2, jGene, cGene,
			buffer2, buffer2, seq,
			&score, &abund, similarityBuffer ) ;
		
		if ( strstr( similarityBuffer, "nan" ) )
			similarity = 0 ;
		else
			sscanf( similarityBuffer, "%lf", &similarity ) ;

		if ( score == 0 )
			continue ;

		// Only IG chains need evolution 
		int igType = GetIgType( vGene, jGene, cGene ) ;
		if ( igType == 0 )
			continue ;
		struct _cdr3 nc ;
		int clusterId ;
		std::string clusterName( clusterIdBuffer ) ;
		if ( clusterNameToId.find( clusterName ) != clusterNameToId.end() )
		{
			clusterId = clusterNameToId[clusterName] ;
		}
		else
		{
			clusterId = clusterIdToName.size() ;
			clusterIdToName.push_back( clusterName ) ;
			clusterNameToId[ clusterName ] = clusterId ;
		}
		nc.clusterId = clusterId ;
		nc.subId = k ;
		nc.abund = abund ;
		nc.similarity = similarity ;
		if ( igType == 2 )
			nc.isotype = 0 ;
		else
			nc.isotype = isotypeRank[std::string(cGene)] ;
		for ( i = 0 ; seq[i] ; ++i )
			if ( seq[i] == 'N' )
				break ;
		if (seq[i] == 'N' )
			continue ;

		if ( nc.abund < minAbund )
			continue ;
		
		nc.seq = strdup( seq ) ;
		cdr3s.push_back( nc ) ;
	}

	std::sort( cdr3s.begin(), cdr3s.end() ) ;

	// Remove repeated cdr3s
	std::vector<struct _cdr3> allCdr3s = cdr3s ;
	cdr3s.clear() ;
	
	int cdr3Cnt = allCdr3s.size() ;
	for ( i = 0 ; i < cdr3Cnt ; )
	{
		for ( j = i + 1 ; j < cdr3Cnt ; ++j )
			if ( !(allCdr3s[j] == allCdr3s[i]) )
				break ;
		//if ( i > 0 )
		//	printf( "%lf %d %d\n", allCdr3s[i].abund, strcmp(allCdr3s[i].seq, allCdr3s[i - 1].seq), allCdr3s[i].isotype) ;
		cdr3s.push_back( allCdr3s[i] ) ;	
		i = j ;
	}
	//printf( "%d %d\n", cdr3Cnt, cdr3s.size() ) ;

	// Process cluster by cluster. Assume the input is already sorted by cluster
	cdr3Cnt = cdr3s.size() ;
	for ( i = 0 ; i < cdr3Cnt ; )
	{
		for ( j = i + 1 ; j < cdr3Cnt ; ++j )
			if ( cdr3s[j].clusterId != cdr3s[i].clusterId )
				break ;
		int n = j - i ;
		if ( n <= 2 )
		{
			i = j ;
			continue ;
		}
		int l ;
		int **dist = new int*[n] ;
		for ( k = 0 ; k < n ; ++k )
			dist[k] = new int[n] ;
		struct _adj *adj = new struct _adj[3 * n] ;
		int adjUsed = 0 ;
		
		// Build the distance matrix
		for ( k = 0 ; k < n ; ++k )
		{
			dist[k][k] = 0 ;
			for ( l = k + 1 ; l < n ; ++l )
				dist[k][l] = dist[l][k] = HammingDistance( cdr3s[k + i].seq, cdr3s[l + i].seq ) ;
			adj[k].next = -1 ;
		}
		adjUsed = k ;
		
		// Build the tree
		int bestRoot = Prim( dist, n, adj, adjUsed, i, cdr3s ) ;

		// Find the root for the tree
		std::vector<int> leaves ;
		for ( k = 0 ; k < n ; ++k )
		{
			if ( adj[ adj[k].next ].next == -1 ) // The algorithm should make sure adj[k].next is not -1
				leaves.push_back( k ) ;
		}

		bool *visited = new bool[n] ;
		memset( visited, false, sizeof( bool ) * n ) ;
		int minMaxRootLeafDist = LongestRootLeafDistMST( bestRoot, adj, n, visited, dist ) ;
		//int bestRootIsotype = cdr3s[i + bestRoot].isotype ; This assignment is wrong when the root is -1.
		int bestRootIsotype = isotypeRank.size() ;
		for (k = 0 ; k < n ; ++k)
		{
			if (cdr3s[i + k].isotype > -1 && cdr3s[i + k].isotype < bestRootIsotype )	
				bestRootIsotype = cdr3s[i + k].isotype ;
		}
		// Adjust the best root.
		for ( k = 0 ; k < n ; ++k )
		{
			if ( cdr3s[i + k].similarity < cdr3s[i + bestRoot].similarity )
				continue ;
			if ( !IsCompatibleLevel( bestRootIsotype, cdr3s[i + k].isotype ) )
				continue ;
			if ( k == bestRoot )
				continue ;

			memset( visited, false, sizeof( bool ) * n ) ;
			int max = LongestRootLeafDistMST( k, adj, n, visited, dist ) ;
			if ( max < minMaxRootLeafDist )
			{
				bestRoot = k ;
				minMaxRootLeafDist = max ;
			}
			else if ( max == minMaxRootLeafDist )
			{
				if ( cdr3s[k + i].abund > cdr3s[bestRoot +i].abund )
					bestRoot = k ;
			}
		}

		// Output the tree
		memset( visited, false, sizeof( bool ) * n ) ;
		PrintMST( stdout, bestRoot, bestRoot, visited, dist, i, adj, cdr3s, clusterIdToName ) ;
	
		delete[] visited ;
		delete[] adj ;
		for ( k = 0 ; k < n ; ++k )
			delete[] dist[k] ;
		delete[] dist ;

		i = j ;
		continue ;
	}
	fclose( fp ) ;
	
	cdr3Cnt = allCdr3s.size() ;
	for ( i = 0 ; i < cdr3Cnt ; ++i )
		free( allCdr3s[i].seq ) ;
	return 0 ;
}
