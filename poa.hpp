#ifndef _MOURISL_POA_HEADER
#define _MOURISL_POA_HEADER

#include <string.h>

#include "SimpleVector.hpp"
#include "defs.h"

struct _POAnode
{
	//2 bits: character|
	short info ; 
	int edgeHeader ;
} ;

struct _POAedge
{
	unsigned char weight ;
	int nodeInd ;
	int next ;
} ;

struct _inEdge
{
	int nodeInd ;
	int next ;

	int poaEdge ;
} ;

extern char nucToNum[26] ;
extern char numToNuc[26] ;

class POA
{
private:
	SimpleVector<struct _POAnode> nodes ;
	SimpleVector<struct _POAedge> edges ;

	
	char *consensus ;
	int consensusLen ; // We will make sure that the nodes with id less than consensusLen is in the consensus.
	
	int AddNode( char c )
	{
		struct _POAnode n ;
		n.edgeHeader = -1 ;
		n.info = ( nucToNum[ c - 'A' ] ) & 3 ;
		
		return nodes.PushBack( n ) - 1 ;
	}
	
	int AddEdge( int u, int v, unsigned char weight = 1 )
	{
		int p = nodes[u].edgeHeader ;

		struct _POAedge e ;
		e.nodeInd = v ;
		e.weight = weight ;
		e.next = p ;
		nodes[u].edgeHeader = edges.Size() ;
		//printf( "addEdge: %d\n", nodes[u].edgeHeader ) ;
		return edges.PushBack( e ) - 1 ;
	}

	int ChangeWeight_ByNode( int u, int v, int weight = 1 )
	{
		int e = nodes[u].edgeHeader ;
		while ( e != -1 )	
		{
			if ( edges[e].nodeInd == v )
				break ;			
			e = edges[e].next ;
		}
		if ( e == -1 )
			return -1 ; 
		else
			edges[e].weight += weight ;
		return e ;
	}

	int ChangeWeight_ByLabel( int u, char c, int weight = 1 )
	{
		int e = nodes[u].edgeHeader ;
		while ( e != -1 )	
		{
			int v = edges[e].nodeInd ;
			if ( ( nodes[v].info & 3 ) == ( nucToNum[ c - 'A' ] & 3 ) )
				break ;
			e = edges[e].next ;
		}
		if ( e == -1 )
			return -1 ; 
		else
			edges[e].weight += weight ;
		return e ;
	}

public:
	POA() { consensus = NULL ; }
	
	POA( const POA &in )
	{
		if ( in.consensus != NULL )
		{
			consensus = strdup( in.consensus ) ;
			consensusLen = in.consensusLen ;

			nodes.Overwrite( in.nodes ) ;
			edges.Overwrite( in.edges ) ;
		}
		else
			consensus = NULL ; 
	}
	
	~POA() 
	{
		if ( consensus != NULL )
			free( consensus ) ;
	}

	int MemoryUsage()
	{
		return nodes.Capacity() * sizeof( struct _POAnode ) + edges.Capacity() * sizeof( struct _POAedge ) + 
			strlen( consensus ) + sizeof( consensusLen ) ;
	}
	
	void Release()
	{
		//if ( consensus != NULL )
		//	free( consensus ) ;
	}

	char *GetConsensus()
	{
		return consensus ;
	}

	int GetConsensusLen()
	{
		return consensusLen ;
	}

	void Initialize( char *seq, int len, bool onlySeq )
	{
		int i ;

		consensus = strdup( seq ) ;
		consensusLen = len ;

		if ( onlySeq ) // the poa does not need update
			return ;
		
		nodes.Reserve( len ) ;
		edges.Reserve( len ) ;
		
		nodes.SetMaxInc( len ) ;			
		edges.SetMaxInc( len ) ;
		//printf( "%d\n", len ) ;
		//nodes[0] is a dummy source node.
		for ( i = 0 ; i <= len ; ++i )
		{
			struct _POAnode n ;
			if ( i > 0 )
				n.info = nucToNum[ seq[i - 1] - 'A' ] ;
			else	
				n.info = 0 ;
			n.edgeHeader = i ;
			nodes.PushBack(n) ;
		}
		nodes[i - 1].edgeHeader = -1 ;

		for ( i = 0 ; i <= len - 1 ; ++i )
		{
			struct _POAedge e ;
			e.nodeInd = i + 1 ;
			e.weight = 0 ; // The initial value set to 0, so we allow the base read itself to be added again. 		
			e.next = -1 ;
			edges.PushBack( e ) ;
		}
	}

	// We have a alignment between consensus[start,...,end] with p,
	// update the POA in that aligned region.
	// We kind of merge the node with the same in-node.
	void UpdatePOA( int start, int end, char *p, int pLen, char *align )
	{
		int i ;
		int posc = 0, posp = 0 ;
		int prevTag = start ;
		int e ;
		const int offset = start + 1 ;
		
		for ( i = 0 ; align[i] != -1 ; ++i )
		{
			if ( align[i] == EDIT_MATCH )	
			{
				if ( ChangeWeight_ByNode( prevTag, posc + offset ) == -1 )
					AddEdge( prevTag, posc + offset ) ;

				prevTag = posc + offset ;
				++posc ;
				++posp ;
			}
			else //if ( align[i] == EDIT_MISMATCH )
			{
				//printf( "hi1: %d %d %c %c\n", prevTag, posp, consensus[ posc + start], p[posp] ) ;
				if ( ( e = ChangeWeight_ByLabel( prevTag, p[ posp ] ) ) == -1 )
				{
					int v = AddNode( p[posp] ) ;
					e = AddEdge( prevTag, v ) ;
				}
				//printf( "e=%d %d\n", e, edges.Size() ) ;
				prevTag = edges[e].nodeInd ;
				
				if ( align[i] == EDIT_MISMATCH )
				{
					++posc ;
					++posp ;
				}
				else if ( align[i] == EDIT_INSERT )
				{
					++posp ;
				}
				else if ( align[i] == EDIT_DELETE )
				{
					++posc ;
				}
				//printf( "hi2\n" ) ;
			}
		}
		// point back to the seed hit
		if ( ChangeWeight_ByNode( prevTag, posc + offset ) == -1 )
			AddEdge( prevTag, posc + offset ) ;
	}

	// Three quick operations.
	// If there is a kmer hit matches the consensus.
	void UpdatePOA_ConsecutiveMatch( int start, int matchLength ) 
	{
		//++start ;
		int i ;
		for ( i = start ; i <= start + matchLength - 1 ; ++i )
			ChangeWeight_ByNode( i, i + 1 ) ;
	}

	// There is a insertion of a sequence between consensus[start,...,start+1]
	void UpdatePOA_Insert( int start, char *p, int plen )
	{
		int prevTag = start ;
		int i, e ;
		for ( i = 0 ; i < plen ; ++i )
		{
			e = ChangeWeight_ByLabel( prevTag, p[i] ) ;
			if ( e == -1 )
			{
				int v = AddNode( p[i] ) ;
				e = AddEdge( prevTag, v ) ;
			}
			prevTag = edges[e].nodeInd ;
		}
	}

	// There is a deletion of consensus[start,...,end]
	void UpdatePOA_Delete( int start, int end )
	{
		// TODO: Do we need an edge to represent this node can be the end of the consensus?
		if ( end == consensusLen - 1 )
			return ;

		int e = ChangeWeight_ByNode( start, end + 2 ) ;
		if ( e == -1 )
			AddEdge( start, end + 2 ) ;
	}

	// Reorganize the poa, update the consensus and modify the kmer tree accordingly 
	void ReOrganize()
	{
		// Find the new consensus
		
		// Reorganize the poa

		// Update the index.
	} 

	void PrintPOA()
	{
		int i, e ;
		int nsize = nodes.Size() ;		
		for ( i = 0 ; i < nsize ; ++i )
		{
			printf( "%d %c:", i, numToNuc[ nodes[i].info & 3 ] ) ;
			e = nodes[i].edgeHeader ;
			while ( e != -1 )
			{
				printf( " (%d %u)", edges[e].nodeInd, edges[e].weight ) ;
				e = edges[e].next ;
			}
			printf( "\n" ) ;
		}
	}

} ;


#endif
