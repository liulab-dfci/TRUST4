// Methods related to alignment and string operations
#ifndef _LSONG_ALIGNALGO_HEADER
#define _LSONG_ALIGNALGO_HEADER

#include "defs.h"

#define EDIT_MATCH 0
#define EDIT_MISMATCH 1
#define EDIT_INSERT 2
#define EDIT_DELETE 3

#define SCORE_MATCH 5
#define SCORE_MISMATCH (-2)
#define SCORE_GAPOPEN (-4)
#define SCORE_GAPEXTEND (-1)
#define SCORE_INDEL (-4)

struct _posWeight
{
	int count[4] ;

	struct _posWeight &operator+=( const struct _posWeight &rhs )
	{
		count[0] += rhs.count[0] ;
		count[1] += rhs.count[1] ;
		count[2] += rhs.count[2] ;
		count[3] += rhs.count[3] ;
	}
} ;

class AlignAlgo
{
private:
	static bool IsBaseEqual( const struct _posWeight &w, const char c )
	{
		int sum = w.count[0] + w.count[1] + w.count[2] + w.count[3] ;
		if ( sum == 0 || c == 'N' || sum < 3 * w.count[ nucToNum[ c - 'A' ] ] )
			return true ;
		return false ;
	}
public:
	static int GlobalAlignment( char *t, int lent, char *p, int lenp, char *align )
	{
		int *m, *e, *f ;
		
		int i, j ;
		int bmax = ( lent + 1 ) ;

		m = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		e = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // insertion (to the text) gap
		f = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // deletion (to the text) gap
		
		m[0] = e[0] = f[0] = 0 ;
		for ( i = 1 ; i <= lenp ; ++i )
		{
			e[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPEXTEND ; 
	
			f[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
			m[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
		}

		for ( j = 1 ; j <= lent ; ++j  )
		{
			f[0 + j] = SCORE_GAPOPEN + j * SCORE_GAPEXTEND ;
			
			e[0 + j] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
			m[0 + j] = SCORE_GAPOPEN + j * SCORE_GAPOPEN ;
		}

		for ( i = 1 ; i <= lenp ; ++i )
		{
			for ( j = 1 ; j <= lent ; ++j )
			{
				int score ;
				// for e
				score = e[ ( i - 1 ) * bmax + j ] + SCORE_GAPEXTEND ;
				score = MAX( score, m[ ( i - 1) * bmax + j ] + SCORE_GAPOPEN + SCORE_GAPEXTEND ) ;
				e[ i * bmax + j ] = score ;

				// for f
				score = f[ i * bmax + j - 1 ] + SCORE_GAPEXTEND ;
				score = MAX( score, m[ i * bmax + j - 1 ] + SCORE_GAPOPEN + SCORE_GAPEXTEND ) ; 
				f[ i * bmax + j ] = score ;

				// for m 
				// Note that the index in the matrix is 1+the offset in the string.
				score = m[ ( i - 1 ) * bmax + j - 1] + ( ( t[j - 1] == p[i - 1] || t[j - 1] == 'N' 
										|| p[i - 1] == 'N' )? SCORE_MATCH : SCORE_MISMATCH ) ;	
				//printf( "%d %d: %d. %d %d %d\n", i, j, score, m[ (i - 1)*bmax + j - 1], e[ i * bmax + j], f[ i * bmax + j] ) ;
				score = MAX( score, e[ i * bmax + j] ) ;
				score = MAX( score, f[ i * bmax + j] ) ;
				m[i * bmax + j ] = score ;
			}
		}
		
		/*printf( "m:\n" ) ;
		for ( i = 0 ; i <= lenp ; ++i )
		{
			for ( j = 0 ; j <= lent ; ++j )
				printf( "%d ", m[i * bmax + j ] ) ;
			printf( "\n" ) ;
		}*/
		int ret = m[ lenp * bmax + lent ] ;
		
		// Trace back.
		int tagi = lenp, tagj = lent ;
		int mat = 0 ; // 0-m,1-e,2-f, which matrix the backtrace is in.
		int tag = 0 ;

		while ( tagi > 0 || tagj > 0 )
		{
			//printf( "%d %d %d\n", tagi, tagj, mat ) ;
			if ( mat == 0 )
			{
				int max = e[tagi * bmax + tagj] ;
				int a = EDIT_INSERT ;

				if ( f[tagi * bmax + tagj] >= m[tagi * bmax + tagj] )
					a = EDIT_DELETE ;
				if ( tagi > 0 && tagj > 0 
					&& ( m[ ( tagi - 1 ) * bmax + tagj - 1] + 
						( ( t[tagj - 1] == p[tagi - 1] || t[tagj - 1] == 'N' || p[tagi - 1] == 'N' ) ? 
								SCORE_MATCH : SCORE_MISMATCH )  == m[tagi * bmax + tagj] ) )
				{
					if ( t[tagj - 1] == p[tagi - 1] || t[tagj - 1] == 'N' || p[tagi - 1] == 'N' )
						a = EDIT_MATCH ;
					else
						a = EDIT_MISMATCH ;
				}

				if ( a == EDIT_MATCH || a == EDIT_MISMATCH )
				{
					align[tag] = a ;
					++tag ;

					--tagi ; --tagj ;
				}
				else if ( a == EDIT_INSERT )
					mat = 1 ;
				else if ( a == EDIT_DELETE )
					mat = 2 ;
			}
			else if ( mat == 1 ) // insertion to the text
			{
				int a = EDIT_INSERT ;
				align[tag] = a ;
				++tag ;

				if ( tagi > 0 )
				{	
					if ( m[ (tagi - 1) * bmax + tagj] + SCORE_GAPOPEN + SCORE_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagi ;
						mat = 0 ;
					}
					else //if ( e[( tagi - 1 ) * bmax + tagj] + EDIT_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagi ;
						mat = 1 ;
					}
				}
				else
				{
					mat = 2 ;
				}
			}
			else if ( mat == 2 ) // deletion to the text
			{
				int a = EDIT_DELETE ;
				align[tag] = a ;
				++tag ;

				if ( tagj > 0 )
				{
					if ( m[ tagi * bmax + tagj - 1] + SCORE_GAPOPEN + SCORE_GAPEXTEND == f[tagi * bmax + tagj] )
					{
						--tagj ;
						mat = 0 ;
					}
					else //if ( e[( tagi - 1 ) * bmax + tagj] + EDIT_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagj ;
						mat = 2 ;
					}
				}
				else
				{
					mat = 1 ;
				}
			}
		}
		align[tag] = -1 ;
		for ( i = 0, j = tag - 1 ; i < j ; ++i, --j )
		{	
			char tmp = align[i] ;
			align[i] = align[j] ;
			align[j] = tmp ;
		}
		delete[] m ;
		delete[] e ;
		delete[] f ;
		
		return ret ;
	}
	
	static int GlobalAlignment_PosWeight( struct _posWeight *tWeights, int lent, char *p, int lenp, char *align )
	{
		int *m, *e, *f ;
		
		int i, j ;
		int bmax = ( lent + 1 ) ;

		m = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		e = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // insertion (to the text) gap
		f = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // deletion (to the text) gap
		
		m[0] = e[0] = f[0] = 0 ;
		for ( i = 1 ; i <= lenp ; ++i )
		{
			e[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPEXTEND ; 
	
			f[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
			m[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
		}

		for ( j = 1 ; j <= lent ; ++j  )
		{
			f[0 + j] = SCORE_GAPOPEN + j * SCORE_GAPEXTEND ;
			
			e[0 + j] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
			m[0 + j] = SCORE_GAPOPEN + j * SCORE_GAPOPEN ;
		}

		for ( i = 1 ; i <= lenp ; ++i )
		{
			for ( j = 1 ; j <= lent ; ++j )
			{
				int score ;
				// for e
				score = e[ ( i - 1 ) * bmax + j ] + SCORE_GAPEXTEND ;
				score = MAX( score, m[ ( i - 1) * bmax + j ] + SCORE_GAPOPEN + SCORE_GAPEXTEND ) ;
				e[ i * bmax + j ] = score ;

				// for f
				score = f[ i * bmax + j - 1 ] + SCORE_GAPEXTEND ;
				score = MAX( score, m[ i * bmax + j - 1 ] + SCORE_GAPOPEN + SCORE_GAPEXTEND ) ; 
				f[ i * bmax + j ] = score ;

				// for m 
				// Note that the index in the matrix is 1+the offset in the string.
				score = m[ ( i - 1 ) * bmax + j - 1] + ( IsBaseEqual( tWeights[j - 1], p[i - 1] ) ? SCORE_MATCH : SCORE_MISMATCH ) ;	
				//printf( "%d %d: %d. %d %d %d\n", i, j, score, m[ (i - 1)*bmax + j - 1], e[ i * bmax + j], f[ i * bmax + j] ) ;
				score = MAX( score, e[ i * bmax + j] ) ;
				score = MAX( score, f[ i * bmax + j] ) ;
				m[i * bmax + j ] = score ;
			}
		}
		
		/*printf( "m:\n" ) ;
		for ( i = 0 ; i <= lenp ; ++i )
		{
			for ( j = 0 ; j <= lent ; ++j )
				printf( "%d ", m[i * bmax + j ] ) ;
			printf( "\n" ) ;
		}*/
		int ret = m[ lenp * bmax + lent ] ;
		
		// Trace back.
		int tagi = lenp, tagj = lent ;
		int mat = 0 ; // 0-m,1-e,2-f, which matrix the backtrace is in.
		int tag = 0 ;

		while ( tagi > 0 || tagj > 0 )
		{
			//printf( "%d %d %d\n", tagi, tagj, mat ) ;
			if ( mat == 0 )
			{
				int max = e[tagi * bmax + tagj] ;
				int a = EDIT_INSERT ;

				if ( f[tagi * bmax + tagj] >= m[tagi * bmax + tagj] )
					a = EDIT_DELETE ;
				if ( tagi > 0 && tagj > 0 
					&& ( m[ ( tagi - 1 ) * bmax + tagj - 1] + 
						( IsBaseEqual( tWeights[tagj - 1], p[tagi - 1] ) ? 
								SCORE_MATCH : SCORE_MISMATCH )  == m[tagi * bmax + tagj] ) )
				{
					if ( IsBaseEqual( tWeights[tagj - 1], p[tagi - 1] ) ) 
						a = EDIT_MATCH ;
					else
						a = EDIT_MISMATCH ;
				}

				if ( a == EDIT_MATCH || a == EDIT_MISMATCH )
				{
					align[tag] = a ;
					++tag ;

					--tagi ; --tagj ;
				}
				else if ( a == EDIT_INSERT )
					mat = 1 ;
				else if ( a == EDIT_DELETE )
					mat = 2 ;
			}
			else if ( mat == 1 ) // insertion to the text
			{
				int a = EDIT_INSERT ;
				align[tag] = a ;
				++tag ;

				if ( tagi > 0 )
				{	
					if ( m[ (tagi - 1) * bmax + tagj] + SCORE_GAPOPEN + SCORE_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagi ;
						mat = 0 ;
					}
					else //if ( e[( tagi - 1 ) * bmax + tagj] + EDIT_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagi ;
						mat = 1 ;
					}
				}
				else
				{
					mat = 2 ;
				}
			}
			else if ( mat == 2 ) // deletion to the text
			{
				int a = EDIT_DELETE ;
				align[tag] = a ;
				++tag ;

				if ( tagj > 0 )
				{
					if ( m[ tagi * bmax + tagj - 1] + SCORE_GAPOPEN + SCORE_GAPEXTEND == f[tagi * bmax + tagj] )
					{
						--tagj ;
						mat = 0 ;
					}
					else //if ( e[( tagi - 1 ) * bmax + tagj] + EDIT_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagj ;
						mat = 2 ;
					}
				}
				else
				{
					mat = 1 ;
				}
			}
		}
		align[tag] = -1 ;
		for ( i = 0, j = tag - 1 ; i < j ; ++i, --j )
		{	
			char tmp = align[i] ;
			align[i] = align[j] ;
			align[j] = tmp ;
		}
		delete[] m ;
		delete[] e ;
		delete[] f ;
		
		return ret ;
	}
	
	
	// Semi-global alignment where one end is matched, the other side is free.
	//   the alignment score should be at least threshold
	static int GlobalAlignment_OneEnd( char *t, int lent, char *p, int lenp, int threshold, char *align )
	{
		int *m, *e, *f ;
		if ( lent == 0 || lenp == 0 )
		{
			align[0] = -1 ;
			return 0 ;
		}

		int i, j ;
		int bmax = ( lent + 1 ) ;

		m = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		e = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // insertion (to the text) gap
		f = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // deletion (to the text) gap
		
		m[0] = e[0] = f[0] = 0 ;
		for ( i = 1 ; i <= lenp ; ++i )
		{
			e[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPEXTEND ; 
	
			f[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
			m[i * bmax + 0] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
		}

		for ( j = 1 ; j <= lent ; ++j  )
		{
			f[0 + j] = SCORE_GAPOPEN + j * SCORE_GAPEXTEND ;
			
			e[0 + j] = SCORE_GAPOPEN + i * SCORE_GAPOPEN ; 
			m[0 + j] = SCORE_GAPOPEN + j * SCORE_GAPOPEN ;
		}

		for ( i = 1 ; i <= lenp ; ++i )
		{
			for ( j = 1 ; j <= lent ; ++j )
			{
				int score ;
				// for e
				score = e[ ( i - 1 ) * bmax + j ] + SCORE_GAPEXTEND ;
				score = MAX( score, m[ ( i - 1) * bmax + j ] + SCORE_GAPOPEN + SCORE_GAPEXTEND ) ;
				e[ i * bmax + j ] = score ;

				// for f
				score = f[ i * bmax + j - 1 ] + SCORE_GAPEXTEND ;
				score = MAX( score, m[ i * bmax + j - 1 ] + SCORE_GAPOPEN + SCORE_GAPEXTEND ) ; 
				f[ i * bmax + j ] = score ;

				// for m 
				// Note that the index in the matrix is 1+the offset in the string.
				score = m[ ( i - 1 ) * bmax + j - 1] + ( ( t[j - 1] == p[i - 1] || t[j - 1] == 'N' 
										|| p[i - 1] == 'N' )? SCORE_MATCH : SCORE_MISMATCH ) ;	
				//printf( "%d %d: %d. %d %d %d\n", i, j, score, m[ (i - 1)*bmax + j - 1], e[ i * bmax + j], f[ i * bmax + j] ) ;
				score = MAX( score, e[ i * bmax + j] ) ;
				score = MAX( score, f[ i * bmax + j] ) ;
				m[i * bmax + j ] = score ;
			}
		}
		
		/*printf( "m:\n" ) ;
		for ( i = 0 ; i <= lenp ; ++i )
		{
			for ( j = 0 ; j <= lent ; ++j )
				printf( "%d ", m[i * bmax + j ] ) ;
			printf( "\n" ) ;
		}*/

		// Locate the ending point.
		int max = threshold ;
		int tagi = 0, tagj = 0 ;

		for ( i = 0 ; i <= lenp ; ++i )
			for ( j = 0 ; j <= lent ; ++j )
			{
				int scoreThreshold = ( i + j ) * 0.5 * 0.75 * SCORE_MATCH ;
				if ( m[i * bmax + j] < scoreThreshold )
					continue ;
				if ( m[i * bmax +j] > max )
				{
					max = m[i * bmax + j] ;
					tagi = i ;
					tagj = j ;
				}
			}
		int ret = m[ tagi * bmax + tagj ] ;
		
		// Trace back.
		int mat = 0 ; // 0-m,1-e,2-f, which matrix the backtrace is in.
		int tag = 0 ;

		while ( tagi > 0 || tagj > 0 )
		{
			//printf( "%d %d %d\n", tagi, tagj, mat ) ;
			if ( mat == 0 )
			{
				int max = e[tagi * bmax + tagj] ;
				int a = EDIT_INSERT ;

				if ( f[tagi * bmax + tagj] >= m[tagi * bmax + tagj] )
					a = EDIT_DELETE ;
				if ( tagi > 0 && tagj > 0 
					&& ( m[ ( tagi - 1 ) * bmax + tagj - 1] + 
						( ( t[tagj - 1] == p[tagi - 1] || t[tagj - 1] == 'N' || p[tagi - 1] == 'N' ) ? 
								SCORE_MATCH : SCORE_MISMATCH )  == m[tagi * bmax + tagj] ) )
				{
					if ( t[tagj - 1] == p[tagi - 1] || t[tagj - 1] == 'N' || p[tagi - 1] == 'N' )
						a = EDIT_MATCH ;
					else
						a = EDIT_MISMATCH ;
				}

				if ( a == EDIT_MATCH || a == EDIT_MISMATCH )
				{
					align[tag] = a ;
					++tag ;

					--tagi ; --tagj ;
				}
				else if ( a == EDIT_INSERT )
					mat = 1 ;
				else if ( a == EDIT_DELETE )
					mat = 2 ;
			}
			else if ( mat == 1 ) // insertion to the text
			{
				int a = EDIT_INSERT ;
				align[tag] = a ;
				++tag ;

				if ( tagi > 0 )
				{	
					if ( m[ (tagi - 1) * bmax + tagj] + SCORE_GAPOPEN + SCORE_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagi ;
						mat = 0 ;
					}
					else //if ( e[( tagi - 1 ) * bmax + tagj] + EDIT_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagi ;
						mat = 1 ;
					}
				}
				else
				{
					mat = 2 ;
				}
			}
			else if ( mat == 2 ) // deletion to the text
			{
				int a = EDIT_DELETE ;
				align[tag] = a ;
				++tag ;

				if ( tagj > 0 )
				{
					if ( m[ tagi * bmax + tagj - 1] + SCORE_GAPOPEN + SCORE_GAPEXTEND == f[tagi * bmax + tagj] )
					{
						--tagj ;
						mat = 0 ;
					}
					else //if ( e[( tagi - 1 ) * bmax + tagj] + EDIT_GAPEXTEND == e[tagi * bmax + tagj] )
					{
						--tagj ;
						mat = 2 ;
					}
				}
				else
				{
					mat = 1 ;
				}
			}
		}
		align[tag] = -1 ;
		for ( i = 0, j = tag - 1 ; i < j ; ++i, --j )
		{	
			char tmp = align[i] ;
			align[i] = align[j] ;
			align[j] = tmp ;
		}
		delete[] m ;
		delete[] e ;
		delete[] f ;
		
		return ret ;
	}
	
	
	
	static double GlobalAlignment_classic( char *t, int lent, char *p, int lenp, char *align ) 
	{
		int i, j ;
		int *score = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		int bmax = ( lent + 1 ) ;
		for ( i = 0 ; i <= lenp ; ++i )
			score[i * bmax + 0] = i ;
		for ( j = 0 ; j <= lent ; ++j )
			score[0 * bmax + j] = j ;

		for ( i = 1 ; i <= lenp ; ++i )
		{	
			for ( j = 1 ; j <= lent ; ++j )
			{
				int min = score[ ( i - 1 ) * bmax + j] + 1 ;
				if ( score[i * bmax + j - 1] + 1 < min )
					min = score[i * bmax + j - 1] + 1 ;
				int diag = score[ ( i - 1 ) * bmax + j - 1] + ( t[j - 1] == p[i - 1] ? 0 : 1 ) ;
				if ( diag < min )
					min = diag ;
				score[i * bmax + j] = min ;
			}
		}

		int tagi = lenp, tagj = lent ;
		int tag = 0 ;
		while ( tagi > 0 || tagj > 0 )
		{
			int min = score[tagi * bmax + tagj] ;
			int a = 0 ;
			if ( tagj > 0 && score[tagi * bmax + tagj - 1] + 1 == min )
				a = EDIT_DELETE ;	
			if ( tagi > 0 && score[ ( tagi - 1 ) * bmax + tagj] + 1 == min )
				a = EDIT_INSERT ; 
			if ( tagj > 0 && tagi > 0 )
			{
				int diff = ( t[tagj - 1] == p[tagi - 1] ? 0 : 1 ) ;
				if ( score[ ( tagi - 1) * bmax + tagj - 1 ] + diff == min )
				{
					if ( diff == 0 )
						a = EDIT_MATCH ;
					else
						a = EDIT_MISMATCH ;
				}
			}

			align[ tag ] = a ;
			++tag ;
			if ( a == EDIT_DELETE )
				--tagj ;
			else if ( a == EDIT_INSERT )
				--tagi ;
			else
			{
				--tagi ;
				--tagj ;
			}
		}
		align[tag] = -1 ;
		for ( i = 0, j = tag - 1 ; i < j ; ++i, --j )
		{	
			char tmp = align[i] ;
			align[i] = align[j] ;
			align[j] = tmp ;
		}
		//printf( "%d\n", score[lenp * bmax + lent] ) ;
		return score[lenp * bmax + lent] ;
	}

	static void VisualizeAlignment( char *t, int lent, char *p, int lenp, char *align )
	{
		int i, k, j;
		int tagt, tagp, taga ;
		taga = 0 ;
		int width = 100 ;
		k = j = 0 ;
		/*for ( i = 0 ; align[i] != -1 ; ++i )
		  printf( "%d", align[i] ) ;
		  printf( "\n" ) ;*/
		while ( align[taga] != -1 )
		{
			for ( i = taga ; i < taga + width && align[i] != -1 ; ++i  )
			{
				if ( align[i] == EDIT_INSERT )
					printf( "-" ) ;
				else
				{
					printf( "%c", t[k] ) ;
					++k ;
				}
			}
			printf( "\n" ) ;
			for ( i = taga ; i < taga + width && align[i] != -1 ; ++i )
			{
				if ( align[i] == EDIT_MATCH )
					printf( "|" ) ;
				else
					printf( " " ) ;
			}
			printf( "\n" ) ;
			for ( i = taga ; i < taga + width && align[i] != -1 ; ++i  )
			{
				if ( align[i] == EDIT_DELETE )
					printf( "-" ) ;
				else
				{
					printf( "%c", p[j] ) ;
					++j ;
				}
			}
			printf( "\n\n" ) ;
			taga = i ;
		}
	}
} ;

#endif