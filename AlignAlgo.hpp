// Methods related to alignment and string operations
#ifndef _LSONG_ALIGNALGO_HEADER
#define _LSONG_ALIGNALGO_HEADER

#include "defs.h"

#define EDIT_MATCH 0
#define EDIT_MISMATCH 1
#define EDIT_INSERT 2
#define EDIT_DELETE 3

#define SCORE_MATCH 2
#define SCORE_MISMATCH (-2)
#define SCORE_GAPOPEN (-4)
#define SCORE_GAPEXTEND (-1)
#define SCORE_INDEL (-4)

#define SCORE_MATCH_LOCAL 1
#define SCORE_MISMATCH_LOCAL (-2)

struct _posWeight
{
	int count[4] ;

	struct _posWeight &operator+=( const struct _posWeight &rhs )
	{
		count[0] += rhs.count[0] ;
		count[1] += rhs.count[1] ;
		count[2] += rhs.count[2] ;
		count[3] += rhs.count[3] ;
	
		return *this ;
	}

	int Sum() const
	{
		return count[0] + count[1] + count[2] + count[3] ;
	}
} ;

class AlignAlgo
{
private:
	static bool IsBaseEqual( const struct _posWeight &w, const char c )
	{
		int sum = w.Sum() ;
		if ( sum == 0 || c == 'N' || sum < 3 * w.count[ nucToNum[ c - 'A' ] ] )
			return true ;
		return false ;
	}
public:
	static double GlobalAlignment_PosWeight( struct _posWeight *tWeights, int lent, char *p, int lenp, char *align ) 
	{
		if ( lent == 0 || lenp == 0 )
		{
			align[0] = -1 ;
			return 0 ;
		}
		else if ( lent == 1 && lenp == 1 )
		{
			if ( IsBaseEqual( tWeights[0], p[0] ) )
			{
				align[0] = EDIT_MATCH ;
				align[1] = -1 ;
				return SCORE_MATCH ;
			}
			else
			{
				align[0] = EDIT_MISMATCH ;
				align[1] = -1 ;
				return SCORE_MISMATCH ;
			}
		}

		int i, j ;
		if ( lent == lenp )
		{
			// Check whether no-indel alignment could work.
			// For them to match but with indel, we need at least one indel in t and one delin in p
			int score = 0 ;
			for ( i = 0 ; i < lent ; ++i )
			{
				if ( IsBaseEqual( tWeights[i], p[i] ) )
				{
					align[i] = EDIT_MATCH ;
					score += SCORE_MATCH ;
				}
				else
				{
					align[i] = EDIT_MISMATCH ;
					score += SCORE_MISMATCH ;
				}
			}
			align[i] = -1 ;

			if ( score >= lent * SCORE_MATCH + 2 * SCORE_INDEL )
				return score ;
		}


		int leftBand = 5 ;
		int rightBand = 5 ;
		if ( lent > lenp )
			rightBand += lent - lenp ;
		else if ( lent < lenp ) // more rows than column.
			leftBand += lenp - lent ;
		
		int negInf = ( lent + 1 ) * ( lenp + 1 ) * SCORE_INDEL ;
		int bmax = ( lent + 1 ) ;

		int *m = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		
		m[0] = 0 ;
		for ( i = 1 ; i <= lenp ; ++i )
		{
			m[i * bmax + 0] = SCORE_INDEL + i * SCORE_INDEL; 
		}

		for ( j = 1 ; j <= lent ; ++j  )
		{
			m[0 + j] = SCORE_INDEL + j * SCORE_INDEL ;
		}

		for ( i = 1 ; i <= lenp ; ++i )
		{
			int start = ( i - leftBand < 1 ) ? 1 : ( i - leftBand ) ;
			int end = ( i + rightBand > lent ) ? lent : ( i + rightBand ) ;

			if ( start > 1 )
			{
				j = start - 1 ;
				m[i * bmax + j] = negInf ;
			}
			if ( end < lent )
			{
				j = end + 1 ;
				m[i * bmax + j] = negInf ;
			}
			
			for ( j = start ; j <= end ; ++j )
			//for ( j = 1 ; j <= lent ; ++j )
			{
				int score ;
				score = m[ ( i - 1 ) * bmax + j - 1] + ( IsBaseEqual( tWeights[j - 1], p[i - 1] )? SCORE_MATCH : SCORE_MISMATCH ) ;	
				//printf( "%d %d: %d. %d %d %d\n", i, j, score, m[ (i - 1)*bmax + j - 1], e[ i * bmax + j], f[ i * bmax + j] ) ;
				score = MAX( score, m[ i * bmax + ( j - 1 ) ] + SCORE_INDEL ) ;
				score = MAX( score, m[ ( i - 1 ) * bmax + j ] + SCORE_INDEL ) ;
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
		int tag = 0 ;
		
		while ( tagi > 0 || tagj > 0 )
		{
			int max = m[tagi * bmax + tagj] ;
			int a = 0 ;
			if ( tagj > 0 && m[tagi * bmax + tagj - 1] + SCORE_INDEL == max )
				a = EDIT_DELETE ;	
			if ( tagi > 0 && m[ ( tagi - 1 ) * bmax + tagj] + SCORE_INDEL == max )
				a = EDIT_INSERT ; 
			if ( tagj > 0 && tagi > 0 )
			{
				int diff = ( IsBaseEqual( tWeights[ tagj - 1], p[tagi - 1] ) ? SCORE_MATCH : SCORE_MISMATCH ) ;
				if ( m[ ( tagi - 1) * bmax + tagj - 1 ] + diff == max )
				{
					if ( diff == SCORE_MATCH )
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
		delete[] m ;
		
		return ret ;
	}
		
	static int GlobalAlignment( char *t, int lent, char *p, int lenp, char *align )
	{
		if ( lent == 0 || lenp == 0 )
		{
			align[0] = -1 ;
			return 0 ;
		}
		else if ( lent == 1 && lenp == 1 )
		{
			if ( t[0] == p[0] || t[0] == 'N' || p[0] == 'N' )
			{
				align[0] = EDIT_MATCH ;
				align[1] = -1 ;
				return SCORE_MATCH ;
			}
			else
			{
				align[0] = EDIT_MISMATCH ;
				align[1] = -1 ;
				return SCORE_MISMATCH ;
			}
		}

		int *m, *e, *f ;

		int leftBand = 5 ;
		int rightBand = 5 ;
		if ( lent > lenp )
			rightBand += lent - lenp ;
		else if ( lent < lenp ) // more rows than column.
			leftBand += lenp - lent ;

		int i, j ;
		int negInf = ( lent + 1 ) * ( lenp + 1 ) * SCORE_GAPOPEN ;
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
			int start = ( i - leftBand < 1 ) ? 1 : ( i - leftBand ) ;
			int end = ( i + rightBand > lent ) ? lent : ( i + rightBand ) ;

			if ( start > 1 )
			{
				j = start - 1 ;
				e[i * bmax + j] = f[i * bmax + j] = m[i * bmax + j] = negInf ;
			}
			if ( end < lent )
			{
				j = end + 1 ;
				e[i * bmax + j] = f[i * bmax + j] = m[i * bmax + j] = negInf ;
			}
			
			for ( j = start ; j <= end ; ++j )
			//for ( j = 1 ; j <= lent ; ++j )
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
	
	static int GlobalAlignment_PosWeight_Affine( struct _posWeight *tWeights, int lent, char *p, int lenp, char *align )
	{
		if ( lent == 0 || lenp == 0 )
		{
			align[0] = -1 ;
			return 0 ;
		}
		else if ( lent == 1 && lenp == 1 )
		{
			if ( IsBaseEqual( tWeights[0], p[0] ) )
			{
				align[0] = EDIT_MATCH ;
				align[1] = -1 ;
				return SCORE_MATCH ;
			}
			else
			{
				align[0] = EDIT_MISMATCH ;
				align[1] = -1 ;
				return SCORE_MISMATCH ;
			}
		}
		//printf( "%d %d\n", lent, lenp ) ;
		int *m, *e, *f ;
		
		int i, j ;
		int bmax = ( lent + 1 ) ;
		int negInf = ( lent + 1 ) * ( lenp + 1 ) * SCORE_GAPOPEN ;
		
		m = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		e = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // insertion (to the text) gap
		f = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ; // deletion (to the text) gap
		int band = 5 ; // asumme lenp==lent.

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
			int start = ( i - band < 1 ) ? 1 : ( i - band ) ;
			int end = ( i + band > lent ) ? lent : ( i + band ) ;

			if ( start > 1 )
			{
				j = start - 1 ;
				e[i * bmax + j] = f[i * bmax + j] = m[i * bmax + j] = negInf ;
			}
			if ( end < lent )
			{
				j = end + 1 ;
				e[i * bmax + j] = f[i * bmax + j] = m[i * bmax + j] = negInf ;
			}
			for ( j = start ; j <= end ; ++j )
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
				int scoreThreshold = ( i + j ) * 0.5 *(  0.8 * SCORE_MATCH + 0.2 * SCORE_MISMATCH ) ;
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

				if ( f[tagi * bmax + tagj] >= max )
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
	
	static int LocalAlignment( char *t, int lent, char *p, int lenp, int &tstart, int &pstart, char *align )
	{
		int i, j ;
		int *m = new int[ ( lenp + 1 ) * ( lent + 1 ) ] ;
		int bmax = ( lent + 1 ) ;
		for ( i = 0 ; i <= lenp ; ++i )
			m[i * bmax + 0] = 0 ;
		for ( j = 0 ; j <= lent ; ++j )
			m[0 * bmax + j] = 0 ;
		
		tstart = 0 ;
		pstart = 0 ;
		for ( i = 1 ; i <= lenp ; ++i )
		{	
			for ( j = 1 ; j <= lent ; ++j )
			{
				int score ;
				score = m[ ( i - 1 ) * bmax + j - 1] + ( t[j - 1] == p[i - 1] ? SCORE_MATCH_LOCAL : SCORE_MISMATCH_LOCAL ) ;	
				score = MAX( score, m[ i * bmax + ( j - 1 ) ] + SCORE_INDEL ) ;
				score = MAX( score, m[ ( i - 1 ) * bmax + j ] + SCORE_INDEL ) ;
				if ( score < 0 )
					score = 0 ;
				m[i * bmax + j ] = score ;
			}
		}

		int tagi = lenp, tagj = lent ;
		int maxScore = 0 ;
		for ( i = 0 ; i <= lenp ; ++i )
			for ( j = 0 ; j <= lent ; ++j )
				if ( m[i * bmax + j] >= maxScore )
				{
					maxScore = m[i * bmax + j] ;
					tagi = i ;
					tagj = j ;
				}
		if ( maxScore == 0 )
		{
			delete[] m ;
			return -1 ;
		}
		int tag = 0 ;
		while ( tagi > 0 || tagj > 0 )
		{
			int max = m[tagi * bmax + tagj] ;
			int a = 0 ;
			if ( max == 0 )
			{
				tstart = tagj ;
				pstart = tagi ;
				break ;
			}
			if ( tagj > 0 && m[tagi * bmax + tagj - 1] + SCORE_INDEL == max )
				a = EDIT_DELETE ;	
			if ( tagi > 0 && m[ ( tagi - 1 ) * bmax + tagj] + SCORE_INDEL == max )
				a = EDIT_INSERT ; 
			if ( tagj > 0 && tagi > 0 )
			{
				int diff = ( t[ tagj - 1] == p[tagi - 1] ? SCORE_MATCH_LOCAL : SCORE_MISMATCH_LOCAL ) ;
				if ( m[ ( tagi - 1) * bmax + tagj - 1 ] + diff == max )
				{
					if ( diff == SCORE_MATCH_LOCAL )
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
		delete[] m ;
		//printf( "%d\n", score[lenp * bmax + lent] ) ;
		return maxScore ;
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

	static int IsMateOverlap( char *fr, int flen, char *sr, int slen, int minOverlap, int &offset, int &bestMatchCnt,
		bool checkTandem = true )
	{
		int i, j, k ;		
		bestMatchCnt = -1 ;
		int offsetCnt = 0 ;
		int overlapSize = -1 ;
		for ( j = 0 ; j < flen - minOverlap ; ++j ) // The overlap start position in first read
		{
			// Whether the overlap works.
			int matchCnt = 0 ;
			bool flag = true ;
			
			double similarityThreshold = 0.95 ;
			if ( flen - j >= 100 )
				similarityThreshold = 0.85 ;
			else if ( flen - j >= 50 )
				similarityThreshold = 0.85 + ( flen - j - 50 ) / 50.0 * 0.1 ;
			
			for ( k = 0 ; j + k < flen && k < slen ; ++k )
			{
				if ( fr[j + k] == sr[k] )
					++matchCnt ;
				if ( matchCnt + ( flen - ( j + k ) - 1 ) < int( ( flen - j ) * similarityThreshold ) )
				{
					flag = false ;
					break ;
				}
			}

			if ( flag ) 
			{
				offset = j ;
				++offsetCnt ;
				overlapSize = k ;
				bestMatchCnt = matchCnt ;
			}
		}

		if ( offsetCnt != 1 )
			return -1 ;
		
		// If the overlap size is near minOverlap, the overlap could still be ambiguous. 
		// i- the repeat size
		if ( checkTandem && overlapSize <= minOverlap * 2 )
		{
			for ( i = 1 ; i <= overlapSize / 2 ; ++i )
			{
				bool tandem = true ;
				for ( j = i ; j + i - 1 < overlapSize ; j += i )
				{
					for ( k = j ; k <= j + i - 1 ; ++k )
					{
						if ( sr[k -j] != sr[k] )
							break ;
					}
					if ( k <= j + i - 1 )
					{
						tandem = false ;
						break ;
					}
				}

				if ( tandem )
					return -1 ;
			}
		}
		
		return overlapSize ;
	}
	
	

	// Find the partial suffix of a that match with prefix of b of length at least minLen.
	// Return: the suffix of a. otherwise -1.
	static int LocatePartialSufPrefExactMatch( char *a, int lenA, char *b, int lenB, int minLen, int &matchLen )
	{
		int i, j, k ;
		int max = 0 ;
		int maxTag = 0 ;
		int secMax = 0 ;
		for ( k = 0 ; k + minLen - 1 < lenA ; ++k )
		{
			for ( i = k, j = 0 ; i < lenA && j < lenB ; ++i, ++j )
			{
				if ( a[i] != b[j] )
					break ;
			}
			/*if ( j + 1 >= minLen )
			{
				matchLen = j + 1 ;
				return k ;
			}*/

			if ( j + 1 > max )
			{
				secMax = max ;
				max = j + 1 ;
				maxTag = k ;
			}
			else if ( j + 1 >= secMax ) // This secMax==max when there are two max value.
				secMax = j + 1 ;
		}

		if ( max >= minLen && max > secMax + 1 )
		{
			matchLen = max ;
			return maxTag ;
		}
		matchLen = 0 ;
		return -1 ;
	}

	// Find the partial suffix of "a" that match with suffix of b of length at least minLen.
	// Return: the partial suffix of "a", otherwise -1.
	static int LocatePartialSufSufExactMatch( char *a, int lenA, char *b, int lenB, int minLen, int &matchLen )
	{
		int i, j, k ;
		// k indicating the end of the match.
		int max = 0 ;
		int maxTag = 0 ;
		int secMax = 0 ;
		for ( k = lenA - 1 ; k >= minLen ; --k )
		{
			for ( i = k, j = lenB - 1 ; i >= 0 && j >= 0 ; --i, --j )
				if ( a[i] != b[j] )
					break ;
			/*if ( k - i >= minLen )
			{
				matchLen = k - i ;
				return i + 1 ; 
			}*/
			if ( k - i > max )
			{
				secMax = max ;
				max = k - i ;
				maxTag = i + 1 ;
			}
			else if ( k - i >= secMax ) // This secMax==max when there are two max value.
				secMax = k - i ;

		}
		if ( max >= minLen && max > secMax + 1 )
		{
			matchLen = max ;
			return maxTag ;
		}

		matchLen = 0 ;
		return -1 ;
	}
} ;

#endif
