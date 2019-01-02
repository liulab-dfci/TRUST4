#ifndef _MOURISL_KMERCODE_HEADER
#define _MOURISL_KMERCODE_HEADER

#include <stdio.h>
#include "defs.h"

extern char nucToNum[26] ;
extern char numToNuc[26] ;

class KmerCode
{
	private:
		int kmerLength ;
		int invalidPos ; // The position contains characters other than A,C,G,T in the code
		uint64_t code ;
		uint64_t mask ;
	public: 
		KmerCode() 
		{
		}
		KmerCode( int kl ) 
		{
			int i ;
			kmerLength = kl ;
			code = 0 ;
			invalidPos = -1 ;

			mask = 0 ;
			for ( i = 0 ; i < kmerLength ; ++i )
			{
				mask = mask << 2 ;
				mask = mask | 3 ;
			}
		}

		KmerCode( const KmerCode& in )
		{
			kmerLength = in.kmerLength ;
			invalidPos = in.invalidPos ;
			mask = in.mask ;
			code = in.code ;
		}

		void Restart() { code = 0ull ; invalidPos = -1 ; } 
		uint64_t GetCode() { return code ; }

		void SetCode( uint64_t in )
		{
			code = in ;
			invalidPos = -1 ;
		}

		uint64_t GetCanonicalKmerCode() 
		{
			int i ;
			uint64_t crCode = 0ull ; // complementary code
			for ( i = 0 ; i < kmerLength ; ++i )
			{
				//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
				//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;

				uint64_t tmp = ( code >> ( 2ull * i ) ) & 3ull ;
				crCode = ( crCode << 2ull ) | ( 3ull - tmp ) ;
			}
			return crCode < code ? crCode : code ;
		}

		uint64_t GetReverseComplementCode()
		{
			uint64_t crCode = 0ull ; // complementary code
			for ( int i = 0 ; i < kmerLength ; ++i )
			{
				//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
				//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;

				uint64_t tmp = ( code >> ( 2ull * i ) ) & 3ull ;
				crCode = ( crCode << 2ull ) | ( 3ull - tmp ) ;
			}
			return crCode ;
		}

		int GetKmerLength() { return kmerLength ; }

		bool IsValid() 
		{
			/*if ( invalidPos != -1 )
			{
				printf( "hi\n") ;
			}*/
			return ( invalidPos == -1 ) ;
		}

		void Append( char c ) 
		{
			if ( invalidPos != -1 )
				++invalidPos ;

			code = ( ( code << 2ull ) & mask ) | 
				( (uint64_t)( nucToNum[ c - 'A' ] & 3 ) ) ;

			//if ( nucToNum[c - 'A'] == -1 )
			if ( c == 'N' )
			{
				invalidPos = 0 ;	
			}
			if ( invalidPos >= kmerLength )
				invalidPos = -1 ;
		}

		void Prepend( char c )
		{
			ShiftRight( 1 ) ;
			if ( nucToNum[c-'A'] == -1 )
			{
				invalidPos = kmerLength - 1 ;
			}
			code = ( code | ( (uint64_t)( nucToNum[c - 'A'] & 3 ) << ( 2ull * ( kmerLength - 1 ) ) ) ) & mask ;
		}

		inline void ShiftRight( int k ) 
		{
			if ( invalidPos != -1 )
				invalidPos -= k ;

			code = ( code >> ( 2ull * k ) ) & ( mask >> ( 2ull * k ) ) ;	

			if ( invalidPos < 0 )
				invalidPos = -1 ;
		}

		bool IsEqual( const KmerCode &in )
		{
			if ( code == in.code )
				return true ;
			else
				return false ;
		}

		KmerCode& operator=( const KmerCode& in ) 
		{
			kmerLength = in.kmerLength ;
			invalidPos = in.invalidPos ;
			mask = in.mask ;
			code = in.code ;

			return *this ;
		}
} ;

#endif 
