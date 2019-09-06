#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: FilterAnnotatedAssembly.pl annot.out > filter_annot.out\n" if ( @ARGV == 0 ) ;

open FP1, $ARGV[0] ;
while ( <FP1> )
{
	my $header ;
	my $seq ;
	chomp ;
	$header = $_ ;
	$seq = <FP1> ;
	chomp $seq ;

	next if ( $header =~ /null/ ) ;# Missing CDR3

	my @cols = split /\s/, $header ;
	if ( $header =~ /\* / ) # More the assembly missing a part, their CDR3 score must be good
	{
		if ( $cols[6] =~ /\):(.+?)=/ )
		{
			if ( $1 >= 100 ) 
			{
				print "$header\n$seq\n" ;
				#print "hello! ", $cols[6], " $1 ", $1>=80, "\n" ; 
			}
			next ; 
		}
		else
		{
			die "Wrong format $header\n" ;
		}
	}
	
	# Obtain the VJC coordinate
	my @vCoord ;
	my @jCoord ;
	my @cCoord ;
	my $cdr3 ;
	if ( $cols[3] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@vCoord = ($1, $2, $3, $4, $5) ;
	}
	else
	{
		die "Wrong format $header\n" ;
	}
	if ( $cols[4] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@jCoord = ($1, $2, $3, $4, $5) ;
	}
	else
	{
		die "Wrong format $header\n" ;
	}
	if ( $cols[5] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@cCoord = ($1, $2, $3, $4, $5 ) ;
	}
	else
	{
		die "Wrong format $header\n" ;
	}
	$cdr3 = ( split /=/, $cols[6] )[1] ;

	#if ( scalar( @vCoord ) < 5 || scalar( @jCoord ) < 5 || scalar( @cCoord ) < 5 ) 
	#{
	#	print( "ERROR ", @vCoord, " ", @jCoord, " ", @cCoord, "\n" ) ;
	#	exit() ;
	#}

	# The VJC coordinate should be correct
	next if ( $vCoord[2] > $jCoord[1] || $jCoord[2] > $cCoord[1] + 6 ) ;

	# Each should end within reasonable range in the annotation.
	next if ( $vCoord[4] < $vCoord[0] - length( $cdr3 ) || $jCoord[3] > length( $cdr3 ) 
		|| $jCoord[4] < $jCoord[0] - 20 || $cCoord[3] > 20 ) ;

	print "$header\n$seq\n" ;
}
close FP1 ;
