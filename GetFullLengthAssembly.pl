#!/usr/bin/env perl

use strict ;
use warnings ;


die "usage: GetFullLengthAssembly.pl annot.fa > full_length_annot.fa\n" if ( @ARGV == 0 ) ;

open FP1, $ARGV[0] ;
while ( <FP1> )
{
	my $header ;
	my $seq ;
	chomp ;
	$header =$_ ;
	$seq = <FP1> ;
	chomp $seq ;

	#next if ( $header =~ /null/ || $header =~ /\* /) ;# Missing components
	# Obtain the VJC coordinate
	my @cols = split /\s/, $header ;

	my @vCoord ;
	my @jCoord ;
	my @cCoord ;
	my @cdr3Coord ;
	my $cdr3 ;
	next if ( $cols[3] eq "*" || $cols[5] eq "*" || $cols[6] eq "*") ;
	if ( $cols[3] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@vCoord = ($1, $2, $3, $4, $5) ;
	}
	else
	{
		die "Wrong format $header\n" ;
	}
	if ( $cols[5] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@jCoord = ($1, $2, $3, $4, $5) ;
	}
	else
	{
		die "Wrong format $header\n" ;
	}
	if ( $cols[6] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@cCoord = ($1, $2, $3, $4, $5) ;
	}
	else
	{
		die "Wrong format $header\n" ;
	}

	next if ( $cols[9] =~ /:0\.00/ ) ;

	if ( $cols[9] =~ /CDR3\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@cdr3Coord = ($1, $2) ;		
	}
	else
	{
		die "Wrong format $header\n" ;
	}
	$cdr3 = ( split /=/, $cols[9] )[1] ;
		
	# The VJC coordinate should be correct
	next if ( $vCoord[2] > $jCoord[1] + 3 || $jCoord[2] > $cCoord[1] + 6 ) ;

	next if ( $vCoord[3] >= 10 || $vCoord[2] < $cdr3Coord[0] ) ;
	next if ( $jCoord[1] > $cdr3Coord[1] || $jCoord[4] < $jCoord[0] - 3 ) ;	
	next if ( $cCoord[3] > 10 ) ;
	next if ( $seq =~ /N/ ) ;
	
	print "$header\n$seq\n" ;
}
close FP1 ;
