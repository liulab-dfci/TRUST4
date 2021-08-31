#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl xxx_cdr3.out xxx_annot.fa > xxx_with_seq_cdr3.out\n" if (@ARGV == 0) ;

my %cdr3Range ;
my %allSeq ;
open FP, $ARGV[1] ;
while (<FP>)
{
	chomp ;
	my $header = substr($_, 1) ;
	my $seq = <FP> ;
	chomp $seq ;
	my @cols = split /\s/, $header ;
	my $seqId = $cols[0] ;
	
	if ( $cols[9] =~ /CDR3\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@{$cdr3Range{$seqId}} = ($1, $2) ;
	}
	else 
	{
		die "Unknown format $header\n" ;
	}
	$allSeq{$seqId} = $seq ;
}
close FP ;

open FP, $ARGV[0] ;
while (<FP>)
{
	chomp ;
	my @cols = split ;
	my $cdr3 = $cols[8] ;
	my $seqId = $cols[0] ;
	my @cdr3Range = @{$cdr3Range{$seqId}} ;
	if ($cdr3Range[0] == 0 && $cdr3Range[1] == 0)
	{
		push @cols, $allSeq{$seqId} ;
	}
	else
	{
		#print $allSeq{$seqId}, " ", $cdr3Range[0], " ", $cdr3Range[1], "\n" ;
		my $tmp = $allSeq{$seqId} ;
		substr($tmp, $cdr3Range[0], $cdr3Range[1] - $cdr3Range[0] + 1, $cdr3) ;
		push @cols, $tmp ;
	}
	print join("\t", @cols), "\n" ;
}
