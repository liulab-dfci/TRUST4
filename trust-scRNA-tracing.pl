#!/usr/bin/env perl

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename ;

die "Usage: ./trust-barcoderep.pl trust_barcode_report.tsv [OPTIONS]\n". 
	"OPTIONS:\n".
	"\t-o STRING: output file prefix (default: trust)\n".
	"\t-s FLOAT: similarity for two CDR3s to be clustered (default: 0.95)\n" 
	if ( @ARGV == 0 ) ;


sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $WD = dirname( abs_path($0) ) ;

my $similarity = 0.95 ;
my $prefix = "trust" ;
my $i ;
for ( $i = 1 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "-o" )
	{
		$prefix = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i + 1] eq "-s" )
	{
		$similarity = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}


# Convert the barcode report file to cdr3 format
open FP1, $ARGV[0] ;
open FPcdr3, ">${prefix}_scrna_igh_cdr3.out" ;
while (<FP1>)
{
	next if ( /^#/ ) ;
	chomp ;
	my @cols = split ;
	next if ( $cols[1] ne "B" || $cols[2] eq "*") ;
	my @igh = split /,/, $cols[2] ;
	
	print FPcdr3 join( "\t", ( $cols[0], 0, $igh[0], $igh[1], $igh[2], $igh[3], "*", "*", $igh[4], 1.00, $igh[6], $igh[7] ) ), "\n" ; 
}
close FP1 ;
close FPcdr3 ;

# Cluster the cells
system_call( "python $WD/trust-cluster.py ${prefix}_scrna_igh_cdr3.out -s $similarity > ${prefix}_scrna_igh_cluster.out" ) ;

# Run evolutionary analysis on this file
system_call( "$WD/clone-evo ${prefix}_scrna_igh_cluster.out -a 0 > ${prefix}_scrna_igh_newick.out" ) ;
system_call( "python $WD/trust-reformatNewick.py ${prefix}_scrna_igh_newick.out > ${prefix}_scrna_igh_tree_cluster.out" ) ;

# Rename the cluster id to barcodes
open FP1, "${prefix}_scrna_igh_cluster.out" ;
my %clusterSubIdToBarcode ;
while (<FP1>)
{
	next if (/^#/) ;
	chomp ;
	my @cols = split ;
	my $key = $cols[0]." ".$cols[1] ;
	$clusterSubIdToBarcode{$key} = $cols[12] ;
}
close FP1 ;

open FP1, "${prefix}_scrna_igh_tree_cluster.out" ;
open FP2, ">${prefix}_scrna_igh_tree.out" ;
while (<FP1>)
{
	if (/^#/)
	{
		print FP2 $_ ;
		next ;
	}
	
	my @cols = split ;
	my $key = $cols[0]." ".$cols[1] ;
	my $parentKey = $cols[3]." ".$cols[4] ;
	
	print FP2 join("\t", ($cols[0], $clusterSubIdToBarcode{$key}, $cols[2], $cols[3], $clusterSubIdToBarcode{$parentKey}, $cols[5])), "\n" ;
}
close FP1 ;
close FP2 ;

print STDERR "Done. Final result is in ${prefix}_scrna_igh_tree.out\n" ;
