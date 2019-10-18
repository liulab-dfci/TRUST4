#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: ./trust-simplerepo.pl xxx_cdr3.out [--junction trust_annot.fa] > trust_report.out\n" if ( @ARGV == 0 ) ;
my %cdr3 ;  

# Copied from http://www.wellho.net/resources/ex.php4?item=p212/3to3
my %DnaToAa = (
		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
		'TCG' => 'S',    # Serine
		'TCT' => 'S',    # Serine
		'TTC' => 'F',    # Phenylalanine
		'TTT' => 'F',    # Phenylalanine
		'TTA' => 'L',    # Leucine
		'TTG' => 'L',    # Leucine
		'TAC' => 'Y',    # Tyrosine
		'TAT' => 'Y',    # Tyrosine
		'TAA' => '_',    # Stop
		'TAG' => '_',    # Stop
		'TGC' => 'C',    # Cysteine
		'TGT' => 'C',    # Cysteine
		'TGA' => '_',    # Stop
		'TGG' => 'W',    # Tryptophan
		'CTA' => 'L',    # Leucine
		'CTC' => 'L',    # Leucine
		'CTG' => 'L',    # Leucine
		'CTT' => 'L',    # Leucine
		'CCA' => 'P',    # Proline
		'CCC' => 'P',    # Proline
		'CCG' => 'P',    # Proline
		'CCT' => 'P',    # Proline
		'CAC' => 'H',    # Histidine
		'CAT' => 'H',    # Histidine
		'CAA' => 'Q',    # Glutamine
		'CAG' => 'Q',    # Glutamine
		'CGA' => 'R',    # Arginine
		'CGC' => 'R',    # Arginine
		'CGG' => 'R',    # Arginine
		'CGT' => 'R',    # Arginine
		'ATA' => 'I',    # Isoleucine
		'ATC' => 'I',    # Isoleucine
		'ATT' => 'I',    # Isoleucine
		'ATG' => 'M',    # Methionine
		'ACA' => 'T',    # Threonine
		'ACC' => 'T',    # Threonine
		'ACG' => 'T',    # Threonine
		'ACT' => 'T',    # Threonine
		'AAC' => 'N',    # Asparagine
		'AAT' => 'N',    # Asparagine
		'AAA' => 'K',    # Lysine
		'AAG' => 'K',    # Lysine
		'AGC' => 'S',    # Serine
		'AGT' => 'S',    # Serine
		'AGA' => 'R',    # Arginine
		'AGG' => 'R',    # Arginine
		'GTA' => 'V',    # Valine
		'GTC' => 'V',    # Valine
		'GTG' => 'V',    # Valine
		'GTT' => 'V',    # Valine
		'GCA' => 'A',    # Alanine
		'GCC' => 'A',    # Alanine
		'GCG' => 'A',    # Alanine
		'GCT' => 'A',    # Alanine
		'GAC' => 'D',    # Aspartic Acid
		'GAT' => 'D',    # Aspartic Acid
		'GAA' => 'E',    # Glutamic Acid
		'GAG' => 'E',    # Glutamic Acid
		'GGA' => 'G',    # Glycine
		'GGC' => 'G',    # Glycine
		'GGG' => 'G',    # Glycine
		'GGT' => 'G',    # Glycine
		);

sub GetChainType
{
	foreach my $g (@_)
	{
		if ( $g =~ /^IGH/ )
		{
			return 0 ;
		}
		elsif ( $g =~ /^IG/ )
		{
			return 0 ;
		}
		elsif ( $g =~ /^TR/ )
		{
			return 2 ;
		}
	}
	return -1 ;
}

my $i ;
my $annotFile = "" ;
for ( $i = 1 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] == "--junction" )
	{
		$annotFile = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown option ", $ARGV[$i], "\n" ;
	}
}

my %junctionInfo ; # V-end, V-del, VD(J)-Ins, D-Ldel, D-start, D-end, D-Rdel, DJ-Ins, J-del, J-start


# read in the input
open FP1, $ARGV[0] ;
my @totalCnt = (0, 0, 0) ;
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	my $key = join( "\t", ( $cols[2], $cols[3], $cols[4], $cols[5], $cols[8] ) ) ;
	if ( defined $cdr3{ $key } )
	{
		my $val = \@{ $cdr3{ $key } } ;
		$val->[0] = $cols[9] if ( $cols[9] > $val->[0] ) ;
		$val->[1] += $cols[10] ;
	}
	else
	{
		@{ $cdr3{ $key } } = ( $cols[9], $cols[10] ) ;
	}
	my $type = GetChainType( $cols[2], $cols[4], $cols[5] ) ;
	$totalCnt[ $type ] += $cols[10] if ( $type != -1 ) ;
}
close FP1 ;

# Output what we collected.
print( "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\n" ) ;
foreach my $key ( sort { $cdr3{$b}[1] <=> $cdr3{$a}[1] } keys %cdr3 )
{
	my @val = @{ $cdr3{ $key } } ;
	my @info = split /\t/, $key ; # V, D, J, C. CDR3
	my $aa = "" ;
	if ( $val[0] == 0 )
	{
		$aa = "partial" ;
	}
	else
	{
		if ( length( $info[4] ) % 3 != 0 )
		{
			$aa = "out_of_frame" ;
		}
		else
		{
			my $len = length( $info[4] ) ;
			my $s = uc( $info[4] ) ;
			for ( my $i = 0 ; $i < $len ; $i += 3 )
			{
				if ( !defined $DnaToAa{ substr( $s, $i, 3 ) } )
				{	
					$aa = "?" ;
				}
				else
				{
					$aa .= $DnaToAa{ substr( $s, $i, 3 ) } ;
				}
			}
		}
	}
	my $freq = 0 ;
	my $type = GetChainType( $info[0], $info[2], $info[3] ) ;
	$freq = $val[1] / $totalCnt[ $type ] if ( $type != -1 ) ;
	printf( "%.2f\t%e\t%s\t%s\t%s\t%s\t%s\t%s\n", $val[1], $freq, $info[4], $aa, $info[0], $info[1], $info[2], $info[3] ) ;
}
