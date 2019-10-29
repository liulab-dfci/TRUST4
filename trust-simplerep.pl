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
my $reportJunctionInfo = 0 ;
for ( $i = 1 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "--junction" )
	{
		$annotFile = $ARGV[$i + 1] ;
		$reportJunctionInfo = 1 ;
		++$i ;
	}
	else
	{
		die "Unknown option ", $ARGV[$i], "\n" ;
	}
}

my %junctionInfo ; # key: assembly id. Values: V-end, V-del, VD(J)-Ins, D-Ldel, D-start, D-end, D-Rdel, DJ-Ins, J-del, J-start

if ( $reportJunctionInfo == 1 )
{
	open FP1, $annotFile ;
	while ( <FP1> )
	{
		next if ( !/^>/ ) ;
		my @cols = split /\s/, $_ ;

		my @vCoord ;
		my @dCoord ;
		my @jCoord ;
		my @cdr3Coord ;

		if ( $cols[3] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
		{
			@vCoord = ($1, $2, $3, $4, $5) ;
		}
		else
		{
			#die "Wrong format $header\n" ;
			next ;
		}
		if ( $cols[4] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
		{
			@dCoord = ($1, $2, $3, $4, $5) ;
		}
		else
		{
			@dCoord = (-1, -1, -1, -1, -1) ;
		}
		if ( $cols[5] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
		{
			@jCoord = ($1, $2, $3, $4, $5) ;
		}
		else
		{
			#die "Wrong format $header\n" ;
			next ;
		}
		
		next if ( $vCoord[2] >= $jCoord[1] ) ;		
		next if ( $cols[9] =~ /:0\.00/ ) ;
		if ( $cols[9] =~ /CDR3\(([0-9]+?)-([0-9]+?)\)/ )
		{
			@cdr3Coord = ($1, $2) ;		
		}
		else
		{
			next ;
		}
		
		next if ( $vCoord[2] < $cdr3Coord[0] || $jCoord[1] > $cdr3Coord[1] ) ;

		my $chainType = substr( $cols[3], 0, 3 ) ;
		my @info ;
		# v end, vdelete 
		push @info, $vCoord[2] - $cdr3Coord[0], $vCoord[0] - $vCoord[4] - 1 ;
		
		# v insert, d-left-delete, d-start, d-end, d-right-delete, j insert
		if ( $chainType eq "IGH" || $chainType eq "TRB" || $chainType eq "TRD" )
		{
			if ( $dCoord[0] == -1 || $dCoord[1] <= $vCoord[2] || $dCoord[2] >= $jCoord[1])
			{
				push @info, "*", "*", "*", "*", "*", "*" ;
			}
			else
			{
				push @info, $dCoord[1] - $vCoord[2] - 1, $dCoord[3], 
						$dCoord[1] - $cdr3Coord[0], $dCoord[2] - $cdr3Coord[0], 
						$dCoord[0] - $dCoord[4] - 1,
						$jCoord[1] - $dCoord[2] - 1 ;
			}
		}
		else
		{
			# already make sure J gene is after V gene in this case.
			push @info, $jCoord[1] - $vCoord[2] - 1, "*", "*", "*", "*", "*" ; 
		}

		# j-left-delete, j-start
		push @info, $jCoord[3], $jCoord[1] - $cdr3Coord[0] ;
		
		#print(substr($cols[0], 1), ":", join(",", @info), "\n") ;
		@{$junctionInfo{ substr( $cols[0], 1 )}} = @info ;
	}
	close FP1 ;
}
# read in the input
open FP1, $ARGV[0] ;
my @totalCnt = (0, 0, 0) ;
my %cdr3AssemblyId ; # Record which assembly is representative for the cdr3
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	my $key = join( "\t", ( $cols[2], $cols[3], $cols[4], $cols[5], $cols[8] ) ) ;
	my $assemblyId = $cols[0] ;
	if ( defined $cdr3{ $key } )
	{
		my $val = \@{ $cdr3{ $key } } ;
		if ( $cols[9] > $val->[0] ) 
		{
			$val->[0] = $cols[9];
			$val->[2] = $assemblyId ;
		}
		$val->[1] += $cols[10] ;
	}
	else
	{
		@{ $cdr3{ $key } } = ( $cols[9], $cols[10], $assemblyId ) ;
	}
	my $type = GetChainType( $cols[2], $cols[4], $cols[5] ) ;
	$totalCnt[ $type ] += $cols[10] if ( $type != -1 ) ;
}
close FP1 ;

# Output what we collected.
print( "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC" ) ;
if ( $reportJunctionInfo == 1 )
{
	print( "\tjunction" )
}
print( "\n" ) ;

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
	printf( "%.2f\t%e\t%s\t%s\t%s\t%s\t%s\t%s", $val[1], $freq, $info[4], $aa, $info[0], $info[1], $info[2], $info[3] ) ;
	if ( $reportJunctionInfo == 1 )
	{
	        if ( defined $junctionInfo{$val[2]} )
		{
			print( "\t", join( ",", @{$junctionInfo{$val[2]}}) ) ;
		}
		else
		{
			print( "\t*" ) ;
		}
	}
	print( "\n" ) ;
}
