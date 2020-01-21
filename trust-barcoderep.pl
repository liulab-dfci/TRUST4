#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: ./trust-barcoderep.pl xxx_cdr3.out [OPTIONS] > trust_barcode_report.tsv\n". 
	"OPTIONS:\n".
	"\t--noPartial: do not including partial CDR3 in report. (default: include partial)\n"
	if ( @ARGV == 0 ) ;

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

sub GetDetailChainType
{
	foreach my $g (@_)
	{
		if ( $g =~ /^IGH/ )
		{
			return 0 ;
		}
		elsif ( $g =~ /^IGK/ )
		{
			return 1 ;
		}
		elsif ( $g =~ /^IGL/ )
		{
			return 2 ;
		}
		elsif ( $g =~ /^TRA/ )
		{
			return 3 ;
		}
		elsif ( $g =~ /^TRB/ )
		{
			return 4 ;
		}
		elsif ( $g =~ /^TRG/ )
		{
			return 5 ;
		}
		elsif ( $g =~ /^TRD/ )
		{
			return 6 ;
		}
	}
	
}

# Use input V, J, C gene to report back the C gene if it is missing.
sub InferConstantGene
{
	my $ret = $_[2] ;
	my $i ;
	
	if ($_[2] ne "*")
	{
		for ( $i = 0 ; $i <= 1 ; ++$i )
		{
			next if ( $_[$i] eq "*" ) ;
			if ( !($_[$i] =~ /^IGH/) )
			{
				$ret = substr($ret, 0, 4 );
				last ;
			}
		}
		
		return $ret ;
	}
	
	for ( $i = 0 ; $i <= 1 ; ++$i )
	{
		next if ( $_[$i] eq "*" ) ;
		if ( $_[$i] =~ /^IGH/ )
		{
			return $ret ;
		}
		my $prefix = substr( $_[$i], 0, 3 ) ;
		
		return $prefix."C" ; 
	}
	return $ret ;
}

my $i ;
my $reportPartial = 1 ;
for ( $i = 1 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "--noPartial" )
	{
		$reportPartial = 0 ;
	}
	else
	{
		die "Unknown option ", $ARGV[$i], "\n" ;
	}
}


# collect the read count for each chain from assembly id.
open FP1, $ARGV[0] ;
my %barcodeChainAbund ; 
my %barcodeChainInfo ;
my %barcodeShownup ;
my @barcodeList ;

# Read in the report, store the information for each barcode.
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	next if ( $reportPartial == 0 && $cols[9] == 0 ) ;

	my $assemblyId = $cols[0] ;
	my $vgene = (split /,/, $cols[2])[0] ;
	my $dgene = (split /,/, $cols[3])[0] ;
	my $jgene = (split /,/, $cols[4])[0] ;
	my $cgene = (split /,/, $cols[5])[0] ;
	$cgene = InferConstantGene( $vgene, $jgene, $cgene ) ;

	my @cols2 = split/_/, $assemblyId ;
	my $barcode = join( "_", @cols2[0..scalar(@cols2)-2] ) ;
	my $key = $barcode."_".GetDetailChainType( $vgene, $jgene, $cgene ) ;
	my $aa ;
	
	if ( !defined $barcodeShownup{ $barcode } )
	{
		$barcodeShownup{ $barcode } = 1 ;
		push @barcodeList, $barcode ;
	}

	if ( $cols[9] == 0 )
	{
		$aa = "partial" ;
	}
	else
	{
		if ( length( $cols[8] ) % 3 != 0 )
		{
			$aa = "out_of_frame" ;
		}
		else
		{
			my $len = length( $cols[8] ) ;
			my $s = uc( $cols[8] ) ;
			for ( my $i = 0 ; $i < $len ; $i += 3 )
			{
				if ( !defined $DnaToAa{ substr( $s, $i, 3 ) } )
				{	
					$aa .= "?" ;
				}
				else
				{
					$aa .= $DnaToAa{ substr( $s, $i, 3 ) } ;
				}
			}
		}
	}

	if ( defined $barcodeChainAbund{ $key } )
	{
		if ( $cols[10] > $barcodeChainAbund{ $key } )
		{
			$barcodeChainAbund{ $key } = $cols[10] ;
			$barcodeChainInfo{ $key } = join( ",", ($vgene, $dgene, $jgene, $cgene, $cols[8], $aa) ) ;
		}
	}
	else
	{
		$barcodeChainAbund{ $key } = $cols[10] ;
		$barcodeChainInfo{ $key } = join( ",", ($vgene, $dgene, $jgene, $cgene, $cols[8], $aa) ) ;
	}
}
close FP1 ;

# Output what we collected.
print( "#barcode\tcell_type\tchain1\tchain2\n" ) ;

foreach my $barcode (@barcodeList )
{
	# Determine type
	my $i ;
	my $mainType ; # 0-IG, 1-TRA/B, 2-TRG/D
	my $cellType ;
	my $max = -1 ;
	my $maxTag = -1 ;
	my $chain1 = "*" ;
	my $chain2 = "*" ;
	for ( $i = 0 ; $i < 7 ; ++$i )
	{
		my $key = $barcode."_".$i ;
		if ( defined $barcodeChainAbund{ $key } && $barcodeChainAbund{ $key } > $max )
		{
			$max = $barcodeChainAbund{ $key } ;
			$maxTag = $i ;
		}
	}
	if ( $maxTag <= 2 )
	{
		$mainType = 0 ;
		my $keyH = $barcode."_0" ;
		my $keyK = $barcode."_1" ;
		my $keyL = $barcode."_2" ;
		$chain1 = $barcodeChainInfo{ $keyH } if ( defined $barcodeChainInfo{ $keyH } ) ;

		if ( defined $barcodeChainInfo{ $keyK } && defined $barcodeChainInfo{ $keyL } )
		{
			if ( $barcodeChainAbund{ $keyK } >= $barcodeChainAbund{ $keyL } )
			{
				$chain2 = $barcodeChainInfo{ $keyK } ;
			}
			else
			{
				$chain2 = $barcodeChainInfo{ $keyL } ;
			}
		}
		elsif ( defined $barcodeChainInfo{ $keyK } )
		{
			$chain2 = $barcodeChainInfo{ $keyK } ;
		}
		elsif ( defined $barcodeChainInfo{ $keyL } )
		{
			$chain2 = $barcodeChainInfo{ $keyL } ;
		}

		$cellType = "B" ;
	}
	else
	{
		my $key1 ;
		my $key2 ;
		if ( $maxTag <= 4 )
		{
			$key1 = $barcode."_3" ;
			$key2 = $barcode."_4" ;
			$cellType = "abT" ;
		}
		elsif ( $maxTag <= 6 )
		{
			$key1 = $barcode."_5" ;
			$key2 = $barcode."_6" ;
			$cellType = "gdT" ;
		}
		$chain1 = $barcodeChainInfo{ $key1 } if ( defined $barcodeChainInfo{ $key1 } ) ;
		$chain2 = $barcodeChainInfo{ $key2 } if ( defined $barcodeChainInfo{ $key2 } ) ;
	}
	print( join( "\t", ($barcode, $cellType, $chain1, $chain2 ) ), "\n" ) ;
}
