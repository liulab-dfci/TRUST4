#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: trust-airr.pl trust_cdr3.out trust_annot.fa > trust_airr.tsv\n" if (@ARGV == 0) ;

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

sub Translate
{
	my $s = $_[0] ;
	my $aa = "" ;
	if (length($s) % 3 != 0) 
	{
		return $aa ;
	}
	else 
	{
		my $len = length($s) ;
		my $s = uc($s) ;
		for (my $i = 0 ; $i < $len ; $i += 3)
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

		return "" if ($aa =~ /\?/ || $aa =~ /_/) ;
		return $aa ;
	}
	return $aa ;
}

sub CoordToCigar
{
	my $cigar = "" ;
	$cigar = $_[3]."S" if ($_[3] > 0) ;
	$cigar .= ($_[4]-$_[3]+1)."M" ;
	$cigar .= ($_[0]-$_[4]-1)."S" if ($_[0]-$_[4]-1 > 0) ;
} 

my %seqCDR3s ;

open FP, $ARGV[0] ;
while (<FP>)
{
	chomp ;
	my @cols = split ;
	next if ($cols[9] <= 0) ;
	push @{$seqCDR3s{$cols[0]}}, $cols[8] ; # cdr3nt
	push @{$seqCDR3s{$cols[0]}}, $cols[10] ; # abundance
	my $fullLength = "F" ;
	$fullLength = "T" if ($cols[12] == 1) ;
	push @{$seqCDR3s{$cols[0]}}, $fullLength ; # full-length
}
close FP ;

print "sequence_id\tsequence\trev_comp\tproductive\tv_call\td_call\tj_call\tc_call\tsequence_alignment\tgermline_alignment\tjunction\tjunction_aa\tv_cigar\td_cigar\tj_cigar\tv_identity\tj_identity\tcomplete_vdj\tconsensus_count\n" ;

open FP, $ARGV[1] ;
while (<FP>)
{
	chomp ;
	my $header = $_ ;
	my $seq = <FP> ;
	chomp $seq ;
	
	my @cols = split /\s/, substr($header,1) ;
	
	my @vCoord ;
	my @dCoord ;
	my @jCoord ;
	my @cCoord ;
	my @cdr3Coord ;
	my $cdr3 ;
	my $vcall = "" ;
	my $vcigar = "" ;
	my $videntity = "" ;
	my $dcall = "" ;
	my $dcigar = "" ;
	my $jcall = "" ; 
	my $jcigar = "" ;
	my $jidentity = "" ;
	my $ccall = "" ;

	if ( $cols[3] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		$vcall = (split /\(/, $cols[3])[0] ;
		@vCoord = ($1, $2, $3, $4, $5) ;
		$vcigar = CoordToCigar(@vCoord) ;
		$videntity = (split /:/, $cols[3])[-1] ;
	}
	else
	{
		@vCoord = (-1, -1, -1, -1, -1) ;
	}

	if ( $cols[4] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		$dcall = (split /\(/, $cols[4])[0] ;
		@dCoord = ($1, $2, $3, $4, $5) ;
		$dcigar = CoordToCigar(@dCoord) ;
	}
	else
	{
		@dCoord = (-1, -1, -1, -1, -1) ;
	}

	if ( $cols[5] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		$jcall = (split /\(/, $cols[5])[0] ;
		@jCoord = ($1, $2, $3, $4, $5) ;
		$jcigar = CoordToCigar(@jCoord) ;
		$jidentity = (split /:/, $cols[3])[-1] ;
	}
	else
	{
		@jCoord = (-1, -1, -1, -1, -1) ;
	}

	if ( $cols[6] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\)/ )
	{
		$ccall = (split /\(/, $cols[6])[0] ;
		@cCoord = ($1, $2, $3, $4, $5) ;
	}
	else
	{
		@cCoord = (-1, -1, -1, -1, -1) ;
	}

	next if ( $cols[9] =~ /:0\.00/ ) ; # partial CDR3

	if ( $cols[9] =~ /CDR3\(([0-9]+?)-([0-9]+?)\)/ )
	{
		@cdr3Coord = ($1, $2) ;		
	}
	else
	{
		next ;
	}
	$cdr3 = ( split /=/, $cols[9] )[1] ;

	my @cdr3s = @{$seqCDR3s{$cols[0]}} ;

#print "sequence_id\tsequence\trev_comp\tproductive\tv_call\td_call\tj_call\tc_call\tsequence_alignment\tgermline_alignment\tjunction\tjunction_aa\tv_cigar\td_cigar\tj_cigar\tv_identity\tj_identity\tcomplete_vdj\tconsensus_count\n" ;
	for (my $i = 0 ; $i < scalar(@cdr3s) ; $i += 3)
	{
		my $cdr3aa = Translate($cdr3s[$i]) ;
		my $productive = "T" ;
		$productive = "F" if ($cdr3aa eq "") ;
		my $outputSeq = $seq ;
		substr($outputSeq, $cdr3Coord[0], $cdr3Coord[1] - $cdr3Coord[0] + 1, $cdr3s[$i]) ;
		print join("\t", ($cols[0]."_".($i/3), $outputSeq, "F", $productive,
			$vcall, $dcall, $jcall, $ccall, "", "", $cdr3s[$i], $cdr3aa, 
		  $vcigar, $dcigar, $jcigar, $videntity, $jidentity, $cdr3s[$i + 2], $cdr3s[$i + 1]) ). "\n" ;
	}
}
close FP ;
