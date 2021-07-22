#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: ./trust-barcoderep-to-10X.pl trust_barcode_report.tsv 10X_report_prefix\n" if (@ARGV == 0) ;

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
	return 7 ;	
}

sub IsProductive
{
	my $aa = $_[0] ;
	return 0 if ($aa eq "partial" || $aa =~ /_/ || $aa =~ /\?/) ;
	return 1 ;
}

my @chainName = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD", "None") ;
open FP, $ARGV[0] ;
my $prefix = $ARGV[1] ;
open FPoutT, ">".$prefix."_t.csv" ;
open FPoutB, ">".$prefix."_b.csv" ;

print FPoutT "barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id\n" ;
print FPoutB "barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id\n" ;

my $header = <FP> ;
while (<FP>)
{
	chomp ;
	my @cols = split ;
        #AACTCTTGTTATCCGA-1	abT	*	TRAV6*01,*,TRAJ43*01,*,TGTGCTCTAGCCGGGGAGGGCATGCGCTTT,CALAGEGMRF,2.80,AACTCTTGTTATCCGA-1_26095,100.00	*	*
	for (my $i = 2 ; $i <= 3 ; ++$i)
	{
		next if ($cols[$i] eq "*") ;
		my @chainCols = split /,/, $cols[$i] ;
		#barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id
		#AAACCTGAGTCGATAA-1,True,AAACCTGAGTCGATAA-1_contig_1,True,551,IGK,IGKV1-5,None,IGKJ1,IGKC,True,True,CQQYNSYSRTF,TGCCAACAGTATAATAGTTATTCTCGAACGTTC,1197,22,clonotype77,clonotype77_consensus_2
		my @outputCols = ($cols[0], "True", $chainCols[7], "True", "None", 
			$chainName[GetDetailChainType($chainCols[0], $chainCols[2],$chainCols[3])],
			$chainCols[0] eq "*" ? "None": $chainCols[0],
			$chainCols[1] eq "*" ? "None": $chainCols[1],
			$chainCols[2] eq "*" ? "None": $chainCols[2],
			$chainCols[3] eq "*" ? "None": $chainCols[3],
			$chainCols[9] == 1 ? "True" : "False", IsProductive($chainCols[5]) ? "True" : "False", 
			$chainCols[5], $chainCols[4], $chainCols[6], $chainCols[6], "None", "None"
			) ;
		if (substr($cols[1], -1) eq "T")
		{
			print FPoutT join(",", @outputCols), "\n" ;
		}
		else
		{
			print FPoutB join(",", @outputCols), "\n" ;
		}
	}
}
close FPoutT ;
close FPoutB ;
close FP ;
