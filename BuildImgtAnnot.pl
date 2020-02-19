#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl species_name(Homo_sapiens or others) > output.fa\n" if ( @ARGV == 0 ) ;

sub system_call
{
	print STDERR "SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $species = $ARGV[0] ;
# Download
#system_call( "wget -np -nd -r -A fasta -P tmp_download http://www.imgt.org//download/V-QUEST/IMGT_V-QUEST_reference_directory/".$species."/" ) ;
#system_call( "cat tmp_download/*.fasta > tmp_download/vquest.fa") ;
system_call( "wget -O IMGT_download.fa http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP" ) ;

# Reformat
#>J00256|IGHJ1*01|Homo sapiens|F|J-REGION|723..774|52 nt|1| | | | |52+0=52| | |
#
open FP1, "IMGT_download.fa" ;
my $prevId = "" ;
my $prevGeneName = "" ;
my $output = 0 ; 
my $skipHeader = 0 ; # For the genes that are split into pieces in the database. 

while ( <FP1> )
{
	if ( /^>/ )
	{
		my @cols = split /\|/, substr( $_, 1 ) ;
		my $s = $cols[2] ;
		$s =~ tr/ /_/ ;
		if ( !($s =~ /$species/) )
		{
			$output = 0 ;
		}
		elsif ( !( $cols[1] =~ /^IG/ ) && !($cols[1] =~ /^TR/) )
		{
			$output = 0 ;
		}
		else
		{
			$output = 1 ;
			if ( $cols[1] eq $prevGeneName )
			{
				$output = 0 if ( $cols[0] ne $prevId ) ; # Could be some middle part of a gene.
				$skipHeader = 1 ;
			}
			else
			{
				$skipHeader = 0 ;
			}
		}
		$prevId = $cols[0] ;
		$prevGeneName = $cols[1] ;
	}
	next if ( $output == 0) ;
	if ( !/^>/ )
	{
		$_ =~ tr/acgtn/ACGTN/ ;
		print $_ ;
	}
	elsif ( $skipHeader == 0 )
	{
		print ">".(split /\|/ )[1]."\n" ;
	}
}
close FP1 ;

unlink "IMGT_download.fa" ;
