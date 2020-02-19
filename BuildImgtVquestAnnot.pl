#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl bcrtcr.fa cgene.list species_name(Homo_sapiens or others)\n" if ( @ARGV == 0 ) ;

sub system_call
{
	print STDERR "SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $species = $ARGV[2] ;
# Download
system_call( "wget -np -nd -r -A fasta -P tmp_download http://www.imgt.org//download/V-QUEST/IMGT_V-QUEST_reference_directory/".$species."/" ) ;
system_call( "cat tmp_download/*.fasta > tmp_download/vquest.fa") ;

# Reformat
#>J00256|IGHJ1*01|Homo sapiens|F|J-REGION|723..774|52 nt|1| | | | |52+0=52| | |
#
open FP1, "tmp_download/vquest.fa" ;
open FPout, ">IMGT+C.fa" ;
while ( <FP1> )
{
	if ( !/^>/ )
	{
		$_ =~ tr/acgtn/ACGTN/ ;
		print FPout $_ ;
		next ;
	}
	print FPout ">".(split /\|/ )[1]."\n" ;
}
close FP1 ;
close FPout ;

# Append C gene from bcrtcr.fa
system_call( "grep -A 1 --no-group-separator -f ".$ARGV[1]." ".$ARGV[0]." >> IMGT+C.fa" ) ;


# remove the temporary download directory
rmdir "tmp_download"
