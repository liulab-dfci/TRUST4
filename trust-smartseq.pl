#!/usr/bin/env perl

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename ;

die "TRUST4 SMART-seq pipeline usage: perl smartseq-process.pl [OPTIONS]:\n".
    "\t-1 STRING: file containing the list of read 1 (or single-end) files\n".
    "\t-2 STRING: file containing the list of read 2 files.\n".
    "\t-f STRING: path to the fasta file coordinate and sequence of V/D/J/C genes\n".
    "\t--ref STRING: path to detailed V/D/J/C gene reference file from IMGT database. (default: not used but recommended)\n".
    "\t-o STRING: prefix of final output files. (default: TRUST)\n".
    "\t-t INT: number of threads (default: 1)\n".
    "\t--representative INT: number of representative for each detected chain (default: 1)\n".
		"\t--trust-path STRING: TRUST4 executable files (default: same as this script)\n"
		#"\t--noclear: do not clear the intermediate results (default: clear)\n"
    if (@ARGV == 0) ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

sub GetPairChainType
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
			return 1 ;
		}
		elsif ( $g =~ /^TRA/ )
		{
			return 2 ;
		}
		elsif ( $g =~ /^TRB/ )
		{
			return 3 ;
		}
		elsif ( $g =~ /^TRG/ )
		{
			return 4 ;
		}
		elsif ( $g =~ /^TRD/ )
		{
			return 5 ;
		}
	}
}
my $WD = dirname( abs_path( $0 ) ) ;

my $i ;
my $readFile1 = "" ;
my $readFile2 = "" ;
my $outputPrefix = "TRUST" ;
my $hasMate = 0 ;
my $trust4Args = "" ;
my $representativeN = 1 ;
for ( $i = 0 ; $i < @ARGV ; ++$i )
{
	if ($ARGV[$i] eq "-1")	
	{
		$readFile1 = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-2")
	{
		$readFile2 = $ARGV[$i + 1] ;
		++$i ;
		$hasMate = 1 ;
	}
	elsif ($ARGV[$i] eq "-o")
	{
		$outputPrefix = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-t" || $ARGV[$i] eq "-f" || $ARGV[$i] eq "--ref")
	{
		$trust4Args .= " ".$ARGV[$i]." ".$ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "--representative")
	{
		$representativeN = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "--trust-path")
	{
		$WD = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}

die "Need to use -1 to specify the list of read 1 files.\n" if ($readFile1 eq "") ;

# Process the fastq file for each cell one by one
open FP1, $readFile1 ;
my $FP2 ;
open $FP2, $readFile2 if ($hasMate) ;
open FPfinalreport, ">${outputPrefix}_report.tsv" ;
open FPfinalannot, ">${outputPrefix}_annot.fa" ;
open FPfinalairr, ">${outputPrefix}_airr.tsv" ;

print FPfinalreport "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length\n" ;
my $cellProcessed = 0 ;

while (<FP1>)
{
	chomp ;
	my $file1 = $_ ;
	my $file2 = "" ;
	
	my $fname = (fileparse($file1))[0] ;
	my $cellPrefix = (split /\./, $fname)[0] ; # use the content before the first "." as identifier.
	if ($hasMate)
	{
		$file2 = <$FP2> ;
		chomp $file2 ;
		system("$WD/run-trust4 $trust4Args -1 $file1 -2 $file2 -o tmp_smartseq --skipMateExtension") ;
	}
	else
	{
		system("$WD/run-trust4 $trust4Args -u $file1 -o tmp_smartseq") ;
	}

	# Summary the output and extracting the assemblies for representative chains.
	open FPreport, "tmp_smartseq_report.tsv" ;
	my $header = <FPreport>	 ;
	my $line = <FPreport> ; 
	if (!$line) 
	{
		print STDERR "WARNING: no assemblies from $cellPrefix.\n" ;
		next ;
	}

	chomp $line ;
	my $report1 = $line ;
	my @cols = split /\t/, $line ;
	my $mainChainType = GetPairChainType($cols[4], $cols[6], $cols[7]) ;
	
	my $report2 = "" ;

	my @representativeCols ;

	@{$representativeCols[0]} = @cols ;
	my $representativeCols1Cnt = 1 ;
	my $representativeCols2Cnt = 0 ;
	my $representativeCnt = 1 ;

	# Go throught the list to find the second representative.
	while (<FPreport>)	
	{
		chomp ;
		$line = $_ ;
		@cols = split /\t/, $line ;
		my $chainType = GetPairChainType($cols[4], $cols[6], $cols[7]) ;
		my $add = 0 ;
		if ($chainType == $mainChainType)
		{
			if ($representativeCols1Cnt < $representativeN)
			{
				$add = 1 ;
				++$representativeCols1Cnt ;
			}
		}
		elsif ( int($chainType/2) == int($mainChainType/2) && $chainType%2 == 1 - ($mainChainType%2))
		{
			if ($representativeCols2Cnt < $representativeN)
			{
				$add = 1 ;
				++$representativeCols2Cnt ;
			}
		}
		if ($add == 1)
		{
				@{$representativeCols[$representativeCnt]} = @cols ;
				++$representativeCnt ;
		}

		last if ($representativeCols1Cnt >= $representativeN && $representativeCols2Cnt >= $representativeN) ;
	}
	
	# Output the representative information
	my %selectedContigs ;
	for ($i = 0 ; $i < $representativeCols1Cnt + $representativeCols2Cnt ; ++$i)
	{
		my @cols = @{$representativeCols[$i]} ;
		my $contigId = $cols[8] ;
		$cols[8] = $cellPrefix."_".$contigId ;
		$selectedContigs{$contigId} = $i if (!defined $selectedContigs{$contigId}) ;
		print FPfinalreport join("\t", @cols), "\n" ;
	}
	close FPreport ;

	# Process the annotation file for the assemblies
	open FPannot, "tmp_smartseq_annot.fa" ;
	while (<FPannot>)
	{
		chomp ;
		my $header = $_ ;
		my $seq = <FPannot> ;
		chomp $seq ;

		my @cols = split / /, $header ;
		my $contigId = substr($cols[0], 1) ;
		if (defined $selectedContigs{$contigId})
		{
			$cols[0] = ">$cellPrefix"."_".$contigId ;
			print FPfinalannot join(" ", @cols), "\n$seq\n" ;
		}
	}
	close FPannot ;

	# Process the AIRR file
	open FPairr, "tmp_smartseq_airr.tsv" ;
	my $lineCnt = 0 ;
	while (<FPairr>)
	{
		chomp ;
		if ($cellProcessed == 0 && $lineCnt == 0)
		{
			print FPfinalairr $_,"\n" ;
		}
		if ($lineCnt == 0)
		{
			++$lineCnt ;
			next ;
		}
		
		++$lineCnt ;
		my $line = $_ ;
		my @cols = split /\t/, $line ;
		my $contigId = (split /_/, $cols[0])[0] ;
		next if (!defined $selectedContigs{$contigId}) ;
		my @matchedCols = @{$representativeCols[ $selectedContigs{$contigId} ]} ;
		if ($matchedCols[2] eq $cols[12])
		{
			$cols[0] = ${cellPrefix}."_".$cols[0] ;
			print FPfinalairr join("\t", @cols), "\n" ;
		}
	}
	close FPairr ;
	# remove the temporary files.
	system("rm ./tmp_smartseq_*") ;

	++$cellProcessed ;
}
close FP1 ;
close $FP2 if ($hasMate) ;

close FPfinalreport ;
close FPfinalannot ;
close FPfinalairr ;
