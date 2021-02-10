#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: ./trust-simplerep.pl xxx_cdr3.out [OPTIONS] > trust_report.out\n". 
	"OPTIONS:\n".
	"\t--decimalCnt: the count column uses decimal instead of truncated integer. (default: not used)\n".
	"\t--barcodeCnt: the count column is the number of barcode instead of read support. (default: not used)\n".
	"\t--junction trust_annot.fa: output junction information for the CDR3 (default: not used)\n".
	"\t--reportPartial: include partial CDR3 in report. (default: no partial)\n".
	"\t--filterTcrError FLOAT: filter TCR CDR3s less than the fraction of representative CDR3 in the consensus. (default: 0.05)\n".
	"\t--filterBcrError FLOAT: filter BCR CDR3s less than the fraction of representative CDR3 in the consensus. (default: 0)\n"
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
			return 1 ;
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
	
	if ($_[2] ne ".")
	{
		$ret = (split /\*/, $ret)[0] ; # Remove the allele id from c gene
		for ( $i = 0 ; $i <= 1 ; ++$i )
		{
			next if ( $_[$i] eq "." ) ;
			if ( !($_[$i] =~ /^IGH/) )
			{
				$ret = substr($ret, 0, 4 );
				last ;
			}
		}
		
		return $ret ;
	}
		
	# For TRA and TRD gene, we don't infer its constant gene.
	if ($_[0] =~ /^TR[AD]/ || $_[1] eq ".")
	{
		return $ret ;
	}

	for ( $i = 1 ; $i >= 0 ; --$i )
	{
		next if ( $_[$i] eq "." ) ;
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
my $annotFile = "" ;
my $reportJunctionInfo = 0 ;
my $tcrErrorFilter = 0.05 ;
my $bcrErrorFilter = 0 ;
my $roundDownCount = 1 ;
my $useBarcodeCnt = 0 ;
my $reportPartial = 0 ;
for ( $i = 1 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "--junction" )
	{
		$annotFile = $ARGV[$i + 1] ;
		$reportJunctionInfo = 1 ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--filterTcrError" )
	{
		$tcrErrorFilter = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--filterBcrError" )
	{
		$bcrErrorFilter = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--decimalCnt" )
	{
		$roundDownCount = 0 ;
	}
	elsif ( $ARGV[$i] eq "--barcodeCnt" )
	{
		$useBarcodeCnt = 1 ;
	}
	elsif ( $ARGV[$i] eq "--reportPartial" )
	{
		$reportPartial = 1 ;
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

# collect the read count for each assembly id.
open FP1, $ARGV[0] ;
my %assemblyMostReads ; 
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	my $assemblyId = $cols[0] ;
	if ( defined $assemblyMostReads{ $assemblyId } )
	{
		$assemblyMostReads{$assemblyId} = $cols[10] if ($cols[10] > $assemblyMostReads{ $assemblyId });
	}
	else
	{
		$assemblyMostReads{$assemblyId} = $cols[10] ;
	}
}
close FP1 ;

# Read in the input
open FP1, $ARGV[0] ;
my @totalCnt = (0, 0, 0) ;
my %cdr3AssemblyId ; # Record which assembly is representative for the cdr3
my %cdr3Barcode ; # record whether a "chain_cdr3_barcode" has showed up before or not 
my %assemblyIdFullLength ; # record whether this assembly id contains full length assembly
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	next if ( $reportPartial == 0 && $cols[9] == 0 ) ;
	
	my $assemblyId = $cols[0] ;
	
	for (my $i = 2 ; $i <= 5 ; ++$i )
	{
		if ($cols[$i] eq "*")
		{
			$cols[$i] = "." ;
		}
	}

	my $vgene = (split /,/, $cols[2])[0] ;
	my $dgene = (split /,/, $cols[3])[0] ;
	my $jgene = (split /,/, $cols[4])[0] ;
	my $cgene = (split /,/, $cols[5])[0] ;
	
	$cgene = InferConstantGene( $vgene, $jgene, $cgene ) ;
	my $key = join( "\t", ( $vgene, $dgene, $jgene, $cgene, $cols[8] ) ) ;
	my $type = GetDetailChainType( $vgene, $jgene, $cgene ) ;
	
	if ($type > 2) # TCR, ignore the low abundance 
	{
		if ($cols[10] < $assemblyMostReads{$assemblyId} * $tcrErrorFilter ) 
		{
			next ;
		}
	}
	elsif ( $type <= 2 ) # BCR
	{
		if ($cols[10] < $assemblyMostReads{$assemblyId} * $bcrErrorFilter ) 
		{
			next ;
		}
		$type = 1 if ($type == 2) ; # Merge IGL to IGK.
	}
	
	# Ignore the CDR3 that is too long, could be from mis-annotation.
	next if (length($cols[8]) >= 180 ) ;

	if ( $useBarcodeCnt )
	{
		my @cols2 = split/_/, $assemblyId ;
		my $barcode = join( "_", @cols2[0..scalar(@cols2)-2] ) ;
		my $tmp = $type."_".$cols[8]."_".$barcode ;
		if ( defined $cdr3Barcode{$tmp} )
		{
			next ;
		}
		else
		{
			$cdr3Barcode{$tmp} = 1 ;
			$cols[10] = 1 ;
		}
	}

	if ( defined $cdr3{ $key } )
	{
		my $val = \@{ $cdr3{ $key } } ;
		if ( $cols[9] > $val->[0] )  # found an assembly with better CDR3 score
		{
			$val->[0] = $cols[9] ;
		}

		if ( $cols[10] > $val->[3] ) 
		{
			$val->[2] = $assemblyId ;
			$val->[3] = $cols[10];
		}
		$val->[1] += $cols[10] ;
	}
	else
	{
		@{ $cdr3{ $key } } = ( $cols[9], $cols[10], $assemblyId, $cols[10] ) ;
	}
	$totalCnt[ $type ] += $cols[10] if ( $type != -1 ) ;

	$assemblyIdFullLength{$assemblyId} = $cols[12] ;
}
close FP1 ;

# Output what we collected.
print( "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length" ) ;
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
					$aa .= "?" ;
				}
				else
				{
					$aa .= $DnaToAa{ substr( $s, $i, 3 ) } ;
				}
			}
		}
	}
	my $freq = 0 ;
	my $type = GetDetailChainType( $info[0], $info[2], $info[3] ) ;
	$type = 1 if ($type == 2) ;
	$freq = $val[1] / $totalCnt[ $type ] if ( $type != -1 ) ;
	if ($roundDownCount == 1 )
	{
		my $cnt = int($val[1]) ;
		next if ($cnt == 0) ;
		printf( "%d\t%e\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d", $cnt, $freq, $info[4], $aa, 
					$info[0], $info[1], $info[2], $info[3], $val[2], $assemblyIdFullLength{$val[2]} ) ;
	}
	else
	{
		printf( "%.2f\t%e\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d", $val[1], $freq, $info[4], $aa, 
					$info[0], $info[1], $info[2], $info[3], $val[2], $assemblyIdFullLength{$val[2]} ) ;
	}
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
