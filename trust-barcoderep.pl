#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: ./trust-barcoderep.pl xxx_cdr3.out [OPTIONS] > trust_barcode_report.tsv\n". 
	"OPTIONS:\n".
	"\t-a xxx_annot.fa: TRUST4's annotation file. (default: not used)\n".
	"\t--noImputation: do not perform imputation for partial CDR3. (default: impute)\n".
	"\t--imputeBCR: perform imputation for BCR partial CDR3. (default: no)\n".
	"\t--reportPartial: include partial CDR3 in report. (default: no partial)\n".
	"\t--chainsInBarcode INT: number of chains in a barcode. (default: 2)\n"
	#"\t--secondary: output secondary chains. (default: no)\n"
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

sub GetCellType
{
	foreach my $g (@_)
	{
		if ( $g =~ /^IG/ )
		{
			return 0 ;
		}
		elsif ( $g =~ /^TR/ )
		{
			return 1 ;
		}
	}
	return -1 ;
}

sub GetDetailChainTypeFromGeneName
{
	my $g = $_[0] ;
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
	
	return -1 ;
}

sub GetDetailChainType
{
	my $i ;
	for ( $i = 1 ; $i <= 2 ; ++$i )
	{
		my $type = GetDetailChainTypeFromGeneName( $_[$i] ) ;
		return $type if ( $type != -1 ) ;
	}
	# In the case of only V gene is available.
	my $type = GetDetailChainTypeFromGeneName( $_[0] ) ;
	return $type ;
}

sub GetChainCellType
{
	if ( $_[0] <= 2 )
	{
		return 0 ; # B
	}
	elsif ( $_[0] <= 4 )
	{
		return 1 ; # abT
	}
	elsif ( $_[0] <= 6 ) 
	{
		return 2 ; # gdT
	}
	return -1 ;
}

# Use input V, J, C gene to report back the C gene if it is missing.
sub InferConstantGene
{
	my $ret = $_[2] ;
	my $i ;
	
	if ($_[2] ne "*")
	{
		$ret = (split /\*/, $ret)[0] ; # Remove the allele id from c gene
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
	
	# For TRA and TRD gene, we don't infer its constant gene.
	if ($_[0] =~ /^TR[AD]/ || $_[1] eq "*")
	{
		return $ret ;
	}

	for ( $i = 1 ; $i >= 0 ; --$i )
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

sub GetAaType 
{
	# 0- partial, 1- non-produtive, 2 - productive. The higher number the better
	my $a = $_[0] ;
	return 0 if ( $a eq "partial" ) ;
	return 1 if ( $a eq "out_of_frame" || $a =~ /_/ ) ;
	return 2 ;
}

# Test whether a is better than b.
sub BetterAA
{
	return GetAaType( $_[0] ) - GetAaType($_[1] ) ;
}

my $i ;
my $reportPartial = 0 ;
my $annotFile = "" ;
my $impute = 1 ;
my $imputeBCR = 0 ;
my $chainsInBarcode = 2 ;
for ( $i = 1 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "--reportPartial" )
	{
		$reportPartial = 1 ;
	}
	elsif ( $ARGV[$i] eq "--noImputation" )
	{
		$impute = 0 ;
	}
	elsif ( $ARGV[$i] eq "--imputeBCR" )
	{
		$imputeBCR = 1 ;
	}
	elsif ( $ARGV[$i] eq "-a" )
	{
		$annotFile = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--chainsInBarcode" )
	{
		$chainsInBarcode = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown option ", $ARGV[$i], "\n" ;
	}
}

# Store whether there is good assemblies that haven't got CDR3 from the annotation file.
my %barcodeChainInAnnot ; 
if ( $annotFile ne "" )
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

		if ( $cols[3] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):([0-9\.]+)/ )
		{
			#print($cols[3], "\t", $1, "\t", $6, "\n") ;
			@vCoord = ($1, $2, $3, $4, $5, $6) ;
		}
		else
		{
			#die "Wrong format $header\n" ;
			@vCoord = (-1, -1, -1, -1, -1, 0)
		}
		
		if ( $cols[5] =~ /\(([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):\(([0-9]+?)-([0-9]+?)\):([0-9\.]+)/ )
		{
			@jCoord = ($1, $2, $3, $4, $5, $6) ;
		}
		else
		{
			#die "Wrong format $header\n" ;
			@jCoord = (-1, -1, -1, -1, -1, 0)
		}

		my $cdr3Score = 0 ;
		if ( $cols[9] =~ /:(.+?)=/ )
		{
			$cdr3Score = $1 ;		
		}
		
		my $chainType = -1 ;
		
		if ( ( $vCoord[2] - $vCoord[1] >= 50 && $vCoord[5] >= 0.95 )
			|| ($cdr3Score > 0 && $vCoord[0] != -1 ) )
		{
			$chainType = GetDetailChainTypeFromGeneName( substr($cols[3], 0, 3) ) ; 		
		}
		elsif ($jCoord[2] - $jCoord[1] >= $jCoord[0] * 0.66 && $jCoord[5] >= 0.95
			|| ($cdr3Score > 0 && $jCoord[0] != -1 ) )
		{
			$chainType = GetDetailChainTypeFromGeneName( substr($cols[5], 0, 3) ) ; 		
		}

		if ( $chainType != -1 )
		{
			my @cols2 = split/_/, substr($cols[0], 1) ;
			my $barcode = join( "_", @cols2[0..scalar(@cols2)-2] ) ;
			my $key = $barcode."_".$chainType ;
			if ( !defined $barcodeChainInAnnot{$key} )
			{
				$barcodeChainInAnnot{$key} = 0 ;
			}
			$barcodeChainInAnnot{ $key } += $cols[2] ;
		}
	}
	close FP1 ;
}

# collect the read count for each chain from assembly id.
open FP1, $ARGV[0] ;
my %barcodeChainAbund ; 
my %barcodeChainRepresentAbund ;
my %barcodeChainRepresent ;
my %barcodeChainOther ;
my %barcodeShownup ;
my %barcodeChainAa ;
my @barcodeList ;
my %barcodeChainPartial ; # only store the information for partial cdr3s for later imputation

# Read in the report, store the information for each barcode.
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	#next if ( $reportPartial == 0 && $cols[9] == 0 ) ;

	my $assemblyId = $cols[0] ;
	my $vgene = (split /,/, $cols[2])[0] ;
	my $dgene = (split /,/, $cols[3])[0] ;
	my $jgene = (split /,/, $cols[4])[0] ;
	my $cgene = (split /,/, $cols[5])[0] ;
	#$cgene = InferConstantGene( $vgene, $jgene, $cgene ) ;

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
	
	my $info = join( ",", ($vgene, $dgene, $jgene, $cgene, $cols[8], $aa, $cols[10], $cols[0], $cols[11], $cols[12]) ) ;

	if ( $aa eq "partial" )
	{
		my $type = GetDetailChainType( $vgene, $jgene, $cgene ) ;
		if ($type == 0 || $type == 4 || $type ==6)
		{
			$type = 0 ;
		}
		else
		{
			$type = 1 ;
		}
		my $k = $barcode."_".$type ;
		@{$barcodeChainPartial{$k}} = () if (!defined $barcodeChainPartial{$k}) ;
		push @{$barcodeChainPartial{$k}}, $info ;
		next if ($reportPartial == 0) ;
	}

	if ( defined $barcodeChainAbund{$key} )
	{
		$barcodeChainAbund{$key} += $cols[10] ;
	}
	else
	{
		$barcodeChainAbund{$key} = $cols[10] ;
	}
	# out_of_frame or partial can only be on others part.
	if ( GetAaType($aa) < 2) 
	{
		@{$barcodeChainOther{$key}} = () if (!defined $barcodeChainOther{$key}) ;
		push @{$barcodeChainOther{$key}}, $info ;
		next ;
	}

	if ( defined $barcodeChainRepresent{ $key } )
	{
		if ( BetterAA($aa, $barcodeChainAa{$key}) > 0 || 
			( $cols[10] > $barcodeChainRepresentAbund{ $key } && BetterAA($aa, $barcodeChainAa{$key}) == 0 ) )
		{
			@{$barcodeChainOther{$key}} = () if (!defined $barcodeChainOther{$key}) ;
			push @{$barcodeChainOther{$key}}, $barcodeChainRepresent{$key} ; # the original representative becomes secondary.
			
			$barcodeChainRepresentAbund{$key} = $cols[10] ;
			$barcodeChainAa{$key} = $aa ;
			$barcodeChainRepresent{ $key } = $info ; 
		}
	}
	else
	{
		$barcodeChainRepresentAbund{$key} = $cols[10] ;
		$barcodeChainAa{$key} = $aa ;
		$barcodeChainRepresent{ $key } = $info ; 
	}
	
	#if ($barcode eq "GTACTTTGTACCAGTT-1")
	#{
	#	print($barcode, " ", $key, " ", $barcodeChainAbund{$key}, " ", $barcodeChainAa{$key}, "\n")
	#}
}
close FP1 ;

# Update barcode chain abundance with annotation. This should not affect the representative.
if ( $annotFile ne "" )
{
	for my $key (keys %barcodeChainAbund)
	{
		if ( defined $barcodeChainInAnnot{$key} )
		{
			$barcodeChainAbund{$key} =  $barcodeChainInAnnot{$key} ;
		}
	}
}

# Generate the output content.
my %barcodeOutput ;
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
	my $secondaryChain1 = "*" ;
	my $secondaryChain2 = "*" ;
	for ( $i = 0 ; $i < 7 ; ++$i )
	{
		my $key = $barcode."_".$i ;
		# gdT should have more stringent criterion.
		last if ($i >= 5 && $maxTag != -1 ) ;

		if ( defined $barcodeChainAbund{ $key } && $barcodeChainAbund{ $key } > $max )
		{
			$max = $barcodeChainAbund{ $key } ;
			$maxTag = $i ;
		}
	}
	
	if ( $maxTag >= 5 && $annotFile ne "" )
	{
		# use annotation file to further filter gdT
		my $tag = 0 ;
		for ( $i = 0 ; $i < 5 ; ++$i )
		{
			$tag |= (1<<$i) if ( defined $barcodeChainInAnnot{$barcode."_".$i} ) ;
		}
		next if ( ($tag&3) == 3 || ($tag&5) == 5 || ($tag&24) == 24 ) ;
	}

	my @otherList = (0, 1, 2, 3, 4, 5, 6) ;
	my @representativeUsed = (0, 0, 0, 0, 0, 0, 0, 0) ;
	if ( $maxTag <= 2 )
	{
		$mainType = 0 ;
		my $keyH = $barcode."_0" ;
		my $keyK = $barcode."_1" ;
		my $keyL = $barcode."_2" ;
		$chain1 = $barcodeChainRepresent{ $keyH } if ( defined $barcodeChainRepresent{ $keyH } ) ;
		$representativeUsed[0] = 1 ;

		if ( defined $barcodeChainRepresent{ $keyK } && defined $barcodeChainRepresent{ $keyL } )
		{
			if ( $barcodeChainAbund{ $keyK } >= $barcodeChainAbund{ $keyL } )
			{
				$chain2 = $barcodeChainRepresent{ $keyK } ;
				$representativeUsed[1] = 1 ;
			}
			else
			{
				$chain2 = $barcodeChainRepresent{ $keyL } ;
				$representativeUsed[2] = 1 ;
			}
		}
		elsif ( defined $barcodeChainRepresent{ $keyK } )
		{
			$chain2 = $barcodeChainRepresent{ $keyK } ;
			$representativeUsed[1] = 1 ;
		}
		elsif ( defined $barcodeChainRepresent{ $keyL } )
		{
			$chain2 = $barcodeChainRepresent{ $keyL } ;
			$representativeUsed[2] = 1 ;
		}
		
		$cellType = "B" ;
	}
	else
	{
		my $key1 ;
		my $key2 ;
		if ( $maxTag <= 4 )
		{
			$key1 = $barcode."_4" ;
			$key2 = $barcode."_3" ;
			$cellType = "abT" ;
			@otherList = (3, 4, 5, 6, 0, 1, 2) ;
			$representativeUsed[3] = 1 ;
			$representativeUsed[4] = 1 ;
		}
		elsif ( $maxTag <= 6 )
		{
			$key1 = $barcode."_6" ;
			$key2 = $barcode."_5" ;
			$cellType = "gdT" ;
			@otherList = (5, 6, 3, 4, 0, 1, 2) ;
			$representativeUsed[5] = 1 ;
			$representativeUsed[6] = 1 ;
		}
		$chain1 = $barcodeChainRepresent{ $key1 } if ( defined $barcodeChainRepresent{ $key1 } ) ;
		$chain2 = $barcodeChainRepresent{ $key2 } if ( defined $barcodeChainRepresent{ $key2 } ) ;
	}

	# Add other chains to the secondary chain field.
	foreach $i (@otherList)
	{
		my $key = $barcode."_".$i ;
		if ( defined $barcodeChainRepresent{ $key } ) 
		{
			my $addition = "" ;
			$addition = $barcodeChainRepresent{$key} if ($representativeUsed[$i] == 0);
			if ($addition ne "")
			{
				$addition = $addition.";".join(";", @{$barcodeChainOther{$key}}) if (defined $barcodeChainOther{$key}) ;
			}
			else
			{
				$addition = join(";", @{$barcodeChainOther{$key}}) if (defined $barcodeChainOther{$key}) ;
			}

			next if ($addition eq "") ;

			if ($i == 0 || $i == 4 || $i == 6)
			{
				if ($secondaryChain1 eq "*")
				{
					$secondaryChain1 = $addition ;
				}
				else
				{
					$secondaryChain1 = $secondaryChain1.";".$addition ;
				}
			}
			else
			{
				if ($secondaryChain2 eq "*")
				{
					$secondaryChain2 = $addition ;
				}
				else
				{
					$secondaryChain2 = $secondaryChain2.";".$addition ;
				}
			}
		}
	}
	next if ( $chain1 eq "*" && $chain2 eq "*" ) ;
	#print( join( "\t", ($barcode, $cellType, $chain1, $chain2, $secondaryChain1, $secondaryChain2 ) ), "\n" ) ;
	if ($chainsInBarcode == 1)
	{
		if ($chain1 eq "*" && $chain2 ne "*")
		{
			$chain1 = $chain2 ;
			$chain2 = "*" ;
			$secondaryChain1 = $secondaryChain2 ;
			$secondaryChain2 = "*" ;
		}
		elsif ($chain1 ne "*" && $chain2 ne "*")
		{
			my $abund1 = (split /,/, $chain1)[6] ;
			my $abund2 = (split /,/, $chain2)[6] ;
			if ($abund2 > $abund1)
			{
				$chain1 = $chain2 ;
				$chain2 = "*" ;
				$secondaryChain1 = $secondaryChain2 ;
				$secondaryChain2 = "*" ;
			}
		}
	}
	@{$barcodeOutput{$barcode}} = ($cellType, $chain1, $chain2, $secondaryChain1, $secondaryChain2) ;
}

sub IsACompatibleToB
{
	my @colsA = split(/,/, $_[0]) ;
	my @colsB = split(/,/, $_[1]) ;
	my $partial = $_[2] ;
	if (GetCellType($colsA[0], $colsA[2], $colsA[3]) != GetCellType($colsB[0], $colsB[2], $colsB[3]))
	{
		return 0 ;
	}
	
	for my $i (0, 2, 3)
	{
		if ($colsA[$i] ne "*" && $colsB[$i] ne "*" && $colsA[$i] ne $colsB[$i])	
		{
			return 0 ;
		}
	}
	my $pattern = $colsA[4] ;
	if ( $partial == 1 && ( $colsB[4] =~ /^$pattern/ || $colsB[4] =~ /$pattern$/ ) )
	{
		return 1 ;
	}
	elsif ( $partial == 0 && ($colsB[4] eq $pattern))
	{
		return 1 ;
	}
	else
	{
		return 0 ;
	}
}

# Go through the list again for imputation.
if ($impute == 1)
{
	my %cdr3ToBarcodes ;
	for my $barcode (keys %barcodeOutput)
	{
		my @cols = @{$barcodeOutput{$barcode}} ;
		next if ($cols[1] eq "*" && $cols[2] eq "*") ;
		my $i ;
		for ($i = 0 ; $i <= 1 ; ++$i)
		{
			if ($cols[$i + 1] ne "*")
			{
				my $cdr3 = (split /,/, $cols[$i + 1])[4] ;
				my $key = $cdr3."_".$i ;
				@{$cdr3ToBarcodes{$key}} = () if (!defined $cdr3ToBarcodes{$key} ) ;
				push @{$cdr3ToBarcodes{$key}}, $barcode ;
			}
		}
	}

	foreach my $barcode (keys %barcodeOutput)
	{
		my @cols = @{$barcodeOutput{$barcode}} ;
		next if ($cols[1] eq "*" && $cols[2] eq "*") ;
		next if ($cols[1] ne "*" && $cols[2] ne "*") ;
		next if ($cols[0] eq "B" && $imputeBCR == 0) ; # only impute for tcr

		my $missingChain = 0 ;
		$missingChain = 1 if ($cols[2] eq "*") ;
		# For imputation, it must have partial CDR3 on the missing chain
		# 	and the chain with CDR3 must show up in some other barcodes.
		next if ( !defined $barcodeChainPartial{ $barcode."_".$missingChain} ) ;
		my $cdr3 = (split(/,/, $cols[2-$missingChain]))[4] ;
		my $candidateOtherBarcode = "" ;
		my $multipleOthers = 0 ;
		foreach my $otherBarcode (@{$cdr3ToBarcodes{$cdr3."_".(1-$missingChain)}})
		{
			my @otherCols = @{$barcodeOutput{$otherBarcode}} ;
			next if ($otherCols[$missingChain + 1] eq "*") ;
			next if (!IsACompatibleToB($cols[2 - $missingChain], $otherCols[2 - $missingChain], 0)) ;
			foreach my $partialInfo (@{$barcodeChainPartial{$barcode."_".$missingChain}})
			{
				if (IsACompatibleToB($partialInfo, $otherCols[$missingChain + 1], 1))
				{
					if ($candidateOtherBarcode ne "")
					{
						if (IsACompatibleToB(${$barcodeOutput{$candidateOtherBarcode}}[$missingChain + 1],
							$otherCols[$missingChain + 1], 0))
						{
							$multipleOthers = 1 ;
						}
					}
					$candidateOtherBarcode = $otherBarcode ;
					last ;
				}
			}
			last if ($multipleOthers == 1) ;
		}
		next if ($candidateOtherBarcode eq "") ;

		# Put in the imputed results.
		my $s = ${$barcodeOutput{$candidateOtherBarcode}}[$missingChain + 1] ;
		# Change the orginated assembly id to "impute".
		@cols = split /,/, $s ;
		$cols[7] = "impute_from_$candidateOtherBarcode" ;
		$s = join(",", @cols) ;
		${$barcodeOutput{$barcode}}[$missingChain + 1] = $s ;
	}
}

# Output the results
print( "#barcode\tcell_type\tchain1\tchain2\tsecondary_chain1\tsecondary_chain2\n" ) ;
foreach my $barcode (keys %barcodeOutput)
{
	print( $barcode, "\t", join( "\t", @{$barcodeOutput{$barcode}}) , "\n" ) ;
}
