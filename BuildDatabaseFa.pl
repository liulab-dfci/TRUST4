#/usr/bin/env perl

use warnings ;
use strict ;

die "usage: perl BuildDatabaseFa.pl reference.fa annotation.gtf interested_gene_name_list > output.fa\n" if ( @ARGV != 3 ) ;

my %genome ;
my %interestedGeneName ;


# Read in the reference genome.
open FP1, $ARGV[0] ;
my $chrom = "" ;
my $seq = "" ;
my $hasChrPrefix = 0 ;
while ( <FP1> )
{
	if ( /^>/ )
	{
		$genome{ $chrom } = $seq if ( $chrom ne "" ) ;
		$seq = "" ;
		$chrom = substr( ( split )[0], 1 ) ;
		if ( $chrom =~ /^c/ )
		{
			$hasChrPrefix = 1 ;
		}
	}
	else
	{
		chomp ;
		$seq .= $_ ;
	}
}
$genome{ $chrom } = $seq if ( $chrom ne "" ) ;
close FP1 ;

# Read in the gene name we interested in.
open FP1, $ARGV[2] ;
while ( <FP1> )
{
	chomp ;
	$interestedGeneName{ $_ } = 1 ;
}
close FP1 ;

# Read in the gtf file and output interested sequences.
open FP1, $ARGV[1] ;
$seq = 0 ;
#chr14	HAVANA	UTR	21712321	21712330	.	+	.	gene_id "ENSG00000211776.2"; transcript_id "ENST00000390424.2"; gene_type "TR_V_gene"; gene_name "TRAV2"; transcript_type "TR_V_gene"; transcript_name "TRAV2-201"; exon_number 1; exon_id "ENSE00001508005.2"; level 2; protein_id "ENSP00000438195.1"; transcript_support_level "NA"; tag "mRNA_end_NF"; tag "cds_end_NF"; tag "basic"; tag "appris_principal_1"; havana_gene "OTTHUMG00000168980.2"; havana_transcript "OTTHUMT00000401875.2";
my $prevTname = "" ;
my $gname = "" ;
my $strand = "." ;
my @range ;

sub OutputGene
{
	my $i ;

	$chrom = $range[0] ;
	my $start = $range[1] ;
	my $end = $range[-1] ;

	my $output = "" ;
	
	if ($strand eq "-" ) 
	{
		$start = $range[-2] ;
		$end = $range[2] ;
	}
	die "Unknown chrom id $chrom " if ( !defined $genome{ $chrom } ) ;
	
	print ">$gname $chrom $start $end $strand\n" ; 
	for ( $i = 0 ; $i < scalar( @range ) ; $i += 3 )
	{
		my $tmp = uc( substr( $genome{ $range[$i] }, $range[$i + 1] - 1, $range[$i + 2] - $range[$i + 1] + 1 ) ) ;
		if ( $strand eq "-" )
		{
			$tmp = reverse( $tmp ) ;
			$tmp =~ tr/ACGT/TGCA/ ;
		}
		$output .= $tmp ;
	}
	print( "$output\n" ) ;
}

while ( <FP1> )
{
	next if ( /^#/ ) ;
	chomp ;
	my @cols = split /\t/ ;
	next if ( $cols[2] ne "exon" ) ;

	my $tname ;	
	if ( $cols[8] =~ /transcript_name \"(.*?)\"/ )
	{
		#print $1, "\n" ; 
		$tname = $1 ;
	}
	else
	{
		die "No transcript_name", $_, "\n" ;
	}


	if ( $tname ne $prevTname )
	{
		if ( (defined $interestedGeneName{ $gname } ) && @range > 0 )
		{
			OutputGene() ;
		}

		$prevTname = $tname ;
		if ( $cols[8] =~ /gene_name \"(.*?)\"/ )
		{
			#print $1, "\n" ; 
			$gname = uc($1) ;
		}
		else
		{
			die "No gene_name: ", $_, "\n" ;
		}
		$strand = $cols[6] ;
		undef @range ;
	}
	
	if ( $hasChrPrefix == 1 && !( $cols[0] =~ /^c/) )
	{
		$cols[0] = "chr".$cols[0] ;
	}
	elsif ( $hasChrPrefix == 0 && $cols[0] =~ /^c/ )
	{
		$cols[0] = substr( $cols[0], 3 ) ;
	}
	push @range, $cols[0], $cols[3], $cols[4] ;			
}

if ( (defined $interestedGeneName{ $gname } ) && @range > 0 )
{
	OutputGene() ;
}


close FP1 ;
