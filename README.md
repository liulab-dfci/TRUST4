TRUST4
=======

Described in: 

TRUST4

	Copyright (C) 2018- and GNU GPL by Li Song, Shirley Liu

Includes portions copyright from: 

	samtools - Copyright (C) 2008-, Genome Research Ltd, Heng Li
	

### What is TRUST4?

Tcr Receptor Utilities for Solid Tissue (TRUST) is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from solid tissues, including tumors. TRUST4 performs de novo assembly on V, J, C genes including the hypervariable complementarity-determining region 3 (CDR3) and reports consensus of BCR/TCR sequences. TRUST4 then realigns the contigs to IMGT reference gene sequences to report the corresponding information. TRUST4 supports both single-end and paired-end sequencing data with any read length. 

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/TRUST4), e.g. with `git clone https://github.com/mourisl/TRUST4.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run TRUST4 without specifying the directory, you can either add the directory of TRUST4 to the environment variable PATH or create a soft link ("ln -s") of the file "run-trust4" to a directory in PATH.

TRUST4 depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib).


### Usage

	Usage: ./run-trust4 [OPTIONS]
		Required:
			-b STRING: path to bam file
			-f STRING: path to the fasta file coordinate and sequence of V/D/J/C genes\n".
		Optional:
			--ref STRING: path to detailed V/D/J/C gene reference file, such as from IMGT database. (default: not used). (recommended) 
			-o STRING: prefix of output files. (default: inferred from file prefix)
			-t INT: number of threads (default: 1)\n".
			--stage INT: start TRUST4 on specified stage (default: 0):\n".
				0: start from beginning (candidate read extraction)\n".
				1: start from assembly\n".
				2: start from annotation\n". 

### Input/Output

The primary input to TURST4 is the alignment of RNA-seq reads in BAM format(-b), the file containing the genomic sequence and coordinate of V,J,C genes(-f), and the reference database sequence containing annotation information, such as IMGT (--ref).

TRUST4 outputs several files. trust_raw.out, trust_final.out are the contigs and corresponding nucleotide weight. trust_annot.out is in fasta format for the annotation of the consensus assembly. And trust_cdr3.out reports the CDR1,2,3 and gene information for each consensus assemblies.

Each header of trust_annot.out is split into fields:

	consensus_id consensus_length average_coverage annotations

"annotaions" also has several field, corresponding to annotation of V,J,C, CDR1, CDR2 and CDR3 respectively. For the annotation of the genes, it follows the pattern 

	gene_name(reference_gene_length):(consensus_start-consensus_end):(reference_start-reference_length):similarity
	
Each type of genes has at most three gene candidate ranked by their similarity. For the annotation of CDRs, it follows the pattern:

	CDRx(consensus_start-consensus_end):score=sequence
	
For CDR1,2, score is similarity. for CDR3, score 0.00 means partial CDR3, score 1.00 means CDR3 with imputed nucleotides and other numbers means the motif signal strength with 100.00 as strongest.

The coordinate is 0-based.

The output trust_cdr3.out is a tsv file. The fields are:

	consensus_id	index_within_consensus	V_gene	J_gene	C_gene	CDR1	CDR2	CDR3	CDR3_score	read_fragment_count


### Practical notes

* Build custom V,J,C gene database (files for -f and --ref)

To generate the file specified by "-f", you need the reference genome of the species you are interested in and corresponding genome annotation GTF file. Then you can use command 
	
	perl BuildDatabaseFa.pl reference.fa annotation.gtf bcr_tcr_gene_name.txt > bcrtcr.fa

to generate the input for "-f". The "bcr_tcr_gene_name.txt" is provided as "human_vdjc.list" in the repository.

Normally, the file specified by "--ref" is downloaded from IMGT website and then supplemented by sequence of constant genes, For example, for human, you can use command

	perl BuildImgtAnnot.pl bcrtcr.fa Cgene.list Homo_sapien > IMGT+C.fa

to generate the input for "--ref". The bcrtcr.fa is the file generated in previous step (for -f). Cgene.list is provided in the repository. The available species name can be found on [IMGT FTP](http://www.imgt.org//download/V-QUEST/IMGT_V-QUEST_reference_directory/).

* Simple report

The default report of TRUST4 is trust_cdr3.out is focus on DNA level. We provide a script "trust-simplerep.pl" that can convert the output to amino acids, collapse repeated entries and normalize the frequency. The command is

	perl trust-simplerep.pl trust_cdr3.out > trust_report.out

The output is a tsv file and columns in the simple report are 

	read_count frequency(proportion of read_count) CDR3_dna CDR3_amino_acids V D J C 
	
The D gene column is a place holder for now and is always "\*" indicating missing. Note that for frequency, the BCR(IG) and TCR(TR) chains are normalized respectively. If you are interested in a subset of chains, you can "grep" those from trust_cdr3.out and run trust-simplerep.pl on the subset.
 
### Example

The directory './example' in this distribution contains one BAM files as input for TRUST4. Run TRUST4 with:

	./run-trust4 -b example/example.bam -f bcrtcr.fa --ref IMGT+C.fa

The run will generate the files TRUST_example_raw.out, TRUST_example_final.out, TRUST_example_annot.fa and TRUST_example_cdr3.out and several fq/fa files.

### Terms of use

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General
Public License along with this program; if not, you can obtain one from
http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
### Support

Create a [GitHub issue](https://github.com/mourisl/TRUST4/issues).
