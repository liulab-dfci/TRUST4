TRUST4
=======

Described in: 

Song, L., Cohen, D., Ouyang, Z. et al. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nat Methods (2021). https://doi.org/10.1038/s41592-021-01142-2

	Copyright (C) 2018-, Li Song, X. Shirley Liu

Includes portions copyright from: 

	samtools - Copyright (C) 2008-, Genome Research Ltd, Heng Li
	

### What is TRUST4?

TRUST4 is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from fluid and solid tissues, including tumors. TRUST4 performs de novo assembly on V, J, C genes including the hypervariable complementarity-determining region 3 (CDR3) and reports consensus contigs of BCR/TCR sequences. TRUST4 then realigns the contigs to IMGT reference gene sequences to identify the corresponding gene and CDR3 details. TRUST4 supports both single-end and paired-end bulk or single-cell sequencing data with any read length. 

### Install

1. Clone the [GitHub repo](https://github.com/liulab-dfci/TRUST4), e.g. with `git clone https://github.com/liulab-dfci/TRUST4.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run TRUST4 without specifying the directory, you can either add the directory of TRUST4 to the environment variable PATH or create a soft link ("ln -s") of the file "run-trust4" to a directory in PATH.

TRUST4 depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib). For MacOS, TRUST4 has been successfully compiled with gcc_darwin17.7.0 and gcc_9.2.0 installed by Homebrew.

TRUST4 is also available form [Bioconda](https://anaconda.org/bioconda/trust4). You can install TRUST4 with `conda install -c bioconda trust4` or use the docker container `docker pull quay.io/biocontainers/trust4:<tag>` (see [trust4/tags](https://quay.io/repository/biocontainers/trust4?tab=tags) for valid values for <tag>). 

### Usage

	Usage: ./run-trust4 [OPTIONS]
		Required:
			-b STRING: path to bam file
			-1 STRING -2 STRING: path to paired-end read files
			-u STRING: path to single-end read file
			-f STRING: path to the fasta file coordinate and sequence of V/D/J/C genes
		Optional:
			--ref STRING: path to detailed V/D/J/C gene reference file, such as from IMGT database. (default: not used). (recommended) 
			-o STRING: prefix of output files. (default: inferred from file prefix)
			--od STRING: the directory for output files. (default: ./)
			-t INT: number of threads (default: 1)
			-k INT: the starting k-mer size for indexing contigs (default: 9)
			--barcode STRING: if -b, bam field for barcode; if -1 -2/-u, file containing barcodes (defaul: not used)
			--barcodeLevel STRING: barcode is for cell or molecule (default: cell)
			--barcodeWhitelist STRING: path to the barcode whitelist (default: not used)
			--barcodeTranslate STRING: path to the barcode translate file (default: not used)
			--UMI STRING: if -b, bam field for 10x Genomics-like UMI; if -1 -2/-u, file containing 10x Genomics-like UMIs (default: not used)
			--readFormat STRING: format for read, barcode and UMI files (example: r1:0:-1,r2:0:-1,bc:0:15,um:16:-1 for paired-end files with barcode and UMI)
			--repseq: the data is from bulk,non-UMI-based TCR-seq or BCR-seq (default: not set)
			--contigMinCov INT: ignore contigs that have bases covered by fewer than INT reads (default: 0)
			--minHitLen INT: the minimal hit length for a valid overlap (default: auto)
			--mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)
			--skipMateExtension: do not extend assemblies with mate information, useful for SMART-seq (default: not used)
			--abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)
			--assembleWithRef: conduct the assembly with --ref file (default: use -f file)\n".
			--noExtraction: directly use the files from provided -1 -2/-u to assemble (default: extraction first)
			--outputReadAssignment: output read assignment results to the prefix_assign.out file (default: no output)
			--stage INT: start TRUST4 on specified stage (default: 0)
				0: start from beginning (candidate read extraction)
				1: start from assembly
				2: start from annotation
				3: start from generating the report table
			--clean INT: clean up files. 0: no clean. 1: clean intermediate files. 2: only keep AIRR files. (default: 0)

### Input/Output

The primary input to TRUST4 is the alignment of RNA-seq reads in BAM format(-b), the file containing the genomic sequence and coordinate of V,J,C genes(-f), and the reference database sequence containing annotation information, such as IMGT (--ref).

An alternative input to TRUST4 is the raw RNA-seq files in fasta/fastq format (-1/-2 for paired; -u for single-end). You still need the files like -f, --ref from above. In this case, you can directly use IMGT's seuqence file for -f. 

TRUST4 outputs several files. trust_raw.out, trust_final.out are the contigs and corresponding nucleotide weight. trust_annot.fa is in fasta format for the annotation of the consensus assembly. trust_cdr3.out reports the CDR1,2,3 and gene information for each consensus assemblies. And trust_report.tsv is a report file focusing on CDR3 and is compatible with other repertoire analysis tool such as VDJTools. 

Each header of trust_annot.fa is split into fields:

	consensus_id consensus_length average_coverage annotations

"annotations" also has several field, corresponding to annotation of V,D,J,C, CDR1, CDR2 and CDR3 respectively. For the annotation of the genes, it follows the pattern 

	gene_name(reference_gene_length):(consensus_start-consensus_end):(reference_start-reference_length):similarity
	
Each type of genes has at most three gene candidate ranked by their similarity. For the annotation of CDRs, it follows the pattern:

	CDRx(consensus_start-consensus_end):score=sequence
	
For CDR1,2, score is similarity. for CDR3, score 0.00 means partial CDR3, score 1.00 means CDR3 with imputed nucleotides and other numbers means the motif signal strength with 100.00 as strongest. The coordinate is 0-based.

The output trust_cdr3.out is a tsv file. The fields are:

	consensus_id	index_within_consensus	V_gene	D_gene	J_gene	C_gene	CDR1	CDR2	CDR3	CDR3_score	read_fragment_count CDR3_germline_similarity complete_vdj_assembly
	
Please note that CDR3_score in trust_cdr3.out has been divided by 100, so 1.00 is the maximum score and 0.01 means imputed CDR3.

The output trust_report.tsv is a tsv file. The fileds are:
	
	read_count	frequency(proportion of read_count)	CDR3_dna	CDR3_amino_acids	V	D	J	C	consensus_id consensus_id_complete_vdj

For frequency, the BCR(IG) and TCR(TR) chains are normalized respectively. In the amino acid sequence, "_" represents stop codon, and "?" represents ambiguous nucleotide "N" in codon.

The output trust_airr.tsv follows [the AIRR format](https://docs.airr-community.org/en/latest/datarep/rearrangements.html). 

### Practical notes

* #### Build custom V,J,C gene database (files for -f and --ref)

Normally, the file specified by "--ref" is downloaded from IMGT website, For example, for human, you can use command

	perl BuildImgtAnnot.pl Homo_sapien > IMGT+C.fa

The available species name can be found on [IMGT FTP](http://www.imgt.org//download/V-QUEST/IMGT_V-QUEST_reference_directory/).

If your input data for TRUST4 is raw FASTQ files, you can use the IMGT file for the "-f" option. If your input data for TRUST4 is BAM files, you need to generate another file for "-f". To do that, you need the reference genome (e.g. hg38 for human, or mm10 for mouse) of the species you are interested in and corresponding genome annotation GTF file (e.g. gencode v35 for human, or gencode mV21 for mouse). Then you can use command 
	
	perl BuildDatabaseFa.pl reference.fa annotation.gtf bcr_tcr_gene_name.txt > bcrtcr.fa

to generate the input for "-f". The "bcr_tcr_gene_name.txt" is provided as "human_vdjc.list" in the repository.

The IMGT+C.fa can also be used to generate "bcr_tcr_gene_name.txt" file with command:

	grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq > bcr_tcr_gene_name.txt  

* #### 10X Genomics data and barcode-based single-cell data

When given barcode, TRUST4 only assembles the reads with the same barcode together. For 10X Genomics data, usually the input is the BAM file from cell-ranger, and you can use "--barcode" to specify the field in the BAM file to specify the barcode: e.g. "--barcode CB".

If your input is raw FASTQ files, you can use "--barcode" to specify the barcode file and use "--readFormat" to tell TRUST4 how to extract barcode information. The "--readFormat" option can also specify the extraction for read1, read2 and UMI. The value for this argument is a comma-separated string, each field in the string is also a semi-comma-splitted string

	[r1|r2|bc|um]:start:end:strand

The start and end are inclusive and -1 means the end of the read. You may use multiple fields to specify non-consecutive segments, e.g. bc:0:15,bc:32:-1. The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction. The strand symbol can be omitted if it is '+' and is ignored on r1 and r2. For example, when the barcode is in the first 16bp of read1, one can use the option `-1 read1.fq.gz -2 read2.fq.gz --barcode read1.fq.gz --read-format bc:0:15,r1:16:-1`.

TRUST4 supports using wildcard in the -1 -2/-u option, so a typical way to run 10X Genomics single-end data is by:

	run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -u path_to_10X_fastqs/*_R2_*.fastq.gz --barcode path_to_10X_fastqs/*_R1_*.fastq.gz --readFormat bc:0:15 --barcodeWhitelist cellranger_folder/cellranger-cs/VERSION/lib/python/cellranger/barcodes/737K-august-2016.txt [other options]

The exact options depend on your 10X Genomics ikit.

Besides, TRUST4 can translate input cell barcodes to another set of barcodes. You can specify the translation file through the option --barcodeTranslate. The translation file is a two-column tsv/csv file with the translated barcode on the first column and the original barcode on the second column. This option also supports combinatorial barcoding, such as SHARE-seq. TRUST4 can translate each barcode segment provided in the second column to the ID in the first column and add "-" to concatenate the IDs in the output.

In the output, the abundance in the report will use the number of barcodes for the CDR3 instead of read count. TRUST4 will also generate the file trust_barcode_report.tsv. In this file, TRUST4 will pick the most abundance pair of chains as the representative for the barcode(cell). The format is:

	barcode	cell_type	IGH/TRB/TRD_information	IGK/IGL/TRA/TRG_information	secondary_chain1_information	secondary_chain2_information

For the chain information it is in CSV format:
	
	V_gene,D_gene,J_gene,C_gene,cdr3_nt,cdr3_aa,read_cnt,consensus_id,CDR3_germline_similarity,consensus_complete_vdj

TRUST4 also converts the barcode report file to the trust_barcode_airr.tsv file to follow the AIRR format.

* #### SMART-Seq data

We provide a wrapper "trust-smartseq.pl" to process the files from platforms like SMART-seq. The user shall give the path to each file in a text file. An example command can be

	perl trust-smartseq.pl -1 read1_list.txt -2 read2_list.txt -t 8 -f hg38_bctcr.fa --ref human_IMGT+C.fa -o TRUST

The script will create two files: TRUST_report.tsv for general summary and TRUST_annot.fa for assemblies. The formats are described above. Each cell's name is inferred by the file name before the first ".".

* #### UMI

For 10x Genomics data, TRUST4 supports UMI-based abundance estimation. You can use --UMI to specify the UMI sequence file or the field in the BAM file. If the sequence contains non-UMI information, you can use --readFormat with keyword "um" to specify the UMI sequence range. 

Note that in 10x Genomics data, UMI plus the cell barcode is the real unique molecular identifier. In other platforms, the UMI can be real unique and be regarded as molecule barcode. You can run trust4 with "--barcode UMIfile --barcodeLevel molecule" to specify UMI as molecule barcode. In this output, the represented chain information is in the chain1 column of the trust_barcode_report.tsv file.

* #### Simple report

The last step of generating simple report can be done with the command:

	perl trust-simplerep.pl trust_cdr3.out > trust_report.out

If you are interested in a subset of chains, you can "grep" those from trust_cdr3.out and run trust-simplerep.pl on the subset.

* #### Annotation only

You can use the "annotator" from TRUST4 to annotate the V,D,J,C genes and CDRs for any given sequences, just like using IgBLAST or IMGT/VQuest. To obtain the annotation in AIRR format for human sequences with eight threads, you can use the command

	./annotator -f human_IMGT+C.fa -a input.fa --fasta -t 8 --needReveserComplement --noImpute --outputFormat 1 > annotation.tsv 
 
### Example

The directory './example' in this distribution contains one BAM files as input for TRUST4. Run TRUST4 with:

	./run-trust4 -b example/example.bam -f hg38_bcrtcr.fa --ref human_IMGT+C.fa

The run will generate the files TRUST_example_raw.out, TRUST_example_final.out, TRUST_example_annot.fa, TRUST_example_cdr3.out, TRUST_example_report.tsv and several fq/fa files in seconds. The results should be the same as the files in the example folder. You can check whether TRUST4 is properly installed by running the command "bash trust-example-test.sh" in the TRUST4 folder using this example data.

The directory also contains two fastq files, and you can run TRUST4 with:

	./run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -1 example/example_1.fq -2 example/example_2.fq -o TRUST_example

Note that the requirement of the hg38_bcrtcr.fa is that it contains genomic coordinates of the V, D, J and C genes which are crucial for extracting candidate reads in alignment BAM files. The coordinate information is not needed for fastq input, so you can use the IMGT reference file for the -f option. This is useful when analyzing species without reference genomes or genome annotations. The example command can be: 

	./run-trust4 -f human_IMGT+C.fa --ref human_IMGT+C.fa -1 example/example_1.fq -2 example/example_2.fq -o TRUST_example

The run will generate the same files as from BAM input 

### Miscellaneous

The evaluation instructions and scripts in TRUST4's manuscript is available at: https://github.com/liulab-dfci/TRUST4_manuscript_evaluation .

### Support

Create a [GitHub issue](https://github.com/liulab-dfci/TRUST4/issues). We will typically respond within a day or two, but it could take longer, e.g. a month, for fixing bugs and adding features.
