PsiCLASS
=======

Described in: 

TRUST4

	Copyright (C) 2018- and GNU GPL by Li Song, Shirley Liu

Includes portions copyright from: 

	samtools - Copyright (C) 2008-, Genome Research Ltd, Heng Li
	

### What is TRUST4?

Tcr Receptor Utilities for Solid Tissue (TRUST) is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from solid tissues, including tumors. TRUST4 performs de novo assembly on V, J, C genes including the hypervariable complementarity-determining region 3 (CDR3) and reports contigs containing BCR/TCR sequences. TRUST4 then realigns the contigs to IMGT reference gene sequences to report the corresponding information. TRUST4 supports both single-end and paired-end sequencing data with any read length. 

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
    			(recommended) --ref STRING: path to detailed V/D/J/C gene reference file, such as from IMGT database. (default: not used).
    			-o STRING: prefix of output files. (default: inferred from file prefix)
 		   	-t INT: number of threads (default: 1)\n".
    			--stage INT: start TRUST4 on specified stage (default: 0):\n".
    				0: start from beginning (candidate read extraction)\n".
    				1: start from assembly\n".
    				2: start from annotation\n". 

### Input/Output

The primary input to TURST4 

### Practical notes

*Build custom V,J,C gene database (files for -f and --ref)* To generate The file specified by "-f"  


*Concise report* The default report of TRUST4 is trust_cdr3.out,  


### Example

The directory './example' in this distribution contains two BAM files, along with an example of a BAM list file. Run PsiCLASS with:

	./run-trust4 -b example/s1.bam,example/s2.bam

The run will generate the files 

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
