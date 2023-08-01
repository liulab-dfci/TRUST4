#!/bin/bash

./run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -b example/example.bam -o example_test

if [ ! -e "example_test_report.tsv" ]; then
	echo "Failed to run TRUST4. Please check whether TRUST4 is compiled with command 'make' or the other files are in the right paths."	
	exit
fi

if [[ $(diff <(sort example_test_report.tsv) <(sort example/TRUST_example_report.tsv) | wc -l) -ge 1 ]]; then
	echo "Results do not match the pregenerated example results. Please check TRUST4's version or the example files are in the folder."	
	exit
fi

rm example_test_*

echo "TRUST4 is ready to use."
