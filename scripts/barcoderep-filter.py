#!/usr/bin/env python3 

import sys
import argparse

parser = argparse.ArgumentParser(description="Filter TRUST4 barcode report file for noise like diffused mRNAs. Output to stdout")
parser.add_argument('-b', help="barcode_report file", required=True, dest="barcode_report")
parser.add_argument('-a', help="annotation file", dest="annot")
parser.add_argument('--highAbund', help="The minimum abundance to be regarded as potential source of diffusion", default=50.0, dest="highAbund")
parser.add_argument('--diffuseFrac', help="The maximum fraction of the diffusion source abundance to be regarded as noise", default=0.02, dest="diffuseFrac")

args = parser.parse_args()

barcodeReport = args.barcode_report
barcodeInfo = {}
highAbundCdr3Info = {}
highAbund = float(args.highAbund)
diffuseFrac = float(args.diffuseFrac)
assembly = {}

fp = open(barcodeReport)
for line in fp:
	if (line[0] == '#'):
		continue
	cols = line.rstrip().split()
	chain1Cols = cols[2].split(',')
	chain2Cols = cols[3].split(',')
#AGAGTGGTCTATCCTA-1	abT	TRBV6-6*01,TRBD2*01,TRBJ2-5*01,*,TGTGCCAGTCTACTTGGGGGGACCCAGTACTTC,CASLLGGTQYF,3.61,AGAGTGGTCTATCCTA-1_48281,100.00,0	TRAV10*01,*,TRAJ34*01,*,TGTGTGGTGAGCGCCCGCACCGACAAGCTCATCTTT,CVVSARTDKLIF,1.00,AGAGTGGTCTATCCTA-1_65563,100.00,0	*	*
	if (len(chain1Cols) > 1 and float(chain1Cols[6]) >= highAbund):
		cdr3 = chain1Cols[4]
		abund = float(chain1Cols[6])
		if (cdr3 not in highAbundCdr3Info):
			highAbundCdr3Info[cdr3] = {}
		highAbundCdr3Info[cdr3][cols[0]] = [0, abund] 
	if (len(chain1Cols) > 1):
		assembly[chain1Cols[7]] = chain1Cols[4]
		chain1Cols[6] = float(chain1Cols[6])

	if (len(chain2Cols) > 1 and float(chain2Cols[6]) >= highAbund):
		cdr3 = chain2Cols[4]
		abund = float(chain2Cols[6])
		if (cdr3 not in highAbundCdr3Info):
			highAbundCdr3Info[cdr3] = {}
		highAbundCdr3Info[cdr3][cols[0]] = [1, abund]
	if (len(chain2Cols) > 1):
		assembly[chain2Cols[7]] = chain2Cols[4]
		chain2Cols[6]	 = float(chain2Cols[6])
	barcodeInfo[cols[0]] = {}
	barcodeInfo[cols[0]]["chain1"] = chain1Cols[:]
	barcodeInfo[cols[0]]["chain2"] = chain2Cols[:]
fp.close()

if (args.annot is not None):
	fp = open(args.annot)
	for line in fp:
		header = line.rstrip()
		seq = fp.readline().rstrip() ;
		assemblyId = header.split()[0][1:]
		if (assemblyId in assembly):
			assembly[assemblyId] = seq
	fp.close()

fp = open(barcodeReport)
for line in fp:
	if (line[0] == '#'):
		print(line.rstrip())
		continue

	cols = line.rstrip().split()
	chain1Cols = cols[2].split(',')
	chain2Cols = cols[3].split(',')

	testAgainst = {}
	if (len(chain1Cols) > 1 and float(chain1Cols[6]) < highAbund and chain1Cols[4] in highAbundCdr3Info):
		cdr3 = chain1Cols[4]
		for bc in highAbundCdr3Info[cdr3].keys():
			if (highAbundCdr3Info[cdr3][bc][0] == 0 and highAbundCdr3Info[cdr3][bc][1] * diffuseFrac > float(chain1Cols[6])):
				testAgainst[bc] = 1
	if (len(chain2Cols) > 1 and float(chain2Cols[6]) < highAbund and chain2Cols[4] in highAbundCdr3Info):
		cdr3 = chain2Cols[4]
		for bc in highAbundCdr3Info[cdr3]:
			if (highAbundCdr3Info[cdr3][bc][0] == 1 and highAbundCdr3Info[cdr3][bc][1] * diffuseFrac > float(chain2Cols[6])):
				testAgainst[bc] = 1
	flagFilter = 0
	for bc in testAgainst:
		testChain1 = barcodeInfo[bc]["chain1"]
		testChain2 = barcodeInfo[bc]["chain2"]
		
		for i in [1, 2]:
			c1 = chain1Cols
			c2 = testChain1
			if (i == 2):
				c1 = chain2Cols
				c2 = testChain2
			if (len(c1) > 1 and len(c2) > 1):
				if (c2[6] * diffuseFrac > float(c1[6]) and assembly[c1[7]] in assembly[c2[7]]):
					flagFilter |= i
			elif (len(c2) > 1):
				flagFilter |= i
	if (flagFilter != 3):
 		print(line.rstrip())
fp.close()




