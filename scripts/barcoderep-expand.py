#!/usr/bin/env python3 

import sys
import argparse

parser = argparse.ArgumentParser(description="Create new entries for the secondary chains in the barocde_report file. This script Will rename the barcode.")
parser.add_argument('-b', help="", required=True, dest="barcode_report")
parser.add_argument('--chain', help="expand chain1 or chain2", default=1, dest="chain")
parser.add_argument('--frac', help="the abundance needs to be more the fraction of the primary chain", default=0.1, dest="frac")

args = parser.parse_args()

barcodeReport = args.barcode_report
chain = int(args.chain)
frac = float(args.frac)

def GetChainType(v, j, c):
	s = ""
	if (c != "*" and c != "."):
		s = c 
	elif (j != "*" and j != "."):
		s = j
	elif (v != "*" and v != "."):
		s = v
	else:
		return 7
	
	if (s[0:3] == "IGH"):
		return 0
	elif (s[0:3] == "IGK"):
		return 1
	elif (s[0:3] == "IGL"):
		return 2
	elif (s[0:3] == "TRA"):
		return 3
	elif (s[0:3] == "TRB"):
		return 4
	elif (s[0:3] == "TRG"):
		return 5
	elif (s[0:3] == "TRD"):
		return 6
	else:
		return 7


def GetCellType(v, j, c, defaultType = "*"):
	chainType = GetChainType(v, j, c)
	if (chainType <= 2):
		return "B"
	elif (chainType <= 4):
		return "abT"
	elif (chainType <= 6):
		return "gdT"
	else:
		return defaultType

#AGAGTGGTCTATCCTA-1	abT	TRBV6-6*01,TRBD2*01,TRBJ2-5*01,*,TGTGCCAGTCTACTTGGGGGGACCCAGTACTTC,CASLLGGTQYF,3.61,AGAGTGGTCTATCCTA-1_48281,100.00,0	TRAV10*01,*,TRAJ34*01,*,TGTGTGGTGAGCGCCCGCACCGACAAGCTCATCTTT,CVVSARTDKLIF,1.00,AGAGTGGTCTATCCTA-1_65563,100.00,0	*	*
fp = open(barcodeReport)
for line in fp:
	if (line[0] == "#"):
		print(line.rstrip())
		continue

	cols = line.rstrip().split() 
	barcode = cols[0]

	# Output the primary chain
	outputCols = cols[:]
	outputCols[0] = barcode + "_0"
	print("\t".join(outputCols))
	
	# Expand the secondary chain
	secondaryEntry = cols[3 + chain]
	if (cols[1 + chain] == "*" or secondaryEntry == "*"):
		continue
	primaryAbund = float(cols[1 + chain].split(',')[6])
	secondaryCols = secondaryEntry.split(";")
	for i in range(2,len(outputCols)):
		outputCols[i] = "*"

	k = 0
	for c in secondaryCols:
		outputCols[0] = barcode + "_" + str(k + 1)
		outputCols[chain + 1] = c
		subCols = c.split(',')
		abund = float(subCols[6])
		outputCols[1] = GetCellType(subCols[0], subCols[2], subCols[3], cols[1])
		if (abund < primaryAbund * frac):
			continue
		print("\t".join(outputCols))
		k += 1
fp.close()
