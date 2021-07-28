#!/usr/bin/env python3

import sys
import argparse
import math

isotypeRanks = {"IGHM":0, "IGHD":1, "IGHG3":2, "IGHG1":3, "IGHA1":4, "IGHG2":5, "IGHG4":6, "IGHE":7, "IGHA2":8, "*":9, ".":9}
isotypeOrder = ["IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2"]
chainOrder = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]

def GetChainType(v, j, c):
	s = ""
	if (c != "*" and c != "."):
		s = c 
	elif (j != "*" and j != "."):
		s = j
	elif (v != "*" and v != "."):
		s = v
	else:
		return -1
	
	if (s[0:3] == "IGH"):
		return (0, isotypeRanks[c])
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
		return -1

def ComputeRichness(rep):
	return len(rep)

def ComputeCPK(rep):
	if (len(rep) == 0):
		return "NA"
	return len(rep) / sum(list(rep.values())) * 1000

def ComputeEntropy(rep):
	if (len(rep) == 0):
		return "NA"
	total = sum(list(rep.values()))
	return sum([-x/total*math.log(x/total) for x in rep.values()]) 

def ComputeClonality(rep):
	if (len(rep) <= 1):
		return "NA"
	return 1 - ComputeEntropy(rep) / math.log(len(rep))

def OutputChain(rep, name):
	outputList = [name]
	outputList.append(sum(list(rep.values())))
	outputList.append(ComputeRichness(rep))
	outputList.append(ComputeCPK(rep))
	outputList.append(ComputeEntropy(rep))
	outputList.append(ComputeClonality(rep))
	print("\t".join([str(x) for x in outputList]))

def ProcessImmuneRepertoire(immrep):
	print("\t".join(["#chain", "Abundance", "Richness", "CPK", "Entropy", "Clonality"]))
	# IGH overall
	tmp = {}
	for i in range(10):
		for c in immrep[(0,i)]:
			if (c not in tmp):
				tmp[c] = 0
			tmp[c] += immrep[(0,i)][c]
	OutputChain(tmp, "IGH")
	# IGH chains
	for i in range(9):
		OutputChain(immrep[(0, i)], isotypeOrder[i])
	for i in range(1, 7):
		OutputChain(immrep[i], chainOrder[i])

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description="Immune repertoire diversity statistics")
	parser.add_argument("-r", help="repertoire file", dest="repfile", required=True)
	parser.add_argument("-f", help="repertoire file format (TRUST4_report, TRUST4_barcode_report,...)", dest="format", default="TRUST4_report")
	parser.add_argument("--ntaa", help="Use nucleotide(nt) or amino acids(aa)", dest="ntaa", default="aa")
	args = parser.parse_args()

	immrep = {}
	for i in range(10):
		immrep[(0, i)] = {}
	for i in range(1, 7):
		immrep[i] = {}

	fp = open(args.repfile)
	if (args.format == "TRUST4_report"):
		for line in fp:
			if (line[0] == "#" or line[0:5] == "count"):
				continue
			cols = line.rstrip().split()
			chainType = GetChainType(cols[4], cols[6], cols[7])
			if ("_" in cols[3] or cols[3] == "partial" or "?" in cols[3] or 
				chainType == -1):
				continue
			if (cols[3] not in immrep[chainType]):
				immrep[chainType][cols[3]] = 0
			immrep[chainType][cols[3]] += int(cols[0])
	elif (args.format == "TRUST4_barcode_report"):
		for line in fp:
			if (line[0] == "#" or line[0:5] == "count"):
				continue
			mainCols = line.rstrip().split()
			for i in [2, 3]:
				if (mainCols[i] == "*"):
					continue
				cols = mainCols[i].split(",")
				chainType = GetChainType(cols[0], cols[2], cols[3])
				if ("_" in cols[5] or cols[5] == "partial" or "?" in cols[5] or
						chainType == -1):
					continue
				if (cols[5] not in immrep[chainType]):
					immrep[chainType][cols[5]] = 0
				immrep[chainType][cols[5]] += 1
	else:
		print("Unknown format " + args.format)
		sys.exit(1)
	fp.close()

	ProcessImmuneRepertoire(immrep) 
