# Li Song: the program to re-cluster the CDR3 regions

#!/usr/bin/env python3

import sys
import math

def GetChainType(v, j, c):
	s = ""
	if (v != "*"):
		s = v
	elif (c != "*"):
		s = c 
	elif (j != "*"):
		s = j
	else:
		return -1
	
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
		return -1

def GetMainGeneName(g):
	return g.split("*")[0] 

def CompatibleGeneAssignment(a, b):
	if ( GetChainType( a[0], a[2], a[3] ) != GetChainType(b[0], b[2], b[3]) ):
		return False 
	
	for i in [0, 2]: # V, J gene assignment
		ga = a[i].split("*")[0]
		gb = b[i].split("*")[0]
		if ( ga != "" and gb != "" and ga != gb ):
			return False
	return True

def GetSimilarity(a, b):
	if (len(a) != len(b)):
		return 0
	sameCnt = 0
	for i in range(len(a)):
		if (a[i] == b[i]):
			sameCnt += 1
	return sameCnt / len(a)

def CompatibleSequence(a, b, similarity):
	if (len(a) != len(b)):
		return False ;
	diffCnt = 0 
	diffMax = len(a) - int(math.ceil(len(a) * similarity))
	for i in range(len(a)):
		if (a[i] != b[i]):
			diffCnt += 1
			if (diffCnt > diffMax):
				return False
	return True

def GetFather(tag, father):
	if ( father[tag] != tag ):
		father[tag] = GetFather( father[tag], father )
	return father[tag]
	
def LargerCluster(rawCdr3List, similarity, prefix, useRepresentative, mode):
	vjCDR3LenList = {}
	clusterNameToId = {}
	clusterIdToName = []
	clusterRepresentativeId = {}
	clusterRepresentativeAbund = {}
	
	if (len(rawCdr3List) == 0):
		return
	cdr3List = sorted(rawCdr3List, key=lambda x:(x[0], x[8]))
	# Select the represenative, allowing same CDR3 shown up in different rows in a cluster to handle scRNA-seq
	abund = cdr3List[0][10]
	for i in range(1, len(cdr3List) + 1):
		prevKey = (cdr3List[i - 1][0], cdr3List[i - 1][8])
		key = "*"
		if (i < len(cdr3List)):
			key = (cdr3List[i][0], cdr3List[i][8])
		if (key == prevKey):
			abund += cdr3List[i][10]
		else:
			cdr3 = cdr3List[i - 1]
			if ( cdr3[0] not in clusterNameToId ):
				clusterNameToId[cdr3[0]] = len(clusterIdToName)
				clusterIdToName.append(cdr3[0])
				clusterRepresentativeId[cdr3[0]] = i - 1
				clusterRepresentativeAbund[cdr3[0]] = abund
			elif (abund > clusterRepresentativeAbund[cdr3[0]] ):
				clusterRepresentativeId[cdr3[0]] = i - 1
				clusterRepresentativeAbund[cdr3[0]] = abund
			if (i < len(cdr3List)):
				abund = cdr3List[i][10]
	#print(cdr3List[clusterRepresentativeId["assemble2"]])
	# Reorganize the CDR3 list by their V, J and CDR3 length.
	i = 0
	for cdr3 in cdr3List:
		if (useRepresentative and clusterRepresentativeId[cdr3[0]] != i):
			i += 1
			continue
		key = (GetMainGeneName(cdr3[2]), GetMainGeneName(cdr3[4]), len(cdr3[8]))
		if (key not in vjCDR3LenList):
			vjCDR3LenList[key] = []
		vjCDR3LenList[key].append(i)
		i += 1
	# Prepare the set-union
	father = []
	for cdr3 in cdr3List:
		father.append(clusterRepresentativeId[cdr3[0]])
	
	# Build up the set-union relation
	if (mode == "aggressive"):
		for key in vjCDR3LenList.keys():
			cdr3IdList = vjCDR3LenList[key]
			size = len( cdr3IdList )
			for i in range(size):
				fi = GetFather(cdr3IdList[i], father)
				for j in range(i + 1, size):
					fj = GetFather(cdr3IdList[j], father) 
					if ( fi != fj  
						and CompatibleSequence( cdr3List[ cdr3IdList[i] ][8], cdr3List[ cdr3IdList[j] ][8], 
						similarity ) ):
						#if (cdr3List[cdr3IdList[j]][8] == "TGCCAACAGTATATTAGTTACTCGTACACTTTT"):
						#	print(i, cdr3List[cdr3IdList[i]])
						#	print(j, cdr3List[cdr3IdList[j]])
						father[fj] = fi
	elif (mode == "center"):
		for key in vjCDR3LenList.keys():
			rawCdr3IdList = vjCDR3LenList[key][:]
			cdr3IdList = sorted(rawCdr3IdList, key=
					lambda x: (clusterRepresentativeAbund[cdr3List[x][0]], cdr3List[x][10]), reverse=True)
			size = len(cdr3IdList)
			for i in range(1, size):
				maxFj = 0
				maxSimilarity = -1
				fi = GetFather(cdr3IdList[i], father)
				for j in range(i):
					fj = GetFather(cdr3IdList[j], father)
					if (fi == fj):
						continue
					s = GetSimilarity(cdr3List[fi], cdr3List[fj])
					if (s > maxSimilarity):
						maxSimilarity = s
						maxFj = fj
				if (maxSimilarity >= similarity):
					father[fi] = maxFj

	# Associate original clusters to larger cluster
	largerClusterToId = []
	largerClusterToClusterName = []
	rootToLargerClusterId = {}
	i = 0 
	for cdr3 in cdr3List:
		root = GetFather(i, father)
		lcId = 0
		if (root not in rootToLargerClusterId):
			rootToLargerClusterId[root] = len(largerClusterToId)
			lcId = len(largerClusterToId)
			largerClusterToId.append([])
			largerClusterToClusterName.append(set({}))
		else:
			lcId = rootToLargerClusterId[root]

		largerClusterToId[lcId].append(i)
		largerClusterToClusterName[lcId].add(cdr3[0])
		i += 1		
	
	# Output the result
	# Output the composition of larger cluster
	#for i in range(len(largerClusterToId)):
	#	print("#\tcluster" + str(i), end = "")
	#	for name in largerClusterToClusterName[i]:
	#		print("\t" + name, end = "")
	#	print("\n", end = "")
	
	# Output the new cdr3 format with new larger cluster Id.
	for i in range(len(largerClusterToId)):
		j = 0
		for cdr3Id in largerClusterToId[i]:
			cdr3List[cdr3Id].append(cdr3List[cdr3Id][0])
			cdr3List[cdr3Id].append(cdr3List[cdr3Id][1])
			cdr3List[cdr3Id][0] = prefix + "_" + str(i) #+ "_" + cdr3List[cdr3Id][0] + "_" + str(cdr3List[cdr3Id][1])
			cdr3List[cdr3Id][1] = j
			print( "\t".join( str(x) for x in cdr3List[cdr3Id] ) )
			j += 1

	return 


if (__name__ == "__main__"):
	if (len(sys.argv) <= 1):
		print("usage: a.py trust_cdr3.out [OPTIONS] > output\n" + 
			"OPTIONS:\n" +
			"\t-s FLOAT: similarity of two CDR3s (default: 0.8)\n" +
			"\t--prefix STRING: prefix to new cluster name (default: cluster)\n" + 
			"\t--center: use the center of the cluster for similarity comparison (default: no)\n"+
			"\t--representative: use representative CDR3 from each contig for cluster (default: no)\n" +
			"\t--format [cdr3, simplerep]: the input format type (default: cdr3)")
		exit(1)

	cdr3List = []
	similarity = 0.8 ;
	prefix = "cluster"
	useRepresentative = False
	mode = "aggressive"
	inputFormat = "cdr3"
	i = 2
	while (i < len(sys.argv)):
		if (sys.argv[i] == "-s"):
			similarity = float(sys.argv[i + 1])
			i += 1
		elif (sys.argv[i] == "--prefix"):
			prefix = sys.argv[i + 1]
			i += 1
		elif (sys.argv[i] == "--representative"):
			useRepresentative = True
		elif (sys.argv[i] == "--center"):
			mode = "center"
		elif (sys.argv[i] == "--format"):
			inputFormat = sys.argv[i + 1]
			if (inputFormat not in ["cdr3", "simplerep"]):
				print("Unknown format: ", inputFormat)
				exit(1)
			i += 1
		else:
			print("Unknown option: ", sys.argv[i])
			exit(1)
		i += 1

	fp = open(sys.argv[1])
	lineCnt = 0
	for line in fp:
		line = line.rstrip()
		cols = line.split("\t")
		if (inputFormat == "cdr3"):
			cols[1] = int(cols[1])
			skip = False
			for g in [2, 4]: # Must have V, J genes.
				if ( cols[g] == "*"):
					skip = True 
			if (float(cols[9]) == 0):
				skip = True
			if ( skip ):
				continue 
			for	g in [2, 3, 4, 5]:
				cols[g] = cols[g].split(",")[0]
			cols[9] = float(cols[9])
			cols[10] = float(cols[10])
			if (cols[9] == 0):
				continue
		elif (inputFormat == "simplerep"):
			if (line[0] == "#"):
				continue
			if ("_" in cols[3] or "?" in cols[3]):
				continue
			for g in [4, 6]: # Must have V, J genes.
				if ( cols[g] == "*"):
					skip = True 
			reformat = [0] * 11
			reformat[0] = "line" + str(lineCnt)
			reformat[1] = 0
			for	g in [4, 5, 6, 7]:
				reformat[g - 2] = cols[g]	
			reformat[6] = reformat[7] = "*"
			reformat[8] = cols[2]
			reformat[9] = 1
			reformat[10] = cols[0]
			cols = reformat[:]
		cdr3List.append(cols)
		lineCnt += 1
	fp.close()
	LargerCluster(cdr3List, similarity, prefix, useRepresentative, mode)
