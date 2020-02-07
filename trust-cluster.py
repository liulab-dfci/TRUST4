# Li Song: the program to re-cluster the CDR3 regions

#!/usr/bin/env python3

import sys

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

def CompatibleSequence(a, b, similarity):
	if (len(a) != len(b)):
		return False ;
	diffCnt = 0 
	diffMax = len(a) - int(len(a) * similarity)
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
	
def LargerCluster(cdr3List, similarity):
	vjCDR3LenList = {}
	clusterNameToId = {}
	clusterIdToName = []
	for cdr3 in cdr3List:
		if ( cdr3[0] not in clusterNameToId ):
			clusterNameToId[cdr3[0]] = len(clusterIdToName)
			clusterIdToName.append(cdr3[0])

	# Reorganize the CDR3 list by their V, J and CDR3 length.
	i = 0
	for cdr3 in cdr3List:
		key = (GetMainGeneName(cdr3[2]), GetMainGeneName(cdr3[4]), len(cdr3[8]))
		if (key not in vjCDR3LenList):
			vjCDR3LenList[key] = []
		vjCDR3LenList[key].append(i)
		i += 1
	# Prepare the set-union
	father = []
	clusterNameFirstId = {}
	i = 0
	for cdr3 in cdr3List:
		father.append(i)
		if (cdr3[0] in clusterNameFirstId ):
			father[i] = clusterNameFirstId[cdr3[0]]
		else:
			clusterNameFirstId[cdr3[0]] = i
		i += 1
	
	# Build up the set-union relation
	for key in vjCDR3LenList.keys():
		cdr3IdList = vjCDR3LenList[key]
		size = len( cdr3IdList )
		for i in range(size):
			for j in range(i + 1, size):
				if ( GetFather( cdr3IdList[i], father ) != GetFather( cdr3IdList[j], father ) 
					and CompatibleSequence( cdr3List[ cdr3IdList[i] ][8], cdr3List[ cdr3IdList[j] ][8], 
					similarity ) ):
					father[ cdr3IdList[j] ] = cdr3IdList[i]
					
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
			cdr3List[cdr3Id][0] = "cluster" + str(i) #+ "_" + cdr3List[cdr3Id][0] + "_" + str(cdr3List[cdr3Id][1])
			cdr3List[cdr3Id][1] = j
			print( "\t".join( str(x) for x in cdr3List[cdr3Id] ) )
			j += 1

	return 


if (__name__ == "__main__"):
	if (len(sys.argv) <= 1):
		print("usage: a.py trust_cdr3.out [OPTIONS]\n" + 
			"OPTIONS:\n" 
			"\t-s FLOAT: similarity of two CDR3s (default: 0.95)") 	
		exit(1)
	cdr3List = []
	similarity = 0.95 ;
	i = 2
	while (i < len(sys.argv)):
		if (sys.argv[i] == "-s"):
			similarity = float(sys.argv[i + 1])
			i += 1
		else:
			print("Unknown option: ", sys.argv[i])
			exit(1)
		i += 1

	fp = open(sys.argv[1])
	for line in fp:
		line = line.rstrip()
		cols = line.split("\t")
		cols[1] = int(cols[1])
		skip = False
		for g in [2, 4]: # Must have V, J genes.
			if ( cols[g] == "*"):
				skip = True 
		if (float(cols[9]) == 0):
			skip = True
		if ( skip ):
			continue 

		for g in [2, 3, 4, 5]:
			cols[g] = cols[g].split(",")[0]
		cols[9] = float(cols[9])
		cols[10] = float(cols[10])
		if (cols[9] == 0):
			continue
		cdr3List.append(cols)
	LargerCluster(cdr3List, similarity)
