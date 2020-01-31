# Li Song: the program to show the 

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

# return the range on the s corresponds 
def GetNodeRange(s, pos):
	i = pos
	while (i < len(s)):
		if (s[i] == "," or s[i] == ')' or s[i] == ';' ):
			break 
		i += 1
	return pos, i - 1 
	
def UpdateTree(s, childList, tree):
	parts = s.split(":") 
	cols = parts[0].split("_")
	abund = float(cols[-1])
	subId = int(cols[-2])
	clusterId = "_".join(cols[0:-2])
	distToParent = -1
	if (len(parts) > 1):
		distToParent = int(parts[1])
	
	node = (clusterId, subId)
	if (node not in tree):
		tree[node] = {"children":[], "abund":abund, "parent":None, "distToParent":distToParent}
	for n in childList:
		tree[node]["children"].append( n )
		tree[n]["parent"] = node
	return node

# Given the Newick representation between [start, end], update the tree structure 
def ParseTree(s, pos, tree):
	brotherList = []
	while (pos < len(s)):
		# Go through the nodes on this level
		childList = []
		if (s[pos] == '('):
			childList, pos = ParseTree(s, pos + 1, tree)
		start, end = GetNodeRange(s, pos)
		newNode = UpdateTree(s[start:end + 1], childList, tree)
		brotherList.append(newNode)
		pos = end + 1
		
		if (s[pos] == ','): # End of this branch 
			pos += 1
			continue
		else:
			pos += 1
			break 
	return brotherList, pos
			
def PrintTree(tree):
	for n in tree.keys():
		parent = tree[n]["parent"]
		if (parent == None):
			parent = (n[0], n[1])
		print( "\t".join([n[0], str(n[1]), "%.2f"%tree[n]["abund"], parent[0], str(parent[1]), str(tree[n]["distToParent"])]) )

if (__name__ == "__main__"):
	if (len(sys.argv) <= 1):
		print("usage: a.py trust_evo_tree.out [OPTIONS]\n" + 
			"OPTIONS:" ) 
		exit(1)
	fp = open(sys.argv[1])
	for line in fp:
		tree = {} # key: (clusterid, subid)
		ParseTree( line.rstrip(), 0, tree )
		PrintTree( tree )

