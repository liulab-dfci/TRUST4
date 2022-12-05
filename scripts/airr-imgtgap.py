#!/usr/bin/env python3

# Add gaps from the IMGT reference to the AIRR sequence and germline alignment field

import argparse
import re

def ParseCigar(cigar):
	cigarFields = re.findall("\d+\w", cigar)
	ret = []
	for f in cigarFields:
		ret.append( (int(f[0:-1]), f[-1]) )
	return ret

def InsertGap(seq, gaps):
	if (len(gaps) > 1):
		subseqs = []
		subseqs.append(seq[0:gaps[0][0]+1])
		for i in range(1, len(gaps)):
			start = gaps[i - 1][0] + 1
			end = gaps[i][0] + 1 # open end

			subseqs.append(seq[start:end])
		subseqs.append( seq[(gaps[-1][0]+1):] )
		
		ret = ""
		for i in range(len(gaps)):
			ret += subseqs[i] + ("."*(gaps[i][1]))
		ret += subseqs[-1]
		return ret 
	else:
		return seq

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description="Add IMGT gaps to the sequence and germline alignment field in the AIRR format. Output to stdout.")
	parser.add_argument("-i", help="IMGT file", dest="imgtFile", required=True)
	parser.add_argument("-a", help="AIRR file", dest="airrFile", required=True)

	args = parser.parse_args()

	# Read in IMGT sequuences
	fp = open(args.imgtFile)
	seq = ""
	gene = ""
	imgtSeq = {}
	imgtSeqWoGap = {}
	imgtSeqGapInfo = {} # a gap after the position in imgtSeqWoGap

	for line in fp:
		line = line.rstrip()
		if (line[0] == ">"):
			if (gene != ""):
				imgtSeq[gene] = seq			
			gene = line[1:].split()[0]
			seq = ""
		else:
			seq += line
	fp.close()
	# Mark the gaps
	for gene, seq in imgtSeq.items(): 
		imgtSeqGapInfo[gene] = []
		matches = re.finditer("(\.+)", seq)
		#imgtSeqGapInfo[gene] = [m.span() for m in matches])

		# Convert the information into gapless sequence
		imgtSeqWoGap[gene] = imgtSeq[gene].replace('.', '')
		psum = 0 # partial sum
		for m in matches:
			span = m.span()
			imgtSeqGapInfo[gene].append((span[0] - psum - 1, span[1] - span[0]))
			psum += span[1] - span[0] 

	# Process the AIRR data 
	fp = open(args.airrFile)
	header = fp.readline()
	header = header.rstrip()
	cols = header.split("\t")
	colNameId = {}
	for i,c in enumerate(cols):
		colNameId[c] = i
	print(header)

	for line in fp:
		cols = line.rstrip().split("\t")
		if (len(cols[colNameId["v_call"]]) >= 4 and len(cols[colNameId["v_cigar"]]) > 0):
			seq = cols[ colNameId["sequence_alignment"] ]
			germline = cols[ colNameId["germline_alignment"]]
			cigarV = cols[ colNameId["v_cigar"] ]
			gene = cols[ colNameId["v_call"] ]
			geneLength = len(imgtSeqWoGap[gene])
			cigarFields = ParseCigar(cigarV)

			seqStart = 0	
			germlineStart = 0
			germlineEnd = 0 # inclusive end
			for i in [0, 1]:
				if (len(cigarFields) <= i):
					continue
				if (cigarFields[i][1] == 'N'):
					germlineStart = cigarFields[i][0]
					
			for i in range(len(cigarFields) - 1, -1, -1):
				if (cigarFields[i][1] == "N"):
					germlineEnd = geneLength - 1 - cigarFields[i][0]
				elif (cigarFields[i][1] == "S"):
					continue
				else:
					break
			tag = 0 # the gap id
			gaps = imgtSeqGapInfo[gene]
			for tag in range(len(gaps)):
				if (gaps[tag][0] >= germlineStart):
					break

			i = 0 # alignment position for both sequence_alignment and germline_alignment
			j = germlineStart # the position within the imgt v gene sequence without gap
			insertGaps = []
			while (i < len(germline) - 1 and j < geneLength and tag < len(gaps)):
				if (germline[i] != "-"):
					if (j == gaps[tag][0]):
						insertGaps.append((i, gaps[tag][1]))
						tag += 1
					j += 1
				i += 1
			cols[ colNameId["sequence_alignment"] ] = InsertGap(seq, insertGaps)
			cols[ colNameId["germline_alignment"] ] = InsertGap(germline, insertGaps)
			#<-if there is v hit
		print("\t".join(cols))
	fp.close()
