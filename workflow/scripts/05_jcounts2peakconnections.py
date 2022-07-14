#!/usr/bin/env python3
# arguments
# arg1 site-to-peak lookup tsv
# arg2 jcounts file
# output
# connected peak clusters ... 1 per line
import sys
import networkx as nx

site2peak=dict()
with open(sys.argv[1],'r') as fh:
	Lines = fh.readlines()
	for line in Lines:
		line=line.strip().split("\t")
		site2peak[line[0]]=line[1]

G = nx.Graph()

with open(sys.argv[2],'r') as fh:
	Lines = fh.readlines()
	for line in Lines:
		line=line.strip().split("\t")
		peakid=line[0]
		strand=peakid.split("_")[-1]
		site1=line[2]+"_"+line[3]+"_"+strand
		site2=line[5]+"_"+line[6]+"_"+strand
		try:
			peak1=site2peak[site1]
		except KeyError:
		#	print(site1)
			continue
		try:
			peak2=site2peak[site2]
		except KeyError:
		#	print(site2)
			continue
		#print(peak1+"\t"+peak2)
		G.add_node(peak1)
		G.add_node(peak2)
		G.add_edge(peak1,peak2)

for i in nx.connected_components(G):
	print("\t".join(list(i)))