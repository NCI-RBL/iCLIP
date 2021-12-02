# Author: Vishal Koparde, PhD
# Date: Jan 2021
import pysam
import sys
import argparse
import os
parser = argparse.ArgumentParser(description='Filter BAM by readids')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (one readid per line)')
args = parser.parse_args()
rids = list(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)
bigdict = dict()

count=0
for read in inBAM.fetch():
	count+=1
	if count%1000000 == 0:
		print("%d reads read!"%(count))
	qn=read.query_name
	if not qn in bigdict:
		bigdict[qn]=list()
	bigdict[qn].append(read)
inBAM.close()

for r in rids:
  try:
    for read in bigdict[r]:
      outBAM.write(read)
  except KeyError:
    continue
outBAM.close()
