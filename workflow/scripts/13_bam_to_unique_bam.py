# Author: Vishal Koparde, PhD
# Date: Jan 2022
import pysam
import argparse

parser = argparse.ArgumentParser(description='Filter BAM by readids')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
args = parser.parse_args()
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)
bigdict = dict()

count=0
for read in inBAM.fetch():
    if read.is_unmapped:
        continue
    count+=1
    if count%1000000 == 0:
        print("%d reads read!"%(count))
    try:
        nh=read.get_tag("NH")
    except KeyError:
        continue
    if nh==1:
        qn=read.query_name
        if not qn in bigdict:
            bigdict[qn]=list()
        bigdict[qn].append(read)
inBAM.close()

for r in bigdict.keys():
  try:
    for read in bigdict[r]:
      outBAM.write(read)
  except KeyError:
    continue
outBAM.close()