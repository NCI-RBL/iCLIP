# Author: Vishal Koparde, PhD
# Date: Apr 2021
import pysam
import sys
import argparse
import os
class Readinfo:
	def __init__(self,read,save_read=False):
		def _get_bitflag(r):
			bitflag=str(r).split("\t")[1]
			return str(bitflag)
		self.bitflag=_get_bitflag(read)
		self.mapq=int(read.mapping_quality)
		self.nhtag=str(read.get_tag("NH"))
		self.ref=read.reference_name
		self.cigar=read.cigarstring
		self.pos=read.reference_start
		if save_read:
			self.alignment=read
		else:
			self.alignment=""
		self.refpos=list()

	def __str__(self):
		return "%s\t%d\t%s\t%s\t%s\t%s"%(self.bitflag,self.mapq,self.nhtag,self.cigar,self.ref,self.pos)
		
	def is_spliced(self):
		cigart=self.alignment.cigartuples
		if 3 in list(map(lambda z:z[0],cigart)):
			return True
		else:
			return False
	
	def get_reference_positions(self):
		self.refpos=set(filter(lambda x:x!=None,self.alignment.get_reference_positions(full_length=True)))
		
def read_bam(bamfilename,save_alignments=False):
	print("Reading BAM:%s"%(bamfilename))
	bam=pysam.AlignmentFile(bamfilename, "rb")
	bigdict = dict()

	count=0
	for read in bam.fetch():
		count+=1
		if count%1000000 == 0:
			print("%d reads read from file: %s!"%(count,bamfilename))
		qn=read.query_name
		if not qn in bigdict:
			bigdict[qn]=dict()
			bigdict[qn]['reads']=list()
			bigdict[qn]['maxmapq']=-1
		bigdict[qn]['reads'].append(Readinfo(read,save_alignments))
		if int(read.mapping_quality) > bigdict[qn]['maxmapq']:
			bigdict[qn]['maxmapq']=int(read.mapping_quality)
	bam.close()
	return bigdict


def find_overlapping_read_mapq(read,readlist):
	mapq=-1
	refpos=set(filter(lambda x:x!=None,read.alignment.get_reference_positions(full_length=True)))
	for r in readlist:
		if r.ref != read.ref:
			continue
		refpos1=set(filter(lambda x:x!=None,r.alignment.get_reference_positions(full_length=True)))
		if len(refpos1.intersection(refpos))>=5:
			mapq=r.alignment.mapping_quality
			break
	return mapq
		

def main():
	parser = argparse.ArgumentParser(
formatter_class=argparse.RawDescriptionHelpFormatter,
description='''
Inputs:

1. inputBAM1 --> Novoalign bam file (splicing unaware) with genome as reference
2. inputBAM2 --> Novoalign bam file (splicing aware) with exon-masked genome + transcriptome as reference
3. inputBAM3 --> Novoalign bam file (splicing aware) with unmasked genome + transcriptome as reference

For inputBAM2 and inputBAM3, transcriptome alignments have already been converted back to genomic coordinates

Outputs:

1. outBAM --> output BAM with corrected MAPQ scores. This uses inputBAM3 as a template for header etc.
2. outTSV --> output with read-by-read mapping metadata eg.

*******
NS500326:331:H53GGBGX5:2:21207:8400:7107:rbc:AGTAT
File1:maxmapq:3
16      3       1       29H34M1I3M      chr1    171377340
File2:maxmapq:3
16      3       1       29H34M1I3M      chr1    171377340
File3:maxmapq:3
16      3       1       29H34M1I3M      chr1    171377340
*******

The tab-delimited columns for each alignment reported are SAM bitflag, MAPQ, NH, CIGAR, chromosome, start position.
In addition to the read-by-read metadata TSV, MAPQ corrected reads are also provided to stdout eg.

*******
MAPQ changed from 3 to 70:NS500326:331:H53GGBGX5:4:11410:3852:5796:rbc:CCCAA
MAPQ changed from 46 to 53:NS500326:331:H53GGBGX5:1:22211:14869:11391:rbc:CGGCG
MAPQ changed from 37 to 38:NS500326:331:H53GGBGX5:2:22103:10319:4954:rbc:AAACC
*******

The recommended way to run this script is something like this:

% python correct_mapq.py --inputBAM1 A.bam --inputBAM2 B.bam --inputBAM3 C.bam --out out.tsv --outBAM C_mapq_updated.bam > out.log 2>&1

If the input BAMs are large with multimappings, then this script does tend to use significantly large amount of memory.
Using --partition=largemem with 1TB of memory request to slurm is recommended.
''')
	parser.add_argument('--inputBAM1', dest='inputBAM1', type=str, required=True,
						help='input BAM file1 (A.bam)')
	parser.add_argument('--inputBAM2', dest='inputBAM2', type=str, required=True,
						help='input BAM file2 (B.bam)')
	parser.add_argument('--inputBAM3', dest='inputBAM3', type=str, required=True,
						help='input BAM file3 (C.bam)')
	parser.add_argument('--outBAM', dest='outBAM', type=str, required=True,
						help='output BAM file (out.bam)')
	parser.add_argument('--out', dest='out', type=str, required=True,
						help='outTSV with read-by-read metadata')
	args = parser.parse_args()

	bigdict1 = read_bam(args.inputBAM1,True)
	bigdict2 = read_bam(args.inputBAM2,True)
	bigdict3 = read_bam(args.inputBAM3,True)
	print("Done reading BAMs!")


	keys=[]
	keys.extend(bigdict1.keys())
	keys.extend(bigdict2.keys())
	keys.extend(bigdict3.keys())
	keys=set(keys)
	dummy="%d\t%d\t%d\t%d\t%d\t%d"%(-1,-1,-1,-1,-1,-1)
	o=open(args.out,'w')
	for key in keys:
		o.write("%s\n"%(key))
		if key in bigdict1:
			o.write("File1:maxmapq:%d\n"%(bigdict1[key]['maxmapq']))
			for i in bigdict1[key]['reads']:
				o.write("%s\n"%(i))
		else:
			bigdict1[key]=dict()
			bigdict1[key]['reads']=list()
			bigdict1[key]['maxmapq']=-1
			o.write("%s\n"%(dummy))
		if key in bigdict2:
			o.write("File2:maxmapq:%d\n"%(bigdict2[key]['maxmapq']))
			for i in bigdict2[key]['reads']:
				o.write("%s\n"%(i))
		else:
			bigdict2[key]=dict()
			bigdict2[key]['reads']=list()
			bigdict2[key]['maxmapq']=-1
			o.write("%s\n"%(dummy))
		if key in bigdict3:
			o.write("File3:maxmapq:%d\n"%(bigdict3[key]['maxmapq']))
			for i in bigdict3[key]['reads']:
				o.write("%s\n"%(i))
		else:
			bigdict3[key]=dict()
			bigdict3[key]['reads']=list()
			bigdict3[key]['maxmapq']=-1
			o.write("%s\n"%(dummy))
		o.write("\n")
	o.close()
	
	print("Done writing %s!"%(args.out))

	inbam = pysam.AlignmentFile(args.inputBAM3, "rb" )
	outbam = pysam.AlignmentFile(args.outBAM, "wb", template=inbam )
	inbam.close()
	
	for rid in bigdict3.keys():
		if bigdict3[rid]['maxmapq']<=0 and bigdict1[rid]['maxmapq']<=0 and bigdict2[rid]['maxmapq']<=0: # read is multimapped and all MAPQs are zero OR maxmapq is -1 ... readid absent in other BAM file(s)
			for r in bigdict3[rid]['reads']:
				outbam.write(r.alignment)
		elif bigdict1[rid]['maxmapq']!=-1 and bigdict2[rid]['maxmapq']!=-1 and bigdict3[rid]['maxmapq']!=-1: # readid present in all 3 files
			if (bigdict1[rid]['maxmapq'] > bigdict3[rid]['maxmapq']) or (bigdict2[rid]['maxmapq'] > bigdict3[rid]['maxmapq']): # better MAPQ may be available
				spliced_reads=[]
				for i,r in enumerate(bigdict3[rid]['reads']):
					if not r.is_spliced:
						outbam.write(r)
					else:
						spliced_reads.append(i)
				for i in spliced_reads:
					read=bigdict3[rid]['reads'][i]
					mq3=read.mapq
					mq1=mq2=-1
					if bigdict1[rid]['maxmapq']>mq3: # file1 MAY have higher quality overlap
						mq1=find_overlapping_read_mapq(read,bigdict1[rid]['reads'])
					if bigdict2[rid]['maxmapq']>mq3: # file2 MAY have higher quality overlap
						mq2=find_overlapping_read_mapq(read,bigdict2[rid]['reads'])
					newmq=mq3
					if mq1>newmq:
						newmq=mq1
					if mq2>newmq:
						newmq=mq2
					read.alignment.mapping_quality=newmq
					if newmq!=mq3:
						print("MAPQ changed from %d to %d:%s"%(mq3,newmq,read.alignment.query_name))
					outbam.write(read.alignment)
			else:
				for r in bigdict3[rid]['reads']:
					outbam.write(r.alignment)
		elif bigdict3[rid]['maxmapq']!=-1 and bigdict2[rid]['maxmapq']!=-1 and bigdict1[rid]['maxmapq']==-1: # readid absent in first file
			if bigdict2[rid]['maxmapq'] > bigdict3[rid]['maxmapq']:
				spliced_reads=[]
				for i,r in enumerate(bigdict3[rid]['reads']):
					if not r.is_spliced:
						outbam.write(r.alignment)
					else:
						spliced_reads.append(i)
				for i in spliced_reads:
					read=bigdict3[rid]['reads'][i]
					mq3=read.mapq
					mq2=-1
					if bigdict2[rid]['maxmapq']>mq3: # file2 MAY have higher quality overlap
						mq2=find_overlapping_read_mapq(read,bigdict2[rid]['reads'])
					newmq=mq3
					if mq2>newmq:
						newmq=mq2
					read.alignment.mapping_quality=newmq
					if newmq!=mq3:
						print("MAPQ changed from %d to %d:%s"%(mq3,newmq,read.alignment.query_name))
					outbam.write(read.alignment)				
			else:
				for r in bigdict3[rid]['reads']:
					outbam.write(r.alignment)
		elif bigdict3[rid]['maxmapq']!=-1 and bigdict1[rid]['maxmapq']!=-1 and bigdict2[rid]['maxmapq']==-1: # readid absent in second file
			if bigdict1[rid]['maxmapq'] > bigdict3[rid]['maxmapq']:
				spliced_reads=[]
				for i,r in enumerate(bigdict3[rid]['reads']):
					if not r.is_spliced:
						outbam.write(r)
					else:
						spliced_reads.append(i)
				for i in spliced_reads:
					read=bigdict3[rid]['reads'][i]
					mq3=read.mapq
					mq1=-1
					if bigdict1[rid]['maxmapq']>mq3: # file2 MAY have higher quality overlap
						mq1=find_overlapping_read_mapq(read,bigdict2[rid]['reads'])
					newmq=mq3
					if mq1>newmq:
						newmq=mq1
					read.alignment.mapping_quality=newmq
					if newmq!=mq3:
						print("MAPQ changed from %d to %d:%s"%(mq3,newmq,read.alignment.query_name))
					outbam.write(read.alignment)				
			else:
				for r in bigdict3[rid]['reads']:
					outbam.write(r.alignment)
		else: # readid absent in first and second file
			for r in bigdict3[rid]['reads']:
				outbam.write(r.alignment)
					
	outbam.close()



if __name__ == "__main__":
    main()
