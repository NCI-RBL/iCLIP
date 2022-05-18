#!/usr/bin/python

import sys
import pysam
import pandas as pd

#https://pysam.readthedocs.io/en/latest/api.html

#args
version_id=sys.argv[1]
bamFile=sys.argv[2]
output_dir=sys.argv[3]

#set bam
bamFP=pysam.Samfile(bamFile, "rb")

#set gap minimum
gap_min=100

#format output df
dict = {'Ratio':[],'Length':[],'Min':[],'Max':[],'id':[]}
df = pd.DataFrame(dict)

for read in bamFP:
      if( not( read.is_unmapped ) ):   #if it's mapped
            #read in cigar file, set position of first alignment value
            cigarLine=read.cigar
            pos=0
            
            #for each alignment, search for gaps (IE N=3 skip)
            for (cigarType,cigarLength) in cigarLine:
                if(cigarType==3):
                    #check the alignment before and after
                    #if gap is not proceeded by a mapping, next
                    check_prev_type=cigarLine[pos-1][0]
                    if(pos+1==len(cigarLine)):
                        check_next_type="N"
                    else:
                        check_next_type=cigarLine[pos+1][0]
 
                    #check there is alignments before and after gap and that the gap length is >gap_min
                    if (check_prev_type==0 and check_next_type==0 and cigarLength>gap_min):
                        #get aligned lengths
                        prev_length=cigarLine[pos-1][1]
                        next_length=cigarLine[pos+1][1]

                        #detemerine numerator vs denominator, calc ratio
                        calc_ratio=min(prev_length,next_length)/max(prev_length,next_length)
                        df.loc[len(df.index)] = ["{:.2f}".format(calc_ratio), cigarLength,min(prev_length,next_length),max(prev_length,next_length),version_id] 
            
                #iterate alignment position
                pos=pos+1

output_file=f"{output_dir}/{version_id}.csv"
df.to_csv(output_file, index=False)