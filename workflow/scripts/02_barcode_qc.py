import matplotlib.pyplot as plt
import HTSeq
import itertools
from collections import Counter
import pandas as pd
import sys
import json

#NOTE
#Current script is set up for two barcoding strategies - if additional are required they 
# must be added below. Current accepted formats include:
#NNNNNCACTGTNNNN = total len 15, seq len 6, position [5:11]
#NNNTGGCNN = total len 9, seq len 4, position [3:7]

#args
sample_manifest = sys.argv[1] #"samples_corrected.tsv"
multiplex_manifest = sys.argv[2] #"multiplex_corrected.tsv"
barcode_input = sys.argv[3]
output_dir = sys.argv[4] #location input barcode counts, output barcode plots and summary
mismatch = sys.argv[5] #2

def CreateSampleDicts(df_m,df_s):
    s_dict={}
    m_dict={}

    #create dict that maps multiplex: sampleid : barcode
    for multi_name in df_s['multiplex'].unique():
        sub_df = df_s[df_s['multiplex']==multi_name]
        s_dict[multi_name]={}

        #each row maps to a sample:barcode and sample:group
        for index,row in sub_df.iterrows():
            sample_name = row['sample']
            s_dict[multi_name][sample_name]=row['barcode']

    #creat dict that maps mulitplex: filename
    #remove .fastq.gz from filename
    m_dict = dict(zip(df_m.multiplex, df_m.file_name.str.replace('.fastq.gz','')))

    return(m_dict,s_dict)

#Read sample/multiplex manifests, create dicts
df_multiplex = pd.read_csv(multiplex_manifest,sep="\t")
df_samples = pd.read_csv(sample_manifest,sep="\t")
(multiplex_dict,samp_dict) = CreateSampleDicts(df_multiplex,df_samples)

#for each multiplexed sample: fastq file
for k,v in multiplex_dict.items():
    #create expected bc list
    bc_exp=[]
    for k2,v2 in samp_dict[k].items():
        bc_exp.append(v2.replace('N',''))

    #check the number of mismatches requested is possible given the expected barcodes
    for a, b in itertools.combinations(bc_exp, 2):
      
        #zip two barcodes and compare each letter
        compare=zip(a,b)
        diff_list=[]
        for i,j in compare:
            #if the letters are a mismatch, add to list
            if i!=j:
                diff_list.append(j)
            #if the mismatch list is not > than the mismatch value requested, print error and fail
        if(len(diff_list)<1 + int(mismatch)):
            print('The number of differences ({}) between barcodes {} and {} is less than or equal to the number of mismatches requested ({})'. format(len(diff_list),a,b,mismatch))
            
            sys.exit('Barcode strategy requires differences between barcodes is greater than mismatch allowance')

    #generate list of possible barcodes variations
    bc_mutants={}
    nuc_list = ['A','T','C','G','N']

    # for each expected barcode
    for bc in bc_exp:
        #vary the string by 1bp, 
        #push str to dict with expected value as the dict key
        for i in range(0,len(list(bc))):
            for nuc in nuc_list:
                tmp = list(bc)
                tmp[i]=nuc
                bc_mutants["".join(tmp)]=bc

                #if more than 1 mismatch is allowed, repeat variation
                if(int(mismatch)==2):
                    for j in range(0,len(list(bc))):
                        for nuc in nuc_list:
                            tmp2 = list("".join(tmp))
                            tmp2[j]=nuc
                            bc_mutants["".join(tmp2)]=bc
    
    # read list of counts:barcodes
    fastq_df = pd.read_csv(barcode_input,header=None)
    fastq_df[0] = fastq_df[0].str.lstrip() #remove leading spaces
    fastq_df[['counts','bc']] = fastq_df[0].str.split(expand=True)#split cols
    fastq_df["counts"] = pd.to_numeric(fastq_df["counts"])
    fastq_df=fastq_df.drop([0], axis=1) #remove merged col

    #for each barcode counted
    bc_dict = {}
    for index, row in fastq_df.iterrows():
        barcode=row['bc']

        #if barcode matches mutant, change barcode to expected value
        if barcode in bc_mutants:
          barcode_update = barcode
          barcode = bc_mutants[barcode_update]

        #if barcode isn't in dict, add it
        if k not in bc_dict: 
            bc_dict[k]={} 

        #create barcode:sample:count 
        #add barcode in dict count + row count from text file
        bc_dict[k][barcode] = bc_dict[k].get(barcode, 0) + row['counts']

    #select top 10 barcodes for each sample
    top_dict = dict(Counter(bc_dict[k]).most_common(10))
    
    #compare observed bc list with expected bc list, write output to text
    bc_obs=[]
    for k2,v2 in top_dict.items():
        bc_obs.append(k2)

    # if the expected list matches with top 10, print data
    # otherwise print error messages    
    check =  all(item in bc_obs for item in bc_exp)
    if check is True:
        file_save = output_dir + '/' + k + '_barcode.txt'
        f = open(file_save,"w+")
        f.write("\n* SampleID {}\n".format(k))
        f.write("\t + Number of mismatches allowed {}\n".format(mismatch))
        f.write("\t + The top barcodes identified {} include the expected barcodes {}\n\n".format(bc_obs, bc_exp))
        f.write("\t + List of top barcodes:counts \n")
        f.write("\t\t + " + json.dumps(top_dict))
        f.close()    
    else :
        file_save = output_dir + k + '_barcode_errors.txt'
        f = open(file_save,"w+")
        f.write("The top barcodes identified {} were not congruent with expected barcode list {}. Review associated img for more information.".format(bc_obs, bc_exp)) 
        f.close()

    # print barplot for top barcodes
    file_save = output_dir + '/' + k + '_barcode.png'
    plt.bar(*zip(*top_dict.items()))
    plt.suptitle('Top Barcodes: ' + k + '\n Number of mismatches allowed: ' + str(mismatch))
    plt.xticks(rotation='45', fontsize=6)
    for k,v in top_dict.items():
        plt.text(x=k , y =v+1 , s=str(v), color = 'black', fontweight='bold')
    plt.tight_layout()
    plt.savefig(file_save)