import pandas as pd
import re
import sys
from datetime import datetime


def check_header(input_df,expected_list,check_file,error_log):
    check =  all(item in list(input_df.columns.values) for item in expected_list)

    if check is False:
        error_log.append("The file {} does not include all required elements:{}\n".format(check_file,set(expected_list).difference(set(list(input_df.columns.values)))))

    return(error_log)

def check_multiplex(input_df,select_file,error_log):
    for select_cols in input_df.columns.values:
        tmp_list = input_df[select_cols].tolist()

        #file_name
        if select_cols == "file_name":
            #check if all filenames are unique
            if len(tmp_list) != len(set(tmp_list)):
                error_log.append("File names must be unique: {}\n".format(set([x for x in tmp_list if tmp_list.count(x) > 1])))

        #check all filenames end in .fastq.gz
            for item in tmp_list:
                if not (item.endswith('.fastq.gz')):
                    error_log.append("All items in file_name column must end in fastq.gz - please update filename: {}\n".format(item))

        #multiplex
        if select_cols == "multiplex":
            #check values are alpha/numeric or _
            regex = re.compile(r'[A-Z]_[a-z][0-9]')
            for item in tmp_list:
                if(regex.search(item) != None):
                    error_log.append("{} values can only contain alpha/numeric values or _ - please review: {}\n".format(select_cols,item))
    return(error_log)

def check_samples(input_df,select_file,error_log):
    for select_cols in input_df.columns.values:
        tmp_list = input_df[select_cols].tolist()

        if select_cols == "sample":
            #check if col values are unique
            if len(tmp_list) != len(set(tmp_list)):
                error_log.append("{} names must be unique in {}: {}".format(select_cols,select_file,set([x for x in tmp_list if tmp_list.count(x) > 1])))
            #check values are alpha/numeric or _
            regex = re.compile(r'[A-Z]_[a-z][0-9]')
            for item in tmp_list:
                if(regex.search(item) != None):
                    error_log.append("{} values can only contain alpha/numeric values or _ - please review: {}".format(select_cols,item))

        if select_cols == "multiplex":
            for item in tmp_list :
                #check samples for unique barcodes
                sub_df = input_df[input_df[select_cols] == item]
                bc_list = sub_df["barcode"].tolist()

                if len(bc_list) != len(set(bc_list)):
                    error_log.append("Barcodes must be unique by sample in {}: sample {} contains duplicate barcodes {}".format(select_file,item,set([x for x in bc_list if bc_list.count(x) > 1])))

                #check values are alpha
                for item in bc_list:
                    if(item.isalpha() != True):
                        error_log.append("Barcode values can only contain alpha characters- please review: {}".format(item))

        #adaptor in samples.tsv
        if select_cols=="adaptor":
          #check values are alpha
          for item in tmp_list:
              if(item.isalpha() != True):
                  error_log.append("{} values can only contain alpha characters- please review: {}".format(select_cols,item))

    return(error_log)

def check_contrasts(input_df,select_file,error_log):
    #check if there are any NA values, otherwise check alpha/num
    check_na = input_df.isnull().values.any()

    if check_na == False:
      for select_cols in input_df.columns.values:
          tmp_list = input_df[select_cols].tolist()
          #check values are alpha/numeric or _
          regex = re.compile(r'[A-Z]_[a-z][0-9]')
          for item in tmp_list:
            if(regex.search(item) != None):
              error_log.append("{} values can only contain alpha/numeric values or _ - please review: {}".format(select_cols,item))
    else:
      error_log.append("Contrast manifest can only contain two comparisons per line. Please review")
    return(error_log)

def check_multiplex_samples(input_df_m,input_df_s,error_log):
    l1 = list(input_df_m.multiplex.unique())
    l2 = list(input_df_s.multiplex.unique())

    l_1to2 = list(set(l1)-set(l2))
    l_2to1 = list(set(l2)-set(l1))

    if len(l_1to2) != 0 or len(l_2to1) != 0:
        error_log.append("Multiplex ID's must be consistent between both files.")

    if len(l_1to2) !=0:
        error_log.append("The following was only found in the multiplex tsv: {}".format(l_1to2))

    if len(l_2to1) !=0:
        error_log.append("The following was only found in the sample tsv: {}".format(l_2to1))

    return(error_log)

def check_MANORM_contrast(input_df_s,input_df_c,error_log,de_log):
    #create list of samples
    l1 = list(input_df_s["sample"].unique())

    #create list of requested sample contrasts
    l2 = list(input_df_c["sample"].unique())

    #create list of requested background contrasts      
    l3 = list(input_df_c["background"].unique())

    #make sure all samples/background are samples found in workflow
    l_2to1 = list(set(l2)-set(l1))
    l_3to1 = list(set(l3)-set(l1))

    #if any samples are missing, print error
    if len(l_2to1) !=0:
        error_log.append("The following sample(s) was/were not found in the sample_manifest tsv but were found in the contrasts manifest: {}. Please review manifests.".format(l_2to1))
      
    #if any background samples are missing, print error
    if len(l_3to1) !=0:
        error_log.append("The following background sample(s) was/were not found in the sample_manifest tsv but were found in the contrasts manifest{}. Please review manifests.".format(l_3to1))
      
    #if all samples and background samples are present, print confirmation of MANORM
    if len(l_2to1) ==0 & len(l_3to1) ==0:
        de_log.append("The following sample to sample comparisons will be run with MANORM:\n{}.".format(input_df_c))

    return(error_log,de_log)

def check_DIFFBIND_contrast(input_df_s,input_df_c,error_log,de_log):
    #create list of groups
    l1 = list(input_df_s["group"].unique())

    #create list of requested group contrasts
    l2 = list(input_df_c["group"].unique())

    #create list of requested background group contrasts      
    l3 = list(input_df_c["background"].unique())

    #make sure all sample groups / background groups are groups found in workflow
    l_2to1 = list(set(l2)-set(l1))
    l_3to1 = list(set(l3)-set(l1))

    #if any sample groups are missing, print error
    if len(l_2to1) !=0:
        error_log.append("The following sample group(s) was/were not found in the sample_manifest tsv but were found in the contrasts manifest: {}. Please review manifests.".format(l_2to1))

    #if any background groups are missing, print error
    if len(l_3to1) !=0:
        error_log.append("The following background group(s) was/were not found in the sample_manifest tsv but were found in the contrasts manifest{}. Please review manifests.".format(l_3to1))

    #if all samples and background samples are present, print confirmation of MANORM
    if len(l_2to1) ==0 & len(l_3to1) ==0:
        de_log.append("The following group to group comparisons will be run with DIFFBIND:")

        #for each row in contrast manifest    
        for index, row in input_df_c.iterrows():
            #create list of samples and background in each group
            l1 = list(input_df_s[(input_df_s['group']==row['group'])]['sample'].unique())
            l2 = list(input_df_s[(input_df_s['group']==row['background'])]['sample'].unique())
        
        #output info to de_log
        de_log.append("- {} to {}".format(row['group'],row['background']))
        de_log.append("-- this includes groups: {}".format(l1))
        de_log.append("-- this includes background: {}".format(l2))
    
    return(error_log,de_log)

#create log
error_log = []
de_log = []

#Check multiplex file
check_file = sys.argv[2]
m_df = pd.read_csv(check_file,sep=",")
m_req = ['file_name','multiplex']
error_log = check_header(m_df,m_req,check_file,error_log)
error_log = check_multiplex(m_df,check_file,error_log)

#Check samples file
check_file = sys.argv[3]
s_df = pd.read_csv(check_file,sep=",")
s_req = ['multiplex','barcode','sample','adaptor']
error_log = check_header(s_df,s_req,check_file,error_log)
error_log = check_samples(s_df,check_file,error_log)

#Check concordance between multiplex and sample files
error_log = check_multiplex_samples(m_df,s_df,error_log)

#check contrast file, if necessary
DE_method = sys.argv[4]
if DE_method == "MANORM":
    #Check contrast file
    check_file = sys.argv[5]
    c_df = pd.read_csv(check_file,sep=",")
    c_req = ['sample','background']
    error_log = check_header(c_df,c_req,check_file,error_log)
    error_log = check_contrasts(c_df,check_file,error_log)
      
    #Check concordance between sample and contrast files
    error_log,de_log = check_MANORM_contrast(s_df,c_df,error_log,de_log)

elif DE_method == "DIFFBIND":
    #Check contrast file
    check_file = sys.argv[5]
    c_df = pd.read_csv(check_file,sep=",")
    c_req = ['group','background']
    error_log = check_header(c_df,c_req,check_file,error_log)
    error_log = check_contrasts(c_df,check_file,error_log)
      
    #Check concordance between sample and contrast files
    error_log,de_log = check_DIFFBIND_contrast(s_df,c_df,error_log, de_log)

# if there is no error, create empty file
# if there is an error, create error log
if len(error_log)==0:
    #create empty no_error file
    new_path = str(sys.argv[1]) + "no_errors.txt"
    open(new_path,"w+").write('\n'.join(error_log))
else:
    date_str = '%s%s_%s' % (sys.argv[1],'contains_errors',datetime.now().strftime('%Y%m%d'))
    new_path =  date_str + '.txt'
    open(new_path, 'w+').write('\n'.join(error_log))

#print out DE log
if len(de_log)!=0:
    new_path = str(sys.argv[1]) + DE_method + "_comparisons.txt"
    open(new_path,"w+").write('\n'.join(de_log))
