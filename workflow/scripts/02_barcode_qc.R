# Purpose:
  # Current script is set up for two barcoding strategies; if additional are required they must be added below. 
  # Current accepted formats include:
  #  - NNNNNCACTGTNNNN = total len 15, seq len 6, position [6:11]
  #  - NNNTGGCNN = total len 9, seq len 4, position [4:7]
# Input
  # to create input, run counts
  # gunzip -c  /data/CCBR_Pipeliner/iCLIP/test/fastq/test_1.fastq.gz | 
  # sed -n '/@/{n;p}' | awk '{print substr($0, 4, 4);}' | sort -n | uniq -c > test.barcode_counts.txt

library(ggplot2)
library(argparse)
library(DescTools)

#set args
parser <- ArgumentParser()
parser$add_argument("-s","--sample_manifest", dest="sample_manifest", required=TRUE, help="/path/to/samples.tsv")
parser$add_argument("-ba","--multiplex_manifest", dest="multiplex_manifest", required=TRUE, help="/path/to/multiplex.tsv")
parser$add_argument("-bu","--barcode_input", dest="barcode_input", required=TRUE, help="/path/to/_barcode_counts.txt")
parser$add_argument("-o","--output_dir", dest="output_dir", required=TRUE, help="dir for barcode plots and summary")
parser$add_argument("-m","--mismatch", dest="mismatch", required=TRUE, help="number of mismatches allowed")
parser$add_argument("-i","--mpid", dest="mpid", required=TRUE, help="sampleid")
parser$add_argument("-q","--qc_dir", dest="qc_dir", required=TRUE, help="dir for failed QC manifests")

args <- parser$parse_args()
sample_manifest = args$sample_manifest
multiplex_manifest = args$multiplex_manifest
barcode_input = args$barcode_input
output_dir = args$output_dir
mismatch = as.integer(args$mismatch)
mpid = args$mpid
qc_dir = args$qc_dir

#test input
testing="N"
if(testing=="Y"){
  sample_manifest  = "~/../../Volumes/RBL_NCI/Wolin/mov10_par_Y_r2_01062023/manifests/samples.tsv"
  multiplex_manifest = "~/../../Volumes/RBL_NCI/Wolin/mov10_par_Y_r2_01062023/manifests/multiplex.tsv"
  barcode_input = "~/../../Volumes/RBL_NCI/Wolin/mov10_par_Y_r2_01062023/sam_testing/barcode_counts.txt"
  output_dir = "~/../../Volumes/RBL_NCI/Wolin/mov10_par_Y_r2_01062023/sam_testing"
  mismatch = 1
  mpid="CLIP"
}

#Read sample/multiplex manifests
df_multiplex = read.table(multiplex_manifest,header=TRUE,sep=",")
df_samples = read.table(sample_manifest,header=TRUE,sep=",")
df_counts = read.table(barcode_input,header=FALSE)
colnames(df_counts)=c("counts","barcode")

#create list of expected barcodes
bc_exp=gsub("N","",subset(df_samples,multiplex==mpid)$barcode)

# determine the length of the barcodes expected and set number of barcodes observed
if (length(bc_exp)<10){
  number_of_barcodes_observed=10
} else if (length(bc_exp)<20) {
  number_of_barcodes_observed=20
} else{
  number_of_barcodes_observed=as.numeric(length(bc_exp))+10
}

#Determine if the variance between any two barcodes are above the threshold of mismatches
string.diff.ex<-function(a="ATTCGAN",b="attTGTT",exclude=c("n","N","?"),ignore.case=TRUE){
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case==TRUE){
    a<-toupper(a)
    b<-toupper(b)
  }
  diff.a<-unlist(strsplit(a,split=""))
  diff.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(diff.a,diff.b)
  for(ex.loop in 1:length(exclude)){
    diff.d<-diff.d[,!(diff.d[1,]==exclude[ex.loop]|diff.d[2,]==exclude[ex.loop])]
  }
  differences<-sum(diff.d[1,]!=diff.d[2,])
  return(differences)
}

#determine all expected barcode comparisons
barcode_combo = as.data.frame(CombSet(bc_exp,2,repl=FALSE))
error_message = ""

#for each comparison determine the number of strings that are diff - 
#if less than the mismatch threshold then document error
for (rowid in rownames(barcode_combo)){
  val = string.diff.ex(barcode_combo[rowid,1],barcode_combo[rowid,2])
  if(val<mismatch+1){
    error_message = paste0(error_message, barcode_combo[rowid,1]," to ", barcode_combo[rowid,2], " has a difference of ",val, " \n")
  }
}

#if there were errors, print message and exit
if(error_message!=""){
  cat(paste0("The number of differences between barocdes is less than, or equal to, the number of mismatches in the following barcodes:\n",error_message))
} else{
  print("Barcode schema passes QC")
}

#Create mutant lists for each 
nuc_list=c("A","C","T","G","N")
for (bc in bc_exp){
  mut_list=""
  
  #for each position of string
  for(bc_pos in 1:length(strsplit(bc,split="")[[1]])){
    
    #vary by nucleotide list
    for(nuc in nuc_list){
      bc_mut = bc
      substr(bc_mut, bc_pos, bc_pos) <- nuc
      
      #save variation to mut_list
      mut_list = c(mut_list,bc_mut)
    }
  }
  #add expected bc
  mut_list=c(mut_list,bc)
  
  #assign list to bc variable name
  assign(bc,mut_list[!duplicated(mut_list)])
}

#Read through input and assign barcodeID
df_counts$barcode_id=df_counts$barcode
for (rowid in 1:nrow(df_counts)){
  #assign barcode observed
  bc_obs = df_counts[rowid,"barcode"]
  
  #for each of the barcodes_expected
  for(bc in bc_exp){

    #if the expected barcode is in that list, assign ID as expected
    #otherwise leave as new barcode
    if(bc_obs %in% get(bc)){
      df_counts[rowid,"barcode_id"]=paste0("*",bc)
    } else{
      next
    }
  }
}

#merge df
df_counts_merged = aggregate(counts~barcode_id, df_counts, sum)

#sort and subset
df_counts_merged = df_counts_merged[order(df_counts_merged$counts,decreasing = TRUE),][1:number_of_barcodes_observed,]

#compare observed bc list with expected bc list, write output to text
diff_lis=setdiff(paste0("*",bc_exp),df_counts_merged$barcode_id)

# if the expected list matches with top number_of_barcodes_observed, print text, otherwise print error messages
# if the expected list matches with top number_of_barcodes_observed, print plot
if(length(diff_lis)==0){
  print("Barcodes observed pass QC")
  #print text file
  file_save = paste0(output_dir,'/',mpid,'_barcode_passed.txt')
  fileConn<-file(file_save)
  line1 = paste0("\n* SampleID ",mpid, " \n")
  line2 = paste0("\t + Number of mismatches allowed ", mismatch, " \n")
  line3 = paste0("\t + The top barcodes identified (", paste(df_counts_merged$barcode_id,collapse = ' '),
                 ") include the expected barcodes (", paste(bc_exp,collapse = ' '), ")")
  line4 = "\t + List of top barcodes:counts \n"
  line5 = ""
  for(rowid in 1:nrow(df_counts_merged)){
    line5=paste0(line5,"\t\t + ", df_counts_merged[rowid,"barcode_id"],": ", df_counts_merged[rowid,"counts"], "\n")
  }
  
  writeLines(paste(line1,line2,line3,line4,line5,sep=""), fileConn)
  close(fileConn)
  
  #print plot
  file_save = paste0(output_dir,'/',mpid,'_barcode.png')
  p = ggplot(df_counts_merged, aes(x = reorder(barcode_id,-counts),counts)) +
    geom_bar(stat ="identity") +
    geom_text(aes(label = counts), vjust = -0.2, size=2)+
    labs(title=paste0("Barcode Counts for ", mpid), x ="Barcode IDs", y = "Barcode Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))
  ggsave(file_save,p)
  
} else{
  file_save = paste0(qc_dir,'/',mpid,'_barcode_errors.txt')
  fileConn<-file(file_save)
  line1 = paste0("\n* SampleID ",mpid, " \n")
  line2 = paste0("\t + Number of mismatches allowed ", mismatch, " \n")
  line3 = paste0("\t + The top barcodes identified (", paste(df_counts_merged$barcode_id,collapse = ' '),
                 ") were not congruent with the expected barcodes (", paste(bc_exp,collapse = ' '), ") \n")
  line4 = "\t + List of top barcodes:counts \n"
  
  line5 = ""
  for(rowid in 1:nrow(df_counts_merged)){
    line5=paste0(line5,"\t\t + ", df_counts_merged[rowid,"barcode_id"],": ", df_counts_merged[rowid,"counts"], "\n")
  }
  
  writeLines(paste(line1,line2,line3,line4,line5,sep=""), fileConn)
  close(fileConn)

  print(paste0("Barcodes failed QC check. Review file located here:",file_save,"."))
  print(paste0("To check the reads go to the /project/01_preprocess/01_fastq/ dir and run `zcat {name_of_sample} | wc -l. Divide this number by 4 to get the total number of reads"))
  print(paste0("If satisfied with this read count change the config/snakemake_config file barcode_qc_flag to IGNORE"))
}
