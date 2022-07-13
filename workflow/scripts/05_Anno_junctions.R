#library
suppressMessages(library(tidyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
suppressMessages(library(parallel))

#set args
parser <- ArgumentParser()
parser$add_argument("-r","--rscript", dest="rscript", required=TRUE, help="rscript for functions")
parser$add_argument("-pt","--peak_type", dest="peak_type", required=TRUE, help="peak type (all, unique)")
parser$add_argument("-pu","--peak_unique", dest="peak_unique", required=TRUE, help="path for unique count peaks")
parser$add_argument("-pa","--peak_all", dest="peak_all", required=TRUE, help="path MM count peaks")
parser$add_argument("-ptot","--peak_total", dest="peak_total", required=TRUE, help="path total count peaks")
parser$add_argument("-jj","--join_junction", dest="join_junction", required=TRUE, help="include junctions (T or F)")
parser$add_argument("-a","--anno_anchor", dest="anno_anchor", required=TRUE, help="5prime or max. Annotation for Spliced peaks will be based on 5` most region or region with most reads")
# parser$add_argument("-ce","--condense_exon", dest="condense_exon", required=TRUE, help="condense exons (T or F)")
parser$add_argument("-rdepth","--read_depth", dest="read_depth", required=TRUE, help="read depth filtering parameter")
parser$add_argument("-de","--demethod", dest="demethod", required=TRUE, help="DE method (MAnorm, DiffBind or none)")
parser$add_argument("-s","--sample_id", dest="sample_id", required=TRUE, help="sample id")
parser$add_argument("-sp","--ref_species", dest="ref_species", required=TRUE, help="reference species")
# parser$add_argument("-anno","--anno_dir", dest="anno_dir", required=TRUE, help="path for annotation dir")
parser$add_argument("-st","--splice_table", dest="splice_table", required=TRUE, help="path for Splice Table")
# parser$add_argument("-reft","--reftable_path", dest="reftable_path", required=TRUE, help="path for reftable")
# parser$add_argument("-g","--gencode_path", dest="gencode_path", required=TRUE, help="path for gencode")
# parser$add_argument("-i","--intron_path", dest="intron_path", required=TRUE, help="path for intron")
# parser$add_argument("-rmsk","--rmsk_path", dest="rmsk_path", required=TRUE, help="path for rmsk")
parser$add_argument("-tmp","--tmp_dir", dest="tmp_dir", required=TRUE, help="path for tmp dir")
parser$add_argument("-od","--out_dir", dest="out_dir", required=TRUE, help="path for output dir")
parser$add_argument("-of","--out_file", dest="out_file", required=TRUE, help="path for output dir")
parser$add_argument("-odm","--out_dir_DEP", dest="out_dir_DEP", required=FALSE, help="path for manorm output dir")
parser$add_argument("-o","--output_file_error", dest="output_file_error", required=FALSE, help="path for output error file")
# parser$add_argument("-ss","--subset", dest="subset", required=T, help="subset")

args <- parser$parse_args()
rscript = args$rscript
peak_type = args$peak_type
peak_unique = args$peak_unique
peak_all = args$peak_all
peak_total = args$peak_total
join_junction = as.logical(args$join_junction)
# condense_exon = as.logical(args$condense_exon)
anno_anchor=args$anno_anchor
read_depth = as.numeric(args$read_depth)
DEmethod = args$demethod
sample_id = args$sample_id
ref_species = args$ref_species
anno_dir = args$anno_dir
splice_table=args$splice_table
# reftable_path = args$reftable_path
# gencode_path = args$gencode_path
# intron_path = args$intron_path
# rmsk_path = args$rmsk_path
tmp_dir = args$tmp_dir
out_dir = args$out_dir
out_file = args$out_file
out_dir_DEP = args$out_dir_DEP
output_file_error = args$output_file_error
# subset=as.numeric(args$subset)

##testing
testing="N"
if(testing=="Y"){
  rm(list=setdiff(ls(), "params"))
  # wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/Sam/novoalign_v4/mm10/"
  wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput"
  setwd(wd)
  wd="."
  
  # peak_type= "ALL"
  # peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/Sam/novoalign_v4/mm10/08_counts/test1/03_allreadpeaks/Ro_Clip_1_test1_uniqueCounts.txt"
  # peak_all = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/Sam/novoalign_v4/mm10/08_counts/test1/03_allreadpeaks/Ro_Clip_1_test1_allFracMMCounts.txt" # peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/13_counts/allreadpeaks/WT_fCLIP_50nt_uniqueCounts.txt"
  # peak_total = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/Sam/novoalign_v4/mm10/08_counts/test1/03_allreadpeaks/Ro_Clip_1_test1_allFracMMCounts.txt" # peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/13_counts/allreadpeaks/WT_fCLIP_50nt_uniqueCounts.txt"
  # splice_table="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/Sam/novoalign_v4/mm10/08_counts/test1/03_allreadpeaks/tmp/KO_connected_peaks.txt"

  peak_type= "ALL"
  peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput/03_peaks/03_counts/FLAG_Ro_fclip_ALLreadpeaks_uniqueCounts.txt"
  peak_all = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput/03_peaks/03_counts/FLAG_Ro_fclip_ALLreadpeaks_FracMMCounts.txt" # peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/13_counts/allreadpeaks/WT_fCLIP_50nt_uniqueCounts.txt"
  peak_total = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput/03_peaks/03_counts/FLAG_Ro_fclip_ALLreadpeaks_totalCounts.txt" # peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/13_counts/allreadpeaks/WT_fCLIP_50nt_uniqueCounts.txt"
  splice_table="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput/04_annotation/02_peaks/FLAG_Ro_fclip_ALL_connected_peaks.txt"
  # 
  subset=5000
  
  reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/config/annotation_config.txt"
  
  # #output
  # out_dir = paste0(wd,"/09_annotations/test1/02_peaks/")
  # out_dir_DEP =paste0(wd,"/14_MAnorm/test1/")
  # 
  # #project annotation files
  # anno_dir = paste0(wd,"/09_annotations/test1/01_project/")
  # ref_path = "/Users/homanpj/Documents/Resources/ref/"
  # tmp_dir= paste0(wd,"/09_annotations/test1/02_peaks/")
  
  #output
  out_dir = paste0(wd,"/04_annotation/02_peaks/")
  out_dir_DEP =paste0(wd,"/05_MAnorm/")
  
  #project annotation files
  anno_dir = paste0(wd,"/04_annotation/01_project/")
  ref_path = "/Users/homanpj/Documents/Resources/ref/"
  tmp_dir= paste0(wd,"/04_annotation/02_peaks/")
  
  
  #feature information
  join_junction = "TRUE"
  anno_anchor="max_total"
  # condense_exon="TRUE"
  read_depth = 3
  DEmethod = "MAnorm"
  ref_species="mm10"
  sample_id = "FLAG_Ro_fclip"
  # sample_id = "Control_Clip_1"
  output_file_error= paste0(wd,"/04_annotation/testoutput/02_peaks/")
  
  # rscript="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/workflow/scripts/09_peak_annotation_functions.R"
  rscript="/Volumes/RBL_NCI/Wolin/Phil/test/testoutput2/09_peak_annotation_functions.R"
  
  if(ref_species == "mm10"){
    gencode_path = paste0(ref_path, "mm10/Gencode_VM23/fromGencode/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.txt")
    intron_path = paste0(ref_path, "mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
    rmsk_path = paste0(ref_path,"mm10/repeatmasker/rmsk_GRCm38.txt")
    
  } else if (ref_species == "hg38"){
    gencode_path = paste0(ref_path,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
    intron_path = paste0(ref_path,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
    rmsk_path = paste0(ref_path,"hg38/repeatmasker/rmsk_GRCh38.txt")
  }
} else if(testing=="SSC"){
  git_base="~/../../Volumes/sevillas2/git/iCLIP/"
  out_base="~/../../Volumes/sevillas2-1/"
  pipe_base="~/../../Volumes/CCBR_Pipeliner/iCLIP/ref/annotations/"
  rbl_base="~/../../Volumes/RBL_NCI-1/"
  
  peak_type= "ALL"
  peak_unique = paste0(rbl_base,"Wolin/6-22-21-HaCaT_fCLIP/12_counts/allreadpeaks/KO_fCLIP_uniqueCounts.txt")
  peak_all = paste0(rbl_base,"Wolin/6-22-21-HaCaT_fCLIP/12_counts/allreadpeaks/KO_fCLIP_allFracMMCounts.txt")
  join_junction = "TRUE"
  # condense_exon="TRUE"
  read_depth = 5
  DEmethod = "MANORM"
  sample_id = "KO-NOPFA"
  ref_species="hg38"
  anno_dir =paste0(out_base, "annotations_test/05_annotation/01_project/")
  reftable_path = paste0(git_base,"config/annotation_config.txt")
  gencode_path = paste0(pipe_base,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
  intron_path = paste0(pipe_base,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
  rmsk_path = paste0(pipe_base,"annotations/hg38/repeatmasker/rmsk_GRCh38.txt")
  tmp_dir = paste0(out_base, "annotations_test/05_annotation/02_peaks/")
  out_dir = paste0(out_base, "annotations_test/05_annotation/02_peaks/")
  out_dir_DEP = paste0(out_base, "annotations_test/06_MAnorm/01_input/")
  output_file_error = paste0(out_base, "annotations_test/05_annotation/read_depth_error.txt")
}

#source R script with functions
source(rscript)

#annotation paths
gencode_transc_path = paste0(anno_dir,"ref_gencode.txt")
lncra_path = paste0(anno_dir,"lncRNA_gencode.txt")
YRNA_path = paste0(anno_dir, "yRNA_repeatmasker.bed")
srpRNA_path = paste0(anno_dir, "srpRNA_repeatmasker.bed")
SKRNA_path = paste0(anno_dir, "7SKRNA_repeatmasker.bed")
scRNA_path = paste0(anno_dir, "scRNA_repeatmasker.bed")
tRNA_path = paste0(anno_dir, "tRNA_repeatmasker.bed")
sncRNA_path = paste0(anno_dir, "sncRNA_gencode.bed")
rRNA_BK00964_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
rRNA_rmsk_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
tRNA_rmsk_path = paste0(anno_dir, "tRNA_repeatmasker.bed")

#set id for files
file_id = paste0(sample_id,"_")
# file_id = sample_id
##########################################################################################
############### unique, all read count input - then merge
##########################################################################################

FtrCount=merge(base::subset(FtrCount_input(peak_unique),select=-c(Geneid)),
               FtrCount_input(peak_all)[,c('ID','Counts')],
               by='ID',suffixes=c("_Unique","_fracMM"),
               all.x=T)
FtrCount=merge(FtrCount,FtrCount_input(peak_total)[,c('ID','Counts')],by='ID',all.x=T)
FtrCount=dplyr::rename(FtrCount,'Counts_total'=Counts)

FtrCount = FtrCount %>%
  dplyr::rename(
    chr = Chr,
    start = Start,
    end = End,
    strand = Strand
  )
FtrCount=FtrCount[duplicated(FtrCount)==F,]


##########################################################################################
############### run filter check
##########################################################################################
if(nrow(FtrCount[FtrCount$Counts_fracMM>=read_depth,])==0){
  text1 = "Based on the read_depth parameters give all peaks are being filtered. Please update read_depth in config file to continue processing"
  writeLines(paste0(text1, "\n", 
                    "max read counts: ", max(FtrCount$Counts_fracMM), "\n",
                    "min read counts: ", min(FtrCount$Counts_fracMM), "\n",
                    "current read count threshold: ", read_depth), output_file_error)
  print("Read depth params error")
  invokeRestart("abort")
} else{
  print("Peaks pass read_depth filter")
}


#########################################################################################
############### splice junctions
##########################################################################################
# FtrCount=FtrCount[1:2000,]

system.time(
  if (join_junction) {  
    print("Running Join Junction")
    ############### Pre-Processing
    #grab peaks file
    if (peak_type=="UNIQUE") {
      FtrCount_fracJCount=read.delim(paste0(peak_unique,".jcounts"), 
                                     header=T,sep="\t",stringsAsFactors = F,comment.char = '#') 
    } else{
      FtrCount_fracJCount=read.delim(paste0(peak_all,".jcounts"), 
                                     header=T,sep="\t",stringsAsFactors = F,comment.char = '#') 
    }    
    
    #rename last Col
    colnames(FtrCount_fracJCount)[ncol(FtrCount_fracJCount)]='counts'
    
    #filter if primary gene is missing
    FtrCount_fracJCount=FtrCount_fracJCount[!is.na(FtrCount_fracJCount$PrimaryGene),]
    
    
    if (nrow(FtrCount_fracJCount)>0) {
      print("Junctions were identitified")
      JoinJunc=T    
      ############################################################    
      # FtrCount_outO=Join_Junction_original(FtrCount,FtrCount_fracJCount)
      FtrCount_out=Join_Junction(FtrCount,splice_table)
      FtrCount_out$FtrCount_trimmed=FtrCount_out$FtrCount_trimmed[duplicated(FtrCount_out$FtrCount_trimmed)==F & !FtrCount_out$FtrCount_trimmed$chr%in%c('chrM'),]
      FtrCount_trimmed = FtrCount_out$FtrCount_trimmed
      print(paste0('AnnotatePeaks: ',nrow(FtrCount_trimmed)))
      print(paste0('SplicePeaks: ',length(FtrCount_out$Junc_peaks)))
      
      write.table(FtrCount_out$PGene_TBL2,paste0(out_dir,file_id,peak_type,"readPeaks_PGene_TBL2.txt"),sep = "\t",row.names = F)
      write.table(FtrCount_out$Junc_peaks,paste0(out_dir,file_id,peak_type,"readPeaks_Junc_peaks.txt"),sep = "\t",row.names = F)
      
      ############################################################
    } else { 
      print("No junctions were identified")
      JoinJunc=F
      FtrCount_trimmed=FtrCount[FtrCount$Counts_fracMM>=read_depth,]
      FtrCount_trimmed$IDmerge=NA
    }
  } else{
    JoinJunc=F
    FtrCount_trimmed=FtrCount[FtrCount$Counts_fracMM>=read_depth,]
    FtrCount_trimmed$IDmerge=NA
  }
)
#remove duplicates, remove chrom M
FtrCount_trimmed=FtrCount_trimmed[duplicated(FtrCount_trimmed)==F & !FtrCount_trimmed$chr%in%c('chrM'),]


#########################################################################################
############### subset table
##########################################################################################
# if (subset>nrow(FtrCount_trimmed)) {
#   subset=nrow(FtrCount_trimmed)
# }else{subset=subset}
# 
# FtrCount_trimmed=FtrCount_trimmed[sample(nrow(FtrCount_trimmed),subset),]

print(paste0("Annotate Rows:",nrow(FtrCount_trimmed)))
##########################################################################################
############### DEP
##########################################################################################
if (DEmethod=='NONE') {
  
  print("Differential Peaks skipped")
  
} else{
  
  DEP_input(FtrCount_trimmed)
} 


write.table(FtrCount_trimmed,paste0(out_file),sep = "\t",row.names = F)










