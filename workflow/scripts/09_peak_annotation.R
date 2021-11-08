#library
library(tidyr)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(data.table)
library(argparse)

#set args
parser <- ArgumentParser()
parser$add_argument("-r","--rscript", dest="rscript", required=TRUE, help="rscript for functions")
parser$add_argument("-pt","--peak_type", dest="peak_type", required=TRUE, help="peak type (all, unique)")
parser$add_argument("-pu","--peak_unique", dest="peak_unique", required=TRUE, help="path for unique peaks")
parser$add_argument("-pa","--peak_all", dest="peak_all", required=TRUE, help="path all peaks")
parser$add_argument("-jj","--join_junction", dest="join_junction", required=TRUE, help="include junctions (T or F)")
parser$add_argument("-ce","--condense_exon", dest="condense_exon", required=TRUE, help="condense exons (T or F)")
parser$add_argument("-rdepth","--read_depth", dest="read_depth", required=TRUE, help="read depth filtering parameter")
parser$add_argument("-de","--demethod", dest="demethod", required=TRUE, help="DE method (MAnorm, DiffBind or none)")
parser$add_argument("-s","--sample_id", dest="sample_id", required=TRUE, help="sample id")
parser$add_argument("-sp","--ref_species", dest="ref_species", required=TRUE, help="reference species")
parser$add_argument("-anno","--anno_dir", dest="anno_dir", required=TRUE, help="path for annotation dir")
parser$add_argument("-reft","--reftable_path", dest="reftable_path", required=TRUE, help="path for reftable")
parser$add_argument("-g","--gencode_path", dest="gencode_path", required=TRUE, help="path for gencode")
parser$add_argument("-i","--intron_path", dest="intron_path", required=TRUE, help="path for intron")
parser$add_argument("-rmsk","--rmsk_path", dest="rmsk_path", required=TRUE, help="path for rmsk")
parser$add_argument("-tmp","--tmp_dir", dest="tmp_dir", required=TRUE, help="path for tmp dir")
parser$add_argument("-od","--out_dir", dest="out_dir", required=TRUE, help="path for output dir")
parser$add_argument("-odm","--out_dir_DEP", dest="out_dir_DEP", required=FALSE, help="path for manorm output dir")
parser$add_argument("-o","--output_file_error", dest="output_file_error", required=FALSE, help="path for output error file")

args <- parser$parse_args()
rscript = args$rscript
peak_type = args$peak_type
peak_unique = args$peak_unique
peak_all = args$peak_all
join_junction = as.logical(args$join_junction)
condense_exon = as.logical(args$condense_exon)
read_depth = as.numeric(args$read_depth)
DEmethod = args$demethod
sample_id = args$sample_id
ref_species = args$ref_species
anno_dir = args$anno_dir
reftable_path = args$reftable_path
gencode_path = args$gencode_path
intron_path = args$intron_path
rmsk_path = args$rmsk_path
tmp_dir = args$out_dir
out_dir = args$out_dir
out_dir_DEP = args$out_dir_DEP
output_file_error = args$output_file_error

#source R script with functions
source(rscript)

##testing
testing="N"
if(testing=="Y"){
  rm(list=setdiff(ls(), "params"))
  wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mESC_clip_2/"
  setwd(wd)
  wd="."
  
  peak_type= "ALL"
  peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mESC_clip_2/03_peaks/03_allreadpeaks/Control_Clip_1_uniqueCounts.txt"
  peak_all = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mESC_clip_2/03_peaks/03_allreadpeaks/Control_Clip_1_allFracMMCounts.txt" # peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/13_counts/allreadpeaks/WT_fCLIP_50nt_uniqueCounts.txt"
  reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/config/annotation_config.txt"
  
  #output
  out_dir = paste0(wd,"/04_annotation/02_peaks/")
  out_dir_DEP =paste0(wd,"/14_MAnorm/")
  
  #project annotation files
  anno_dir = paste0(wd,"/04_annotation/01_project/")
  ref_path = "/Users/homanpj/Documents/Resources/ref/"
  tmp_dir= paste0(wd,"/04_annotation/02_peaks/")
  pipline_dir="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/"
  
  #feature information
  join_junction = "TRUE"
  condense_exon="TRUE"
  read_depth = 3
  DEmethod = "DiffBind"
  ref_species="mm10"
  sample_id = "Control_Clip_1"
  output_file_error= paste0(wd,"/04_annotation/02_peaks/")
  
  
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
  condense_exon="TRUE"
  read_depth = 5
  DEmethod = "MANORM"
  sample_id = "KO-NOPFA"
  ref_species="hg38"
  anno_dir =paste0(out_base, "annotations_test/04_annotation/01_project/")
  reftable_path = paste0(git_base,"config/annotation_config.txt")
  gencode_path = paste0(pipe_base,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
  intron_path = paste0(pipe_base,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
  rmsk_path = paste0(pipe_base,"annotations/hg38/repeatmasker/rmsk_GRCh38.txt")
  tmp_dir = paste0(out_base, "annotations_test/04_annotation/02_peaks/")
  out_dir = paste0(out_base, "annotations_test/04_annotation/02_peaks/")
  out_dir_DEP = paste0(out_base, "annotations_test/05_demethod/01_input/")
  output_file_error = paste0(out_base, "annotations_test/04_annotation/read_depth_error.txt")
}

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
##########################################################################################
############### unique, all read count input - then merge
##########################################################################################

FtrCount=merge(subset(FtrCount_input(peak_unique),select=-c(Geneid)),
               FtrCount_input(peak_all)[,c('ID','Counts')],
               by='ID',suffixes=c("_Unique","_fracMM"),
               all.x=T)
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
    FtrCount_out=Join_Junction(FtrCount,FtrCount_fracJCount)
      FtrCount_out$FtrCount_trimmed=FtrCount_out$FtrCount_trimmed[duplicated(FtrCount_out$FtrCount_trimmed)==F & !FtrCount_out$FtrCount_trimmed$chr%in%c('chrM'),]
      FtrCount_trimmed = FtrCount_out$FtrCount_trimmed
############################################################
  } else { 
    print("No junctions were identified")
    JoinJunc=F
    FtrCount_trimmed=FtrCount[FtrCount$Counts_fracMM>=read_depth,]
  }
  } else{
    JoinJunc=F
    FtrCount_trimmed=FtrCount[FtrCount$Counts_fracMM>=read_depth,]
  }
)
#remove duplicates, remove chrom M
FtrCount_trimmed=FtrCount_trimmed[duplicated(FtrCount_trimmed)==F & !FtrCount_trimmed$chr%in%c('chrM'),]

##########################################################################################
############### DEP
##########################################################################################
if (DEmethod=='NONE') {
  
  print("Differential Peaks skipped")

} else{
  
  DEP_input(FtrCount_trimmed)
} 

##########################################################################################
############### reference, annotation table
##########################################################################################

gencode_Anno_RNA_comb=Ref_Prep('gencode')
Anno_RNA_comb=Ref_Prep('Combined')


##########################################################################################
############### call peaks
##########################################################################################
print("Call peaks")

#prep peaks
peaks_Oppo=FtrCount_trimmed[,c('chr','start','end','ID','ID','strand')] 
peaks_Oppo$strand=gsub("\\+","pos",peaks_Oppo$strand) 
peaks_Oppo$strand=gsub("\\-","+",peaks_Oppo$strand)
peaks_Oppo$strand=gsub("pos","-",peaks_Oppo$strand)

#run functions
peaks_SameAnno = peak_calling(FtrCount_trimmed[,c('chr','start','end','ID','ID','strand')],"Same_")
peaks_OppoAnno = peak_calling(peaks_Oppo,"Oppo_")

##########################################################################################
############### INTRON EXON ANNOTATION
##########################################################################################
### Identify if CLIP peak overlaps with Intron or Exonic region   

#   Using GTF file from GENCODE v23 -mm10
#   Using GTF file from GENCODE v32 -hg38

# Peaks were annotated by whether they overlap with Host gene intron/exon region
# Intron coordinates were calculated from GTF file.
# 
# A second column was added to idenify if the peak also overlapped with the 5'UTR 3'UTR or CDS (Column: Featrue 2)
#########################################################################################3
print("Intron/exon annotation")

intron_exon=IntronExon_prep(gencode_path,intron_path,gencode_transc_path)

##########################################################################################
############### Intron Exon Calling
##########################################################################################
print("IE Calling")

intronexon_Same = IE_calling(peaks_SameAnno,'Same',"Same_")
intronexon_Opposite = IE_calling(peaks_OppoAnno,'Oppo',"Oppo_")

##########################################################################################
############### IDENTIFY PEAKS IN REPEAT REGIONS   
##########################################################################################
# Annotate all repeat regions/Classes identified in Repeatmasker Annotation file (UCSC Table browser)  
# Data was not filtered based on any of the identified Repeats.  
#  1) LINE/SINE   
#  2) LTR   
#  3) DNA   
#  4) Satalites   
#  5) Simple Repeats   
#  6) Low Complexity   
#  7) Other   
#  8) Unknown 
print("RMSK")

#read in rmsk
rmsk_GRCm38=rmsk_prep(rmsk_path)
  
rpmsk_Same = rpmsk_calling(intronexon_Same,"Same_")
rpmsk_Opposite = rpmsk_calling(intronexon_Opposite,"Oppo_")

##########################################################################################
############### Assigning Clip peak attributes by strand
##########################################################################################
# Not all Peaks overlap with a single feature so peak assignments were assigned by priority:
# ncRNA > Protein coding : Exonic > repeats > Pseudogene > Antisense Feature > 
#Protein Coding : Intronic > lncRNA > no Feature  
# All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table
print("Add peak attributes")

#merge attributes
peak_attrib_Oppo = peak_attributes(rpmsk_Opposite,"Oppo_")
PeaksdataOut = merge(peak_attributes(rpmsk_Same,"Same_"),
                     peak_attrib_Oppo[,colnames(peak_attrib_Oppo)[!colnames(peak_attrib_Oppo) %in% 
                                                                    c("chr","start","end","strand" )]],
                     by='ID')


##########################################################################################
############### Assigning Clip peak attributes - merged strands
########################################################################################## 

PeaksdataOut=Assign_CLIP_attributes(PeaksdataOut)
  

##########################################################################################
############### Merge peak annotations with peak junctions
##########################################################################################
Peaksdata2_anno = merge(FtrCount_trimmed,
                        PeaksdataOut[,!colnames(PeaksdataOut)%in%c('chr','start','end','strand')],
                        by='ID')


##########################################################################################
############### Add annotations to junctions
##########################################################################################  
if (JoinJunc==T) {
  print("Add annotations to Junctions")
  
  Peaksdata2_anno=Join_Junction_Combine(Peaksdata2_anno,read_depth,FtrCount_out)
  
}else{
  
  #### collapse exon
  if (condense_exon==T) {
    print("Collapse exon")
    
Peaksdata2_anno=Collapse_Exon(Peaksdata2_anno)
  }  
} 


##########################################################################################
## correct lnLc annotation - keep lnLc annotation above to hepl track source of lnc annotation
##########################################################################################
Peaksdata2_anno=correct_lnLc(Peaksdata2_anno)

##########################################################################################
### Write Output
##########################################################################################
#write out for junction annotation 
write.table(Peaksdata2_anno,paste0(out_dir,file_id,'peakannotation_complete.txt'),sep = "\t")

