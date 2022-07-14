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
# parser$add_argument("-pu","--peak_unique", dest="peak_unique", required=TRUE, help="path for unique count peaks")
# parser$add_argument("-pa","--peak_all", dest="peak_all", required=TRUE, help="path MM count peaks")
# parser$add_argument("-ptot","--peak_total", dest="peak_total", required=TRUE, help="path total count peaks")
# parser$add_argument("-jj","--join_junction", dest="join_junction", required=TRUE, help="include junctions (T or F)")
parser$add_argument("-a","--anno_anchor", dest="anno_anchor", required=TRUE, help="5prime or max. Annotation for Spliced peaks will be based on 5` most region or region with most reads")
# parser$add_argument("-ce","--condense_exon", dest="condense_exon", required=TRUE, help="condense exons (T or F)")
parser$add_argument("-rdepth","--read_depth", dest="read_depth", required=TRUE, help="read depth filtering parameter")
# parser$add_argument("-de","--demethod", dest="demethod", required=TRUE, help="DE method (MAnorm, DiffBind or none)")
parser$add_argument("-s","--sample_id", dest="sample_id", required=TRUE, help="sample id")
parser$add_argument("-sp","--ref_species", dest="ref_species", required=TRUE, help="reference species")
# parser$add_argument("-anno","--anno_dir", dest="anno_dir", required=TRUE, help="path for annotation dir")
# parser$add_argument("-st","--splice_table", dest="splice_table", required=TRUE, help="path for Splice Table")
# parser$add_argument("-reft","--reftable_path", dest="reftable_path", required=TRUE, help="path for reftable")
# parser$add_argument("-g","--gencode_path", dest="gencode_path", required=TRUE, help="path for gencode")
# parser$add_argument("-i","--intron_path", dest="intron_path", required=TRUE, help="path for intron")
# parser$add_argument("-rmsk","--rmsk_path", dest="rmsk_path", required=TRUE, help="path for rmsk")
parser$add_argument("-tmp","--tmp_dir", dest="tmp_dir", required=TRUE, help="path for tmp dir")
parser$add_argument("-od","--out_dir", dest="out_dir", required=TRUE, help="path for output dir")
parser$add_argument("-of","--out_file", dest="out_file", required=TRUE, help="Output File name")
# parser$add_argument("-odm","--out_dir_DEP", dest="out_dir_DEP", required=FALSE, help="path for manorm output dir")
# parser$add_argument("-o","--output_file_error", dest="output_file_error", required=FALSE, help="path for output error file")
# parser$add_argument("-ss","--subset", dest="subset", required=T, help="subset")

args <- parser$parse_args()
rscript = args$rscript
peak_type = args$peak_type
# peak_unique = args$peak_unique
# peak_all = args$peak_all
# peak_total = args$peak_total
# join_junction = as.logical(args$join_junction)
# condense_exon = as.logical(args$condense_exon)
anno_anchor=args$anno_anchor
read_depth = as.numeric(args$read_depth)
# DEmethod = args$demethod
sample_id = args$sample_id
ref_species = args$ref_species
# anno_dir = args$anno_dir
# splice_table=args$splice_table
# reftable_path = args$reftable_path
# gencode_path = args$gencode_path
# intron_path = args$intron_path
# rmsk_path = args$rmsk_path
tmp_dir = args$tmp_dir
out_dir = args$out_dir
out_file = args$out_file
# out_dir_DEP = args$out_dir_DEP
# output_file_error = args$output_file_error
# subset=as.numeric(args$subset)


test=F
if(test){
  rm(list=setdiff(ls(), "params"))
  
  wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput2/"
  setwd(wd)
  wd="."
  
  sample_id = "FLAG_Ro_noPFA"
  
  FtrCount_trimmed=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions.txt")
  transc_SS=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions_transcripts_SameStrand.txt")
  transc_OS=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions_transcripts_OppoStrand.txt")
  EI_SS=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions_IntronExon_SameStrand.txt")
  EI_OS=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions_IntronExon_OppoStrand.txt")
  RMSK_SS=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions_RMSK_SameStrand.txt")
  RMSK_OS=paste0(wd,"/04_annotation/02_peaks/",sample_id,"_ALLreadPeaks_AllRegions_RMSK_OppoStrand.txt")
  
  
  
  peak_type= "ALL"
  anno_dir = paste0(wd,"/04_annotation/01_project/")
  ref_path = "/Users/homanpj/Documents/Resources/ref/"
  tmp_dir= paste0(wd,"/04_annotation/tmp/")
  
  reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/config/annotation_config.txt"
  
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
  read_depth = 3
  ref_species="mm10"
  anno_strand="SameStrand"
  
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
  
}

FtrCount_trimmed=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions.txt")
transc_SS=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions_transcripts_SameStrand.txt")
transc_OS=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions_transcripts_OppoStrand.txt")
EI_SS=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions_IntronExon_SameStrand.txt")
EI_OS=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions_IntronExon_OppoStrand.txt")
RMSK_SS=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions_RMSK_SameStrand.txt")
RMSK_OS=paste0(out_dir,sample_id,"_ALLreadPeaks_AllRegions_RMSK_OppoStrand.txt")


#source R script with functions
source(rscript)

# #annotation paths
# gencode_transc_path = paste0(anno_dir,"ref_gencode.txt")
# lncra_path = paste0(anno_dir,"lncRNA_gencode.txt")
# YRNA_path = paste0(anno_dir, "yRNA_repeatmasker.bed")
# srpRNA_path = paste0(anno_dir, "srpRNA_repeatmasker.bed")
# SKRNA_path = paste0(anno_dir, "7SKRNA_repeatmasker.bed")
# scRNA_path = paste0(anno_dir, "scRNA_repeatmasker.bed")
# tRNA_path = paste0(anno_dir, "tRNA_repeatmasker.bed")
# sncRNA_path = paste0(anno_dir, "sncRNA_gencode.bed")
# rRNA_BK00964_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
# rRNA_rmsk_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
# tRNA_rmsk_path = paste0(anno_dir, "tRNA_repeatmasker.bed")

#set id for files
file_id = paste0(sample_id,"_")
# file_id = sample_id



########################################
### read in files
########################################
FtrCount_trimmed=fread(paste0(FtrCount_trimmed), 
           header=T, sep="\t",stringsAsFactors = F,data.table=F)

transc_SS=fread(paste0(transc_SS), 
                header=T, sep="\t",stringsAsFactors = F,data.table=F)
transc_OS=fread(paste0(transc_OS), 
                header=T, sep="\t",stringsAsFactors = F,data.table=F)

EI_SS=fread(paste0(EI_SS), 
            header=T, sep="\t",stringsAsFactors = F,data.table=F)
EI_OS=fread(paste0(EI_OS), 
            header=T, sep="\t",stringsAsFactors = F,data.table=F)

RMSK_SS=fread(paste0(RMSK_SS), 
              header=T, sep="\t",stringsAsFactors = F,data.table=F)
RMSK_OS=fread(paste0(RMSK_OS), 
              header=T, sep="\t",stringsAsFactors = F,data.table=F)
  
  
########################################
## Merge files
########################################
# clms=c('chr','start','end','ID.1')
# transc=merge(transc_SS,transc_OS[,!colnames(transc_OS)%in%clms],by='ID',suffixes=c('_Same','_Oppo'))
# EI=merge(EI_SS,EI_OS[,!colnames(EI_OS)%in%clms],by='ID',suffixes=c('_Same','_Oppo'))
# RMSK=merge(RMSK_SS,RMSK_OS[,!colnames(RMSK_OS)%in%clms],by='ID',suffixes=c('_Same','_Oppo'))
# 
# peak_in=merge(transc,EI[,!colnames(EI)%in%clms],by='ID',suffixes=c('_Same','_Oppo'))
#   
#   ### protein coding peaks with not Correctly assigned exon overlap --taken from 09_peak_annotation_functions.R , IE_calling
#   peak_in[is.na(peak_in[,paste0('Same_feature')]) &
#             peak_in[,paste0('Same_gene_type')]%in%'protein_coding',paste0('Same_feature')]='exon'
#   
#   peak_in[is.na(peak_in[,paste0('Oppo_feature')]) &
#           peak_in[,paste0('Oppo_gene_type')]%in%'protein_coding',paste0('Oppo_feature')]='exon'
# 
#   
# peak_in=merge(peak_in,RMSK[,!colnames(EI)%in%clms],by='ID',suffixes=c('_Same','_Oppo'))
#   

########################################
## Merge files
########################################

clms=c('chr','start','end','strand','ID.1')
peak_inSS=merge(transc_SS,EI_SS[,!colnames(EI_SS)%in%clms],by='ID')
peak_inSS=merge(peak_inSS,RMSK_SS[,!colnames(RMSK_SS)%in%clms],by='ID')

peak_inOS=merge(transc_OS,EI_OS[,!colnames(EI_OS)%in%clms],by='ID')
peak_inOS=merge(peak_inOS,RMSK_OS[,!colnames(RMSK_OS)%in%clms],by='ID')

# intersect(colnames(transc_SS),colnames(EI_SS))
# intersect(colnames(transc_OS),colnames(EI_OS))

  ### protein coding peaks with not Correctly assigned exon overlap --taken from 09_peak_annotation_functions.R , IE_calling
  peak_inSS[is.na(peak_inSS[,paste0('Same_feature')]) &
            peak_inSS[,paste0('Same_gene_type')]%in%'protein_coding',paste0('Same_feature')]='exon'

  peak_inOS[is.na(peak_inOS[,paste0('Oppo_feature')]) &
              peak_inOS[,paste0('Oppo_gene_type')]%in%'protein_coding',paste0('Oppo_feature')]='exon'

##########################################################################################
############### Assigning Clip peak attributes by strand
##########################################################################################
# Not all Peaks overlap with a single feature so peak assignments were assigned by priority:
# ncRNA > Protein coding : Exonic > repeats > Pseudogene > Antisense Feature > 
#Protein Coding : Intronic > lncRNA > no Feature  
# All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table
print("Add peak attributes")

#merge attributes
peak_attrib_Oppo = peak_attributes(peak_inOS,"Oppo_")
PeaksdataOut = merge(peak_attributes(peak_inSS,"Same_"),
                     peak_attrib_Oppo[,colnames(peak_attrib_Oppo)[!colnames(peak_attrib_Oppo) %in% 
                                                                    c("chr","start","end","strand",'ID.1' )]],
                     by='ID')


##########################################################################################
############### Assigning Clip peak attributes - merged strands
########################################################################################## 

system.time({
PeaksdataOut=Assign_CLIP_attributes(PeaksdataOut)
})

##########################################################################################
############### Merge peak annotations with peak junctions
##########################################################################################
Peaksdata2_anno = merge(FtrCount_trimmed,
                        PeaksdataOut[,!colnames(PeaksdataOut)%in%c('chr','start','end','strand')],
                        by='ID')


##########################################################################################
############### Add annotations to junctions
##########################################################################################  
JoinJunc=file.exists(paste0(out_dir,file_id,peak_type,"readPeaks_PGene_TBL2.txt"))

if (JoinJunc==T) {
  print("Add annotations to Junctions")
  
  PGene_TBL2=read.table(paste0(out_dir,file_id,peak_type,"readPeaks_PGene_TBL2.txt"), 
             header=T, sep="\t")
  Junc_peaks=read.table(paste0(out_dir,file_id,peak_type,"readPeaks_Junc_peaks.txt"), 
             header=T, sep="\t")
  
  system.time({
  Peaksdata2_anno=Join_Junction_Combine(Peaksdata2_anno,read_depth,Junc_peaks,PGene_TBL2,anno_anchor)
  })
  
}else{
  
  Peaksdata2_anno=Peaksdata2_anno
  # #### collapse exon
  # if (condense_exon==T) {
  #   print("Collapse exon")
  #   
  #   Peaksdata2_anno=Collapse_Exon(Peaksdata2_anno)
  # }  
} 

##########################################################################################
## correct lnLc annotation - keep lnLc annotation above to hepl track source of lnc annotation
##########################################################################################
Peaksdata2_anno=correct_lnLc(Peaksdata2_anno)

##########################################################################################
### Write Output
##########################################################################################
#write out for junction annotation 
# write.table(Peaksdata2_anno,paste0(out_dir,file_id,peak_type,'readPeaks_annotation_complete.txt'),sep = "\t")
write.table(Peaksdata2_anno,paste0(out_file),sep = "\t",row.names = F)
