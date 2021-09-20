#library
library(tidyr)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(data.table)
library(argparse)

#set args
parser <- ArgumentParser()
parser$add_argument("-s","--samplename", dest="samplename", required=TRUE, help="Sampleid for sample DiffBind comparison")
parser$add_argument("-b","--background", dest="background", required=TRUE, help="Sampleid for background DiffBind comparison")
parser$add_argument("-pt","--peak_type", dest="peak_type", required=TRUE, help="peak type (all, unique)")
parser$add_argument("-jj","--join_junction", dest="join_junction", required=TRUE, help="include junctions (T or F)")
parser$add_argument("-m","--samplemanifest", dest="samplemanifest", required=TRUE, help="path to sample_manifest.tsv")
parser$add_argument("-pos","--pos_DB", dest="pos_DB", required=TRUE, help="Output from DiffBind for Positive peaks")
parser$add_argument("-neg","--neg_dB", dest="neg_DB", required=TRUE, help="Output from DiffBind for Negative peaks")
parser$add_argument("-Sanno","--anno_dir_sample", dest="anno_dir_sample", required=TRUE, help="path for dir with annotations for each sample (peakannotation_final_table.txt)")
parser$add_argument("-reft","--reftable_path", dest="reftable_path", required=TRUE, help="path for reftable")
parser$add_argument("-anno","--anno_dir_project", dest="anno_dir_project", required=TRUE, help="path for annotation dir")
parser$add_argument("-sp","--ref_species", dest="ref_species", required=TRUE, help="reference species")
parser$add_argument("-g","--gencode_path", dest="gencode_path", required=TRUE, help="path for gencode")
parser$add_argument("-i","--intron_path", dest="intron_path", required=TRUE, help="path for intron")
parser$add_argument("-rmsk","--rmsk_path", dest="rmsk_path", required=TRUE, help="path for rmsk")
parser$add_argument("-fs","--function_script", dest="function_script", required=TRUE, help="path for 09_annotations_function script")
parser$add_argument("-od","--out_dir", dest="out_dir", required=TRUE, help="path for output dir")
parser$add_argument("-tmp","--tmp_dir", dest="tmp_dir", required=TRUE, help="path for tmp dir")
parser$add_argument("-o","--output_file_error", dest="output_file_error", required=FALSE, help="path for output error file")

args <- parser$parse_args()

samplename = args$samplename
background = args$background
peak_type = args$peak_type
join_junction = as.logical(args$join_junction)

samplemanifest=args$samplemanifest
pos_DB = args$pos_DB
neg_DB = args$neg_DB
anno_dir_sample = args$anno_dir_sample

reftable_path = args$reftable_path
anno_dir_project = args$anno_dir_project

ref_species = args$ref_species
gencode_path = args$gencode_path
intron_path = args$intron_path
rmsk_path = args$rmsk_path

function_script = args$function_script
out_dir = args$out_dir
tmp_dir = args$tmp_dir
output_file_error = args$output_file_error

##testing
testing="Y"
if(testing=="Y"){
  rm(list=setdiff(ls(), "params"))
  wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/2-2-21_fCLIP_duplicates/"
  setwd(wd)
  wd="."
  
  samplemanifest=paste0(wd,'/sample_manifest.tsv')
  samplename='WT' # input from contrasts table
  background='KO'
  peak_type= "ALL"
  
  pos_DB = paste0(wd,"/06_DEP/02_analysis/WTvsKO_DiffBind/WTvsKO_DiffBind_P.txt")
  neg_DB = paste0(wd,"/06_DEP/02_analysis/WTvsKO_DiffBind/WTvsKO_DiffBind_N.txt")
  anno_dir_sample=paste0(wd,'/05_annotation/')
  
  reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/config/annotation_config.txt"
  anno_dir_project = paste0(wd,"/05_annotation/01_project/")
  out_dir = paste0(wd,"/06_DEP/02_analysis/WTvsKO_DiffBind/")
  tmp_dir= paste0(wd,"/06_DEP/02_analysis/WTvsKO_DiffBind/")
  
   ref_path = "/Users/homanpj/Documents/Resources/ref/"
   ref_species="hg38"
   join_junction=T
   
   function_script="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/"

  
  if(ref_species == "mm10"){
    gencode_path = paste0(ref_path, "mm10/Gencode_VM23/fromGencode/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.txt")
    intron_path = paste0(ref_path, "mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
    rmsk_path = paste0(ref_path,"mm10/repeatmasker/rmsk_GRCm38.txt")
    
  } else if (ref_species == "hg38"){
    gencode_path = paste0(ref_path,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
    intron_path = paste0(ref_path,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
    rmsk_path = paste0(ref_path,"hg38/repeatmasker/rmsk_GRCh38.txt")
  }
} else if (testing=="SSC"){
  #not a param
  ref_path = "~/../../Volumes/CCBR_Pipeliner-1/iCLIP/ref/annotations/"
  data_dir="~/../../Volumes/data"
  
  samplename='WT'
  background='KO'
  peak_type= "ALL"
  samplemanifest=paste0(data_dir,"/diffbind/sample_manifest.tsv")
  pos_DB =paste0(data_dir,"/diffbind/05_demethod/02_analysis/WT_vs_KO/WT_vs_KO_DIFFBINDTable_N.txt")
  neg_DB =paste0(data_dir,"/diffbind/05_demethod/02_analysis/WT_vs_KO/WT_vs_KO_DIFFBINDTable_P.txt")
  anno_dir_sample=paste0(data_dir,"/diffbind/05_demethod/04_annotation/")
  reftable_path =paste0(data_dir,"/diffbind/config/annotation_config.txt")
  anno_dir_project = paste0(data_dir,"/04_annotation/01_project/")
  function_script="~/sevillas2-1/git/iCLIP/workflow/scripts/09_peak_annotation_functions.R"
  out_dir = paste0(data_dir,"/05_demethod/02_analysis/WTvsKO/")
  tmp_dir= paste0(data_dir,"/05_demethod/02_analysis/WTvsKO/")
  
  ref_species="hg38"
  join_junction=T
  
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

#source annotations_function.R script
source(function_script)

#annotation paths
gencode_transc_path = paste0(anno_dir_project,"ref_gencode.txt")
lncra_path = paste0(anno_dir_project,"lncRNA_gencode.txt")
YRNA_path = paste0(anno_dir_project, "yRNA_repeatmasker.bed")
srpRNA_path = paste0(anno_dir_project, "srpRNA_repeatmasker.bed")
SKRNA_path = paste0(anno_dir_project, "7SKRNA_repeatmasker.bed")
scRNA_path = paste0(anno_dir_project, "scRNA_repeatmasker.bed")
tRNA_path = paste0(anno_dir_project, "tRNA_repeatmasker.bed")
sncRNA_path = paste0(anno_dir_project, "sncRNA_gencode.bed")
rRNA_BK00964_path = paste0(anno_dir_project, "rRNA_repeatmasker.bed")
rRNA_rmsk_path = paste0(anno_dir_project, "rRNA_repeatmasker.bed")
tRNA_rmsk_path = paste0(anno_dir_project, "tRNA_repeatmasker.bed")

#set id for files
##########################################################################################
############### unique, all read count input - then merge
##########################################################################################
file_id=paste0(samplename,'vs',background,'_')

samples=read.table(samplemanifest,sep = '\t',header = T)
samples=rbind(samples[samples$group%in%samplename,],samples[samples$group%in%background,])


DB_N=read.delim(neg_DB, header=T, sep="\t",
                  stringsAsFactors = F,comment.char = '#')
DB_P=read.delim(pos_DB, header=T, sep="\t",
                stringsAsFactors = F,comment.char = '#')

DB=rbind(DB_N,DB_P)
DB=DB[duplicated(DB)==F,]

# ggplot(DB,aes(x=Fold,y=-log10(p.value)))+geom_point()
DB = DB %>%
  dplyr::rename(
    chr = seqnames)




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
DB_Oppo=DB[,c('chr','start','end','ID','ID','strand')] 
DB_Oppo$strand=gsub("\\+","pos",DB_Oppo$strand) 
DB_Oppo$strand=gsub("\\-","+",DB_Oppo$strand)
DB_Oppo$strand=gsub("pos","-",DB_Oppo$strand)

#run functions
DB_SameAnno = peak_calling(DB[,c('chr','start','end','ID','ID','strand')],"Same_")
DB_OppoAnno = peak_calling(DB_Oppo,"Oppo_")

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

intronexon_Same = IE_calling(DB_SameAnno,'Same',"Same_")
intronexon_Opposite = IE_calling(DB_OppoAnno,'Oppo',"Oppo_")

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
DB_attrib_Oppo = peak_attributes(rpmsk_Opposite,"Oppo_")
DBdataOut = merge(peak_attributes(rpmsk_Same,"Same_"),
                     DB_attrib_Oppo[,colnames(DB_attrib_Oppo)[!colnames(DB_attrib_Oppo) %in% 
                                                                    c("chr","start","end","strand" )]],
                     by='ID')


##########################################################################################
############### Assigning Clip peak attributes - merged strands
########################################################################################## 

DBdataOut=Assign_CLIP_attributes(DBdataOut)
  

##########################################################################################
############### Add annotations to junctions
##########################################################################################  
# if (JoinJunc==T)
if (F) {
  print("Add annotations to Junctions")
  
  DBdata2_anno=Join_Junction_Combine(DBdata2_anno,read_depth,FtrCount_out)
  
# }else{
  
  # #### collapse exon
  # if (condense_exon==T) {
  #   print("Collapse exon")
  #   
  #   DBdata2_anno=Collapse_Exon(DBdata2_anno)
  # }  
} 


##########################################################################################
## correct lnLc annotation - keep lnLc annotation above to hepl track source of lnc annotation
##########################################################################################

DBdataOut=correct_lnLc(DBdataOut)

#########################################################################################
############### Merge DiffBind and Annotations
##########################################################################################


DB[,paste0('Conc_',samplename)]=rowMeans(DB[,samples[samples$group%in%samplename,'sample']])
DB[,paste0('Conc_',background)]=rowMeans(DB[,samples[samples$group%in%background,'sample']])
DB=rename(DB,
          Length=width)

DBdata2_anno = merge(DB[,c('ID','chr','start','end','strand','Length','Fold','p.value',
                           paste0('Conc_',samplename),
                           paste0('Conc_',background),
                           c(samples$sample))],
                     DBdataOut[,!colnames(DBdataOut)%in%c('chr','start','end','strand')],
                     by='ID')
DBdata2_anno=rename(DBdata2_anno,
                    FC=Fold,
                    P_value=p.value)

for(x in samples$sample){
  colnames(DBdata2_anno)=gsub(x,paste0(x,'_NormCounts'),colnames(DBdata2_anno))
}

colnames(DBdata2_anno)=gsub(paste0('Conc_',samplename),paste0(samplename,'_AvgNormCounts'),colnames(DBdata2_anno))
colnames(DBdata2_anno)=gsub(paste0('Conc_',background),paste0(background,'_AvgNormCounts'),colnames(DBdata2_anno))


#########################################################################################
############### splice junctions
##########################################################################################
## 1) read in samplename annotations files
##    - get all peaks from IDmerge
## 2) identify Concensus peaks that overlap spliced peaks from IDmerge
## 3) combine Concensus peaks based on overlap with IDmerge peaks

if (join_junction) {
  

files=list.files(paste0(anno_dir_sample))
# smpl=paste(
  # paste0('(^',samples[samples$group%in%samplename,'sample'],'.*.txt',')'),
  # collapse="|")
# files[grep(smpl,files)]
files=paste0(samples[samples$group%in%samplename,'sample'],'_peakannotation_final_table.txt')

# check to make sure samples are selected
if (length(samples[samples$group%in%samplename,'sample'])!=length(files)){print('Could not find Sample Files');Join_Junction_Error
}else{
    
# get all spliced read locations  #sid= "WT1_fCLIP_peakannotation_final_table.txt" chr1:11080536-11080549_-,chr1:11080763-11080797_-
  splicedreads_Table=as.data.frame(c())
  for (sid in files) {
  anno=  fread(paste0(wd,'/05_annotation/',sid), header=T, sep="\t",stringsAsFactors = F,data.table=F,na.strings = "")
  anno$sample=sid
  # anno=paste(anno[is.na(anno$IDmerge)==F,c('ID','IDmerge')],collapse = ",")
  # splicedreads=paste0(splicedreads,",",anno)
  if('IDmerge'%in%colnames(anno)){
  splicedreads_Table=rbind(splicedreads_Table,anno[is.na(anno$IDmerge)==F,c('ID','IDmerge','sample')])
  }else{next}
      }
  if(nrow(splicedreads_Table[is.na(splicedreads_Table$ID)==F,])>0){

  
  splicedreads=paste0(splicedreads_Table$IDmerge,',',splicedreads_Table$IDmerge)
  splicedreads= unlist(str_split(splicedreads,pattern=','))%>%unique()
  
  splicedreads=separate(as.data.frame(splicedreads),1,sep="_",into=c('chr','strand'),remove = F)%>%separate(col='chr',sep=":|-",into=c("chr","start","end"))
  colnames(splicedreads)[1]="ID"
  splicedreads=splicedreads[duplicated(splicedreads)==F,]
  


##################################
## find spliced peaks that overlap with consensus peaks

DBdataOut.GR=GRanges(seqnames = as.character(DBdataOut$chr), 
                     ranges=IRanges(start = as.numeric(DBdataOut$start), 
                                    end = as.numeric(DBdataOut$end)),
                     strand = DBdataOut$strand,
                     ID=DBdataOut$ID)

splicedreads.GR=GRanges(seqnames = as.character(splicedreads$chr), 
                        ranges=IRanges(start = as.numeric(splicedreads$start), 
                                       end = as.numeric(splicedreads$end)),
                        strand = splicedreads$strand,
                        ID=splicedreads$ID)


xo=as.data.frame(GenomicRanges::findOverlaps(DBdataOut.GR, 
                                             splicedreads.GR, 
                                             type = "any",
                                             ignore.strand=F))
qh=as.data.frame(DBdataOut.GR[xo$queryHits],row.names = NULL)
colnames(qh)=paste0(colnames(qh),'_DBdata')

sh=as.data.frame(splicedreads.GR[xo$subjectHits],row.names = NULL)
colnames(sh)=paste0(colnames(sh),'_splicedreads')

DBdata_spliced=cbind(qh,sh)

DBdata_spliced_short=DBdata_spliced
system.time({
  DBdata_spliced_combOut=as.data.frame(c())
# for (x in 1:length(unique(DBdata_spliced_short$ID_DBdata))) {
while (nrow(DBdata_spliced_short)>0) {
  
  ## Idntify spliced read id for single concensus peak
  ID_single=DBdata_spliced_short[1,'ID_splicedreads']
  # ID_single=DBdata_spliced_short[x,'ID_splicedreads']
  
  IDmerge_single=splicedreads_Table[grepl(gsub("_","_\\\\",ID_single),splicedreads_Table$IDmerge),]
  
  # Identify peaks connected to single peak
    IDsearch=unique(unlist(str_split(paste(IDmerge_single$IDmerge,collapse = ','),pattern = ",")))

  
  # search all single sample peaks to identify Consensus peaks connected by splicing
  IDmerge_comb=DBdata_spliced_short[DBdata_spliced_short$ID_splicedreads%in%IDsearch, ]
  IDmerge_comb=IDmerge_comb[duplicated(IDmerge_comb)==F,]
  IDmerge_BDdata = sort(unique(IDmerge_comb$ID_DBdata))
  
  ## enter loop to look for additional peaks are connected to Peaks identified in IDmerege(Concensus peaks)
  ## this happens when samples have different peak coordinates combined into a single peak
  ### test if ('chr17:41502936-41502957_-'%in%IDmerge_BDdata) {xxxxx}

  IDsearchCurrent=""
  IDmerge_BDdataCurrent=""
  # while(length(setdiff(match_current,match_original))>0){
  while(length(setdiff(IDsearch,IDsearchCurrent))>0&
        length(setdiff(IDmerge_BDdata,IDmerge_BDdataCurrent))>0){

    IDsearchCurrent = IDsearch
    IDmerge_BDdataCurrent = IDmerge_BDdata
    
    # get all peaks that overlap with Consensus peaks idnetified above (IDmerge_BDdata)
    IDmerge_comb=DBdata_spliced_short[DBdata_spliced_short$ID_DBdata%in%IDmerge_BDdata, ]
    IDmerge_comb=IDmerge_comb[duplicated(IDmerge_comb)==F,]
     
    # get individule sample peaks connected to Consensus peak 
    IDsearch=unique(IDmerge_comb$ID_splicedread)

    # get Peak merge data from imported spliced reads
    IDsearch=splicedreads_Table[grepl(gsub("_","_\\\\",paste(IDsearch,collapse = "|")),splicedreads_Table$IDmerge),]
    IDsearch=unique(unlist(str_split(paste(IDsearch$IDmerge,collapse = ','),pattern = ",")))
    
    # look to see if new merge peaks overlap with additional Consensus peaks 
    IDmerge_comb=DBdata_spliced_short[DBdata_spliced_short$ID_splicedreads%in%IDsearch, ]
    IDmerge_comb=IDmerge_comb[duplicated(IDmerge_comb)==F,]
    IDmerge_BDdata = sort(unique(IDmerge_comb$ID_DBdata))
  }
    
  IDmerge_BDdata = separate(as.data.frame(IDmerge_BDdata),IDmerge_BDdata,sep = "_",into=c("chr","strand"),remove = F)%>%
    separate(col=chr,sep = ":|-",into=c("chr","start","end"),remove = F)
  
    
  d2=DBdata2_anno[DBdata2_anno$ID%in%IDmerge_BDdata$IDmerge_BDdata,]
  if (unique(IDmerge_BDdata$strand)=="+") {
    d5=d2[d2$start%in%min(d2$start),] 
    dmax=d2[d2$P_value==min(d2$P_value),]
    # dmax=d2[d2$Fold==max(d2$Fold),] ### highest FC
    # dmax=d2[apply(d2[,paste0('Conc_',c(samplename,background))],1,max)%in%max(apply(d2[,paste0('Conc_',c(samplename,background))],1,d2)),] ### highest Count
    if (nrow(dmax)>1) {dmax=dmax[dmax$start%in%min(dmax$start),]}
    
  } else if (unique(IDmerge_BDdata$strand)=="-") {
    d5=d2[d2$end%in%max(d2$end),]
    dmax=d2[d2$P_value==min(d2$P_value),]
    # dmax=d2[d2$Fold==max(d2$Fold),] ### highest FC
    # dmax=d2[apply(d2[,paste0('Conc_',c(samplename,background))],1,max)%in%max(apply(d2[,paste0('Conc_',c(samplename,background))],1,d2)),] ### highest Count
    if (nrow(dmax)>1) {dmax=dmax[dmax$start%in%max(dmax$start),]}
  }
  
  ### collapse all columns 
  d3= ( apply(d2 ,2, function(x){paste(unique(x[!is.na(x)]),collapse = ',')}))
  d3[d3==""]<-NA
  cols=c(paste0(c(samplename,background),'_AvgNormCounts'),paste0(samples$sample,'_NormCounts'),'Length')
  d3=as.data.frame(t(d3))
  d3[1,cols]=t(as.data.frame(colSums(d2[,cols])))   
  d3[,cols]=as.numeric(d3[,cols])
  d3$IDmerge=d3$ID
  if(grepl(',',d3$ID)==F){ERROR_No_conneccted_peaks}
  
  ## Select summary annotation based on Max(dmax) expression peak or 5'(d5) peak   
  danno=dmax # d5 dmax
  cols_AnnotSelect=c('FC','P_value','Same_Comb_type_exon','Same_Comb_type_ncRNA',
                     'Oppo_Comb_type_exon','Oppo_Comb_type_ncRNA','Comb_type_exon_Oppo')
  d3[,cols_AnnotSelect]= (danno[,cols_AnnotSelect])
  
  # Select ID for 5' most read or max counts
  d3$start=d5$start
  d3$end=d5$end
  d3$ID=paste0(d5$chr,":",d5$start,"-",d5$end,"_",d5$strand)   
  

  DBdata_spliced_combOut=rbind(DBdata_spliced_combOut,d3[1,])
  
  ### remove all linked peaks form list (shorten run time)
  DBdata_spliced_short=DBdata_spliced_short[DBdata_spliced_short$ID_splicedreads%in%IDsearch==F, ]
} 
})

## Checkpoint
if(nrow(DBdata_spliced_combOut[grepl(',',DBdata_spliced_combOut$IDmerge)==F,])>0){Error_not_Correctly_Joininging_Spliced_Peaks}


DBdata_spliced_combOut=DBdata_spliced_combOut[duplicated(DBdata_spliced_combOut)==F,]

DBdata_spliced_combOut=DBdata_spliced_combOut[is.na(DBdata_spliced_combOut$ID)==F,]
Junc_peaks=unique(unlist(str_split(paste(DBdata_spliced_combOut$IDmerge,collapse = ','),pattern = ",")))

DBdata2_anno$IDmerge=NA
DBdata2_anno=DBdata2_anno[DBdata2_anno$ID%in%Junc_peaks==F,]
DBdata2_anno=(rbind(DBdata2_anno,DBdata_spliced_combOut))

  }else{print('No spliced reads')}
}



}

## remove peaks with normalized reads <0
DBdata2_anno=DBdata2_anno[rowSums(DBdata2_anno[,paste0(c(samplename,background),'_AvgNormCounts')])>0,]


##########################################################################################
### Save complete  Annotation table
##########################################################################################

write.table(DBdata2_anno,paste0(out_dir,file_id,'diffBind_complete.txt'),sep = "\t",row.names = F)

##########################################################################################
############### Select Columns and rename for final output 
##########################################################################################




if ('IDmerge'%in% colnames(DBdata2_anno)) {
  DBdata3_anno=DBdata2_anno[,c("ID","IDmerge",'Length','FC','P_value',
                     #'chr',"start","end",#"strand",
                                     #'Counts_Unique','Counts_fracMM',
                                     paste0(c(samplename,background),'_AvgNormCounts'),paste0(samples$sample,'_NormCounts'),
                                     "Same_gene_name_comb","Same_ensembl_gene_id","Same_type_comb","Same_Repeat_comb",'Same_feature',"Same_Comb_type_exon","Same_Comb_type_ncRNA",
                     "Oppo_gene_name_comb","Oppo_ensembl_gene_id","Oppo_type_comb","Oppo_Repeat_comb",'Oppo_feature',"Oppo_Comb_type_exon","Oppo_Comb_type_ncRNA","Comb_type_exon_Oppo")]
} else {
  DBdata3_anno=DBdata2_anno[,c("ID",'Length','FC','P_value',
                     #'chr',"start","end",#"strand",
                      #               'Counts_Unique','Counts_fracMM',
                                      paste0(c(samplename,background),'_AvgNormCounts'),paste0(samples$sample,'_NormCounts'),
                                     "Same_gene_name_comb","Same_ensembl_gene_id","Same_type_comb","Same_Repeat_comb",'Same_feature',"Same_Comb_type_exon","Same_Comb_type_ncRNA",
                     "Oppo_gene_name_comb","Oppo_ensembl_gene_id","Oppo_type_comb","Oppo_Repeat_comb",'Oppo_feature',"Oppo_Comb_type_exon","Oppo_Comb_type_ncRNA","Comb_type_exon_Oppo")]
}
DBdata3_anno=rename(DBdata3_anno,
                       # Counts_Unique='Counts_Unique',
                       # Counts_Multimappers_Scaled='Counts_fracMM',
                       `Same Strand: Host_gene_name`="Same_gene_name_comb",
                       `Same Strand: Host_gene_ensembl_id`="Same_ensembl_gene_id",
                       `Same Strand: RNA_type`="Same_type_comb",
                       `Same Strand: Repeat_Type`="Same_Repeat_comb",
                       `Same Strand: Intron Exon`='Same_feature',
                       `Same Strand: RNA_type_Classification`="Same_Comb_type_exon",
                       `Same Strand: ncRNA_Classification`="Same_Comb_type_ncRNA",
                       `Opposite Strand: Host_gene_name`="Oppo_gene_name_comb",
                       `Opposite Strand: Host_gene_ensembl_id`="Oppo_ensembl_gene_id",
                       `Opposite Strand: RNA_type`="Oppo_type_comb",
                       `Opposite Strand: Repeat_Type`="Oppo_Repeat_comb",
                       `Opposite Strand: Intron Exon`='Oppo_feature',
                       `Opposite Strand: RNA_type_Classification`="Oppo_Comb_type_exon",
                       `Opposite Strand: ncRNA_Classification`="Oppo_Comb_type_ncRNA",
                       `Both Strand: RNA_type_Classification`="Comb_type_exon_Oppo"
)



##########################################################################################
### Write Output
##########################################################################################

#write out for junction annotation 
write.table(DBdata3_anno,paste0(out_dir,file_id,'diffBind_final_table.txt'),sep = "\t",row.names = F)


# print(session_info())
print(nrow(DBdata3_anno))

