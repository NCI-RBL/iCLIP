library(data.table)
library(dplyr)
library(tidyr)
#library(GenomicFeatures)
#library(rtracklayer)
# library(VariantAnnotation,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(GenomicRanges,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library("pheatmap")
# library(ggplot2,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library("viridis",quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(edgeR,quietly = T,verbose = F)
# library('GenomicFeatures',quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library('rtracklayer',quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(matrixStats,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(plyr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(tidyr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(fitdistrplus,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(stringr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
# library(reshape)
# library(stringi)
# library(plotly)
# library(GenomicRanges)
# library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
peak_type = args[1] #params$PeakIdnt
peak_unique = args[2] #peak_unique
peak_all = args[3] #params$PeaksFracMM
join_junction = args[4]
read_depth = args[5]
DEmethod = args [6]
sample_id = args[7]
nt_merge = args[8]
out_dir = args[9]


#testing information
testing="Y"
if(testing=="Y"){
  peak_type = "all"
  peak_unique = "/Volumes/data/iCLIP/marco/13_counts/WT_fCLIP_50_unique.txt"
  peak_all = "/Volumes/data/iCLIP/marco/13_counts/WT_fCLIP_50_all.txt"
  join_junction = "TRUE"
  read_depth=5
  DEmethod='MAnorm'
  sample_id = "KO"
  nt_merge = "50nt"
  out_dir = "/Volumes/data/iCLIP/marco/inprogress/"
  peaks_in = "/Volumes/data/iCLIP/marco/inprogress/peaks_KO_fCLIP.txt"
  ref_dir = "/Volumes/iCLIP/ref/annotations/"
  ref_species = "mm10" ### need better name, match to snakemake
  
  ### fix need to figure out how to keep all this info - maybe a config? dict?
  alias_path = paste0(ref_dir,ref_species,"/",ref_species,".chromAlias.txt")
  if(ref_species == "mm10"){
    gencode_path = paste0(ref_dir, "mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "/mm10/NCBI_RefSeq/GCF_000001635.26_GRCm38.p6_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"/mm10/Gencode_VM23/fromUCSC/KnownCanonical/KnownCanonical_GencodeM23_GRCm38.txt")
    intron_path = paste0(ref_dir, "/mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
    rmsk_path = paste0(ref_dir,"/mm10/repeatmasker/rmsk_GRCm38.txt")
    soyeong_flag = "Y" #Y or N
  } else if (ref_species == "hg38"){
    gencode_path = paste0(ref_dir,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"/hg38/Gencode_V32/fromUCSC/KnownCanonical/KnownCanonical_GencodeM32_GRCh38.txt")
    intron_path = paste0(ref_dir,"/hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
    rmsk_path = paste0(ref_dir,"/hg38/repeatmasker/rmsk_GRCh38.txt")
    soyeong_flag = "N" #always now (for now)
  } else{
    quit()
  }
}

##########################################################################################
############### Annotation info
##########################################################################################

###phil why are we naming NCBI (CM###) when it's GENCODE
#read in from reference files 
alias_anno=fread(alias_path, header=T, sep="\t",
                 stringsAsFactors = F,data.table=F,skip = "#",fill=TRUE)

#rename cols
colnames(alias_anno)=c('chr','alias1','aliasNCBI','Refseq')

#create second NCBI col, replace with [_] of chr col
alias_anno$aliasNCBI2=alias_anno$aliasNCBI
alias_anno[-grep('_',alias_anno$chr),'aliasNCBI2']=alias_anno[-grep('_',alias_anno$chr),'chr']

###phil - what is this doing?
alias_anno[grep('_',alias_anno$aliasNCBI2),'aliasNCBI2']=alias_anno[grep('_',alias_anno$aliasNCBI2),'alias1']

#copy refseq col
###phil - what is this doing?
alias_anno$Refseq2=alias_anno$Refseq
alias_anno[(grepl('_',alias_anno$Refseq2)|(alias_anno$Refseq2%in%""))==F,'Refseq2']=alias_anno[(grepl('_',alias_anno$Refseq2)|(alias_anno$Refseq2%in%""))==F,'aliasNCBI']

###phil - why is this different between the two types? 
###phil - why are we making the aliasNCBI2 column? (alias_anno[1:20,])
###phil - refseq2 column is the same as NCBI?
# if (species=='mm10'){
#   alias_anno=fread(paste0(ref_dir,"/mm10/mm10.chromAlias.txt"), header=T, sep="\t",
#               stringsAsFactors = F,data.table=F,skip = "#",fill=TRUE)
#   colnames(alias_anno)=c('chr','alias1','aliasNCBI','Refseq')
#   alias_anno$aliasNCBI2=alias_anno$aliasNCBI
#   alias_anno[-grep('_',alias_anno$chr),'aliasNCBI2']=alias_anno[-grep('_',alias_anno$chr),'chr']
#   
# } else if (species=='hg38') {
#   alias_anno=fread(paste0(Ref,"/hg38/hg38.chromAlias.txt"), header=T, sep="\t",
#               stringsAsFactors = F,data.table=F,skip = "#",fill=TRUE)
#   colnames(alias_anno)=c('chr','alias2','aliasNCBI',"Refseq")
#   alias_anno$aliasNCBI2=alias_anno$aliasNCBI
#   alias_anno[-grep('_',alias_anno$chr),'aliasNCBI2']=alias_anno[-grep('_',alias_anno$chr),'chr']
#   alias_anno[grep('_',alias_anno$aliasNCBI2),'aliasNCBI2']=alias_anno[grep('_',alias_anno$aliasNCBI2),'alias2']
#   alias_anno$Refseq2=alias_anno$Refseq
#   alias_anno[(grepl('_',alias_anno$Refseq2)|(alias_anno$Refseq2%in%""))==F,'Refseq2']=alias_anno[(grepl('_',alias_anno$Refseq2)|(alias_anno$Refseq2%in%""))==F,'aliasNCBI']
# } else{
#   quit(status = 1)
# }

##########################################################################################
############### Peak info
##########################################################################################
#read in peaks file (for testing, will be called in RMD?) 
peaks = read.csv(peaks_in)

#merge peaks with ref
peaks_alias = merge(peaks,alias_anno[,c('chr','aliasNCBI2')],by.x='chr',by.y='aliasNCBI2',all.x=T)

###phil why are we doing both? 
#clean
peaks_alias=peaks_alias[,colnames(peaks_alias)%in%"chr"==F]
peaks_alias=peaks_alias[is.na(peaks_alias$chr)==F,]

#rename col
colnames(peaks_alias)=gsub('chr.y','chr',colnames(peaks_alias))

###phil this col is the same as ID - is it necessary? 
#add IDmod
peaks_alias$IDmod=paste0(peaks_alias$chr,":",peaks_alias$start,"-",peaks_alias$end)

###phil why are we switching positive and negative?
peaks_oppo=peaks 
peaks_oppo$strand=gsub("\\+","pos",peaks_oppo$strand) 
peaks_oppo$strand=gsub("\\-","+",peaks_oppo$strand)
peaks_oppo$strand=gsub("pos","-",peaks_oppo$strand)

##########################################################################################
############### GENCODE ANNOTATION
##########################################################################################
ref_gencode = fread(gencode_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

#rename cols
ref_gencode = ref_gencode %>% 
  rename(
    chr = seqname,
    ensembl_gene_id = gene_id,
    external_gene_name = gene_name
  )

#remove version
removeVersion <- function(ids){
  return(unlist(lapply(stringr::str_split(ids, "[.]"), "[[",1)))
}
ref_gencode$transcript_id=removeVersion(ref_gencode$transcript_id)
ref_gencode$ensembl_gene_id=removeVersion(ref_gencode$ensembl_gene_id)

#merge with alias
ref_gencode = merge(ref_gencode,
                    alias_anno[,c('chr','aliasNCBI2')],
                    by.x='chr',by.y='aliasNCBI2',all.x=T)

#remove chr col
ref_gencode = select(ref_gencode, -chr)

#rename chr.y col
colnames(ref_gencode) = gsub('chr.y','chr',colnames(ref_gencode))

#remove rows with missing data
ref_gencode = ref_gencode[!is.na(ref_gencode$chr), ]

#remove TEC genes
ref_gencode = subset(ref_gencode, transcript_type != "TEC") 
ref_gencode = subset(ref_gencode, gene_type != "TEC") 

#combine types into more general categories
ref_gencode$gene_type_ALL = ref_gencode$gene_type

## combine all Pseudogenes
ref_gencode[grep('pseudogene',ref_gencode$gene_type),'gene_type']='pseudogene'

# subset snoRNA before merge
ref_gencode_sno=subset(ref_gencode, gene_type == 'snoRNA')
ref_gencode_sno$gene_type = 'snoRNA_refGencode'
ref_gencode_sno$gene_type_ALL = ref_gencode_sno$gene_type

###phil why are we merging these all into ncRNA? Can we leave as separate types?
## combine all ncRNA
ReplaceType <-function(df_in,col_in,list_in,replace_id){
  for (id in list_in){
    df_in[,col_in]=gsub(id,replace_id,df_in[,col_in])
  } 
  return(df_in)
}

type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA","tRNA","ribozyme")
ref_gencode = ReplaceType(ref_gencode,"gene_type", type_list,"ncRNA")

#create subsets
###phil FTR is never used
ref_gencode_t = subset(ref_gencode, feature == 'transcript')
ref_gencode_e = subset( ref_gencode, feature == "exon")
#ref_gencode_FTR = ref_gencode[ref_gencode$feature%in%c("3UTR","5UTR",'CDS','UTR'),]

##########################################################################################
############### REFSEQ ANNOTATION
##########################################################################################
###phil the ref being read has NA in all cols past protein_in so data just gets removd with 
#filtering 
#head /data/RBL_NCI/iCLIP/ref/CLIP_Anno/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt 
ref_refseq = fread(refseq_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

ref_refseq = ref_refseq %>% 
  rename(
    gene_type = gene_biotype,
    ensembl_gene_id = gene_id,
    external_gene_name = gene
  )

###phil why remove duplicated cols?
#remove dup cols
#ref_refseq=ref_refseq[,duplicated(colnames(ref_refseq))==F]

#merge with alias
ref_refseq = merge(ref_refseq, 
                   alias_anno[,c('chr','Refseq2')],
                   by.x='seqname',by.y='Refseq2',all.x=T)

#remove rows with missing data
###phil chr col is all NA so everything is deleted
ref_refseq = ref_refseq[!is.na(ref_refseq$chr), ]

#subset for genes only
ref_refseq = subset(ref_refseq, feature=="gene")

##remove TEC genes
###phil - why only remove TEC in geneypte and not also in transcript? (like above)
#ref_refseq=ref_refseq[(ref_refseq$gene_type%in%'TEC')==F,]  
ref_refseq = subset(ref_refseq, gene_type == "TEC")   
ref_refseq = subset(ref_refseq, transcript_type == "TEC") 

#combine types into more general catagories
ref_refseq$gene_type_ALL=ref_refseq$gene_type

#remove pseudogenes
ref_refseq[grep('pseudogene',ref_refseq$gene_type),'gene_type']='pseudogene'

#### refseq - Select contigs not in Gencode
###phil we haven't selected for transcript only here (like we did in gencode, so the
#list will be inherently different)
#select contigs not in gencode  
ref_refseq_only = ref_refseq[ref_refseq$chr %in% unique(ref_gencode_t$chr)==F,]

###phil this annotation is never used?
## combine all ncRNA
type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA","tRNA","ribozyme")
ref_refseq_only = ReplaceType(ref_refseq_only,"gene_type", type_list,"ncRNA")

###phil why is there such a difference between the two?
###phil gene_type col name is the same as the change
# if (species=='mm10') {
#     ref_refseq=fread(paste0(Ref,"/mm10/NCBI_RefSeq/GCF_000001635.26_GRCm38.p6_genomic.gtf.txt"), header=T, sep="\t",stringsAsFactors = F,data.table=F)
#     colnames(ref_refseq)[colnames(ref_refseq)%in%'gene_biotype']='gene_type'
#     ref_refseq=ref_refseq[,duplicated(colnames(ref_refseq))==F]
#   } else if (species=='hg38') {
#     ref_refseq=fread(paste0(Ref,"/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt"), header=T, sep="\t",stringsAsFactors = F,data.table=F)
#     colnames(ref_refseq)[colnames(ref_refseq)%in%'gene_biotype']='gene_type'

###phil cols can't be duplicated?
#     ref_refseq=ref_refseq[,duplicated(colnames(ref_refseq))==F]
#     ref_refseq =merge(ref_refseq,alias_anno[,c('chr','Refseq2')],by.x='seqname',by.y='Refseq2',all.x=T)
#     ref_refseq=ref_refseq[is.na(ref_refseq$chr)==F,]
#     colnames(ref_refseq)[colnames(ref_refseq)%in%'gene_id']='ensembl_gene_id'
#     colnames(ref_refseq)[colnames(ref_refseq)%in%'gene']='external_gene_name'
#     colnames(ref_refseq)[colnames(ref_refseq)%in%'gene_type']='gene_type'
#     ref_refseq=ref_refseq[(ref_refseq$feature)%in%'gene',]
#     
#     ###phil never use this code again
#     mm10Refseq_snoRNA=ref_refseq[ref_refseq$gene_type%in%'snoRNA',]
#     mm10Refseq_snoRNA$gene_type='snoRNA_RefSeq'
#     mm10Refseq_snoRNA$gene_type_ALL=mm10Refseq_snoRNA$gene_type
#     
#     #### refseq - Select contigs not in Gencode
#     ref_refseq=ref_refseq[ref_refseq$chr%in%unique(ref_gencode$chr)==F,]
#     ref_refseq=ref_refseq[(ref_refseq$gene_type%in%'TEC')==F,]  
#     ref_refseq$gene_type_ALL=ref_refseq$gene_type
#     ref_refseq[grep('pseudogene',ref_refseq$gene_type),'gene_type']='pseudogene'
#     
#     ref_refseq$gene_type=gsub('miRNA','ncRNA',ref_refseq$gene_type) ### fix
#     ref_refseq$gene_type=gsub('miscRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('misc_RNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('piRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('rRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('siRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('snRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('snoRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('tRNA','ncRNA',ref_refseq$gene_type)
#     ref_refseq$gene_type=gsub('ribozyme','ncRNA',ref_refseq$gene_type)
#   }
########### UPDATE ###############


##########################################################################################
############### canonical paths
##########################################################################################
canonical=fread(canonical_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

#remove version
canonical$transcript=removeVersion(canonical$transcript)
canonical$protein=removeVersion(canonical$protein)

#merge canonical with annotation
canonical = merge(canonical,
                  alias_anno[,c('chr','aliasNCBI2')],by.x='#chrom',by.y='aliasNCBI2',all.x=T)

#remove #chrom column
canonical = select(canonical, -c("#chrom"))

##########################################################################################
############### introns
##########################################################################################

introns=fread(intron_path, 
              header=F, sep="\t",stringsAsFactors = F,data.table=F, 
              col.names = c('chr','start','end','attribute','V5','strand'))

#remove Non-Chromosome contigs
introns=introns[grepl("_",introns$chr)==F,]

#split attribute col
introns=separate(introns,
                 attribute,
                 into=c('transcript_id','feature','exon_number','level','chr2','intronnumber','dir'),
                 remove = T,sep = "_")

#split out transcript id
introns$transcript_id = removeVersion(introns$transcript_id)

###phil why are we adding 1 to start and exon? 
introns$start=introns$start+1
introns$exon_number=as.numeric(introns$exon_number)+1

#merge with alias
introns =merge(introns,
               alias_anno[,c('chr','aliasNCBI2')],by.x='chr',by.y='aliasNCBI2',all.x=T)

#remove chromosome col
introns = select(introns, -c("chr"))
colnames(introns)=gsub('chr.y','chr',colnames(introns))

##########################################################################################
############### binding exon, intron
##########################################################################################
intron_exon=rbind(ref_gencode_e[,c('chr','feature','start','end','strand','transcript_id','exon_number')],
                  introns[,c('chr','feature','start','end','strand','transcript_id','exon_number')])
intron_exon$ID=paste0(intron_exon$chr,':',intron_exon$start,'-',intron_exon$end)

###phil why is this flag set? none of the following code is used with flag  
calcIntron=0
if (calcIntron==1){
  gtf <- makeTxDbFromGFF(paste0(Ref,"/hg38/Gencode_V32/fromUCSC/KnownGene/hg38.KnownGene.gtf")) #change me!
  exons <- exonsBy(gtf, by="gene")
  
  #make introns
  exons <- reduce(exons)
  exons <- exons[sapply(exons, length) > 1]
  
  introns <- lapply(exons, function(x) {
    
    #Make a "gene" GRange object
    gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)), end=max(end(x))),
                 strand=strand(x)[1])
    db = disjoin(c(x, gr))
    ints = db[countOverlaps(db, x) == 0]
    
    #Add an ID
    if(as.character(strand(ints)[1]) == "-") {
      ints$exon_id = c(length(ints):1)
    } else {
      ints$exon_id = c(1:length(ints))
    }})
  
  introns <- GRangesList(introns)
  as.data.frame(introns)
}

#merge exon/intron with ref_gencode_t
###phil - we merged ref_gencode_exons with the introns from file
#now we're merging the ref_gencode_transcript... why all the merges
intron_exon=merge(intron_exon,
                  ref_gencode_t[,c('transcript_id','gene_type')],by='transcript_id',all.x=T)

##########################################################################################
############### REPEATMASKER
##########################################################################################
# get  repeat regions and extra small RNA locations
###phil why is this being preset
#newRmsk=1
#if (newRmsk==1) {
rmsk_GRCm38=fread(rmsk_path, header=T, sep="\t",stringsAsFactors = F,
                  data.table=F)

#dont need second variable
#rmsk_GRCm38all=rmsk_GRCm38
#}
# 
# if(species=="mm10" & soyeong_flag=="Y"){
#   rRNA_BK00964=SYImport(rRNA_BK00964_path,col_list,c('V1','V2','V3','V4','V5','V6'),'rRNA')
#   rRNA_BK00964$strand="*"

# sncRNA_sy=SYImport(sncRNA_sy_path)
# sncRNA_sy=separate(sncRNA_sy,col = "V9",into = c("gene","trans",'x'),sep = ";")
# sncRNA_sy$gene=gsub("gene_id ","",sncRNA_sy$gene);sncRNA_sy$trans=gsub("transcript_id ","",sncRNA_sy$trans)
# sncRNA_sy=as.data.frame(sncRNA_sy)
# sncRNA_sy=sncRNA_sy[,((colnames(sncRNA_sy)%in%'x')==F)]
# colnames(sncRNA_sy)=c('chr','source','gloc','start','end','V6','strand','V8','Gene_Ensemble','Transcript_id')
# sncRNA_sy$type='sncRNA'
# colnames(sncRNA_sy)[colnames(sncRNA_sy)%in%'Gene_Ensemble']='name'
#}

# yRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_YRNA.bed")
# srpRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_srpRNA.bed")
# tRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_tRNA.bed")
# scRNA_sy_path  = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_scRNA.bed")
# SKRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_7SKRNA.bed")
# ###### rRNA pseudogenes and 5S	mm10_rRNA.bed and mouse_gencode_rRNA.gtf.1
# rRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_rRNA.bed") 
# ###### rRNA BK00964.3 'https://www.ncbi.nlm.nih.gov/nuccore/NR_046233.1'
# ###phil is there a reason this is in a different dir?
# rRNA_BK00964_path = paste0(ref_dir,"/mm10/AdditionalAnno/BK000964.3_TPA_rRNA_repeats2.bed.txt")
# ##### sncRNA (mrp RNA, Rpph1, vault RNA)	mm10_snc.gtf
# sncRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_snc.gtf")
# snRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mouse_gencode_snRNA.gtf.1")
# snoRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mouse_gencode_sno.gtf.1")
# ###### miRNA		mouse_gencode_miRNA.gtf.1 (Gencode)
# miRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mouse_gencode_miRNA.gtf.1")
# ###### rRNA pseudogenes and 5S	mm10_rRNA.bed and mouse_gencode_rRNA.gtf.1 (Gencode)
# rRNA_gtf_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mouse_gencode_rRNA.gtf.1")
# ###### linc RNA	mouse_gencode_linc.gtf.1 (Repeatmasker)
# lincRNA_sy_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mouse_gencode_linc.gtf.1")
# SY_LTR_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_LTR.bed")
# SY_DNA_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_DNA.bed")
# SY_sat_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_sat.bed")
# SY_SR_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_simple.bed")
# SY_LC_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_LC.bed")
# SY_other_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_other.bed")
# SY_unknown_path = paste0(ref_dir,"/mm10/AdditionalAnno/from_Soyeong/mm10_annotation/mm10_unknown.bed")

gencodeLNCRNA<-function(){
  ###phil in this subset we're not just using the transcripts (ref_gencode_t) but everything... why?
  lncRNA_ref_gencode_exon = subset(ref_gencode, gene_type_ALL == "lncRNA" & feature == "exon")
  
  #select exons found in intron list
  lncRNA_ref_gencode_exon = lncRNA_ref_gencode_exon[lncRNA_ref_gencode_exon$transcript_id 
                                                    %in% unique(introns$transcript_id),]
  
  #select introns found in exon list
  lncRNA_ref_gencode_intron = introns[introns$transcript_id 
                                      %in% unique(lncRNA_ref_gencode_exon$transcript_id),]
  
  #merge intron list and exon list
  ###phil why are we merging these when the list above is made with only those introns 
  #found in exon list? 
  #since _exon has duplicate transcriptIDs (14K ids, 46K rows), this merge is creating
  #duplicate rows for each duplicate id (exon = 46K, intron = 32K, merge = 154K)
  lncRNA_ref_gencode_intron = merge(lncRNA_ref_gencode_intron,
                                    lncRNA_ref_gencode_exon[,
                                                            c('transcript_id','ensembl_gene_id',
                                                              'external_gene_name','gene_type',
                                                              'gene_type_ALL','transcript_type',
                                                              'transcript_name','score')],
                                    by='transcript_id',all.x=T)
  
  #select anything that is not an intron (in intron list)
  lncRNA_ref_gencode_notInt = ref_gencode[!(ref_gencode$transcript_id
                                            %in% unique(introns$transcript_id)),]
  
  #This gets overridden in the next command
  #subset ref for lncRNA
  #lncRNA_ref_gencode = subset(ref_gencode, gene_type_ALL == 'lncRNA')
  
  #bind exon (all ref_genes, subset for lncrna and exon, subset for transcript 
  #ID's found in intron list) + intron (introns found in exon list)
  ###phil why are we merging these together when it seems that they are already merges of 
  #one another? this is the opposite of the intron merge (merged exon into intron) above
  lncRNA_ref_gencode=rbind(lncRNA_ref_gencode_exon[,c('chr','start','end','strand',
                                                      'ensembl_gene_id','external_gene_name',
                                                      'gene_type','gene_type_ALL','feature',
                                                      'transcript_id','transcript_type',
                                                      'transcript_name','score')],
                           lncRNA_ref_gencode_intron[,c('chr','start','end','strand',
                                                        'ensembl_gene_id','external_gene_name',
                                                        'gene_type','gene_type_ALL','feature',
                                                        'transcript_id','transcript_type',
                                                        'transcript_name','score')])
  
  #already done this with subses in original exon and intron list
  #lncRNA_ref_gencode$gene_type_ALL='lincRNA'
  
  #add exon/intron to type
  lncRNA_ref_gencode$gene_type_ALL = paste0(lncRNA_ref_gencode$gene_type_ALL,
                                            '-',lncRNA_ref_gencode$feature)
  
  #merge above list with not introns
  lncRNA_ref_gencode = rbind(lncRNA_ref_gencode[,c('chr','start','end','strand',
                                                   'ensembl_gene_id','external_gene_name',
                                                   'gene_type','gene_type_ALL','feature',
                                                   'transcript_id','transcript_type',
                                                   'transcript_name','score')],
                             lncRNA_ref_gencode_notInt[,c('chr','start','end','strand',
                                                          'ensembl_gene_id','external_gene_name',
                                                          'gene_type','gene_type_ALL','feature',
                                                          'transcript_id','transcript_type',
                                                          'transcript_name','score')])
  head(lncRNA_ref_gencode)
  return(lncRNA_ref_gencode)
}


#select rows in ref that are not in above list
###phil why are we doing this? 
#ref_gencode_Q = ref_gencode_t[!(ref_gencode_t$transcript_id 
#                                %in% unique(lncRNA_ref_gencode$transcript_id)),]

#subset
#sRNA_ref_gencode = subset(ref_gencode_Q, transcript_type == 'sRNA')
#misc_RNA_ref_gencode = subset(ref_gencode_Q, transcript_type == 'misc_RNA')
#ribozyme_RNA_ref_gencode = subset(ref_gencode_Q, transcript_type =='ribozyme')

gencodeAnno<-function(rowid){
  transcript_list = c("miRNA","snoRNA")
  #introns
  if(rowid=="intron"){
    df_sub = introns
    content = paste0(unique(df_sub$feature),collapse = ', ')
  } else if (rowid %in% transcript_list){
    df_sub = subset(ref_gencode_t, transcript_type == rowid)
    content = paste0(unique(df_sub$transcript_type),collapse = ', ')
  } else if (rowid == "lncRNA"){
    df_sub = gencodeLNCRNA() %>%
      subset(gene_type_ALL == "lncRNA")
    content = paste0(unique(df_sub$transcript_type),collapse = ', ')
  }
  
  write.table(df_sub, file=paste0(outdir,"/annotation/",rowid,".bed"), 
              sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
  
  
  return(list(df_sub,content))
}

rmskAnno<-function(rowid){
  family_list = c("rRNA","snRNA", "srpRNA","tRNA")
  class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','Unknown','Other')
  
  sycol=c('chr','start','end','name','swScore','strand')
  mm10col=c('chr','start','end','transcript_name','score','strand')
  rmskcol=c('chr','start','end','name','swScore','strand')
  
  if(rowid %in% family_list){
    df_sub = subset(rmsk_GRCm38, repFamily == rowid) %>%
      select(c('genoName','genoStart','genoEnd','repName','swScore','strand')) %>%
      setnames(c('chr','start','end','name','swScore','strand'))
  } else if (rowid %in% class_list){
    df_sub = subset(rmsk_GRCm38, repClass == rowid)
  } else if(rowid== "7SK RNA"){
    df_sub = subset(rmsk_GRCm38, repFamily == 'RNA') %>%
      select(c('genoName','genoStart','genoEnd','repName','swScore','strand')) %>%
      setnames(c('chr','start','end','name','swScore','strand'))
    df_sub$type='7SKRNA'
  } else if (rowid == "LINE SINE"){
    df_sub = subset(rmsk_GRCm38, repClass == c("LINE","SINE"))
  } else if (rowid == "rRNA_DNA"){
    df_sub = subset(rmsk_GRCm38, repClass == "DNA")
  } else if (rowid == "scRNA"){
    df_sub = subset(rmsk_GRCm38, repFamily == 'scRNA' & grepl("HY",rmsk_GRCm38$repName)==FALSE) %>%
      select(c('genoName','genoStart','genoEnd','repName','swScore','strand')) %>%
      setnames(c('chr','start','end','name','swScore','strand'))
    df_sub$type='scRNA'
  } else if (rowid == "yRNA"){
    df_sub = subset(rmsk_GRCm38, repFamily == 'scRNA' & grepl("HY",rmsk_GRCm38$repName)==TRUE) %>%
      select(c('genoName','genoStart','genoEnd','repName','swScore','strand')) %>%
      setnames(c('chr','start','end','name','swScore','strand'))
    df_sub$type='yRNA'
    
    write.table(df_sub,file=paste0(out_dir,"/annotation/",rowid,".bed"), 
                sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
    
  }
  content = paste0(unique(df_sub$name),collapse = ", ")
  return(list(df_sub,content))
}
   
SYAnno<-function(rowid,ref_species){
  #read bedfile
  file_name = ref_table[rowid,"SY_1"]
  df_sub = read.table(paste0(ref_dir,"mm10/additional_anno/",file_name))
  
  #if the file is a .gtf.1 then filter
  if (grepl(".gtf.1",file_name)){
    df_sub = subset(df_sub, V3 == 'gene')
    contents = ""
    
    #if the file is a .gtf then filter
  } else if (grepl(".gtf",file_name)){
    df_sub = separate(df_sub,col = "V9",into = c("gene","trans",'x'),sep = ";")
    df_sub$gene = gsub("gene_id ","",df_sub$gene)
    df_sub$trans = gsub("transcript_id ","",df_sub$trans)
    df_sub = df_sub[,((colnames(df_sub)%in%'x')==F)]
    contents = ""
    
    #all other files 
  } else{
    contents = paste0(unique(df_sub$V4),collapse = ', ')
  }
  
  write.table(df_sub,
              file=paste0(outdir,"/annotation/",rowid,".bed"), 
              sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)

  return(list(df_sub,content))
  
}

###phil this entire df is set - nothing is new related to the project is being 
#integrated - if this is the case then maybe we should create it once and read the entire
#thing in?
##########################################################################################
############### Create annotation_output table 
##########################################################################################
###phil -  removed rRNA_gencode and rRNA_rmsk - this is the only type split into two
###fix - eventaully do not prest this table
#rnames=c('45S','chrM','yRNA','snRNA','snoRNA','srpRNA','tRNA',
#        '7SK RNA','scRNA','sncRNA','miRNA','rRNA',
#         'rRNA_DNA','lncRNA','lincRNA','mRNA','LINE SINE','LTR','DNA',
 #        'satellites','Simple Repeats','Low Complexity','other',
#         'unknown','Introns')
#cnames=c('contents','source','Description','notes','Annotation_file')
annotation_output=data.frame()

###phil - review the reftable to fill in missing data
#input data from annotation df
reftable_path = paste0("/Volumes/sevillas2/git/iCLIP/workflow/scripts/reftable.txt")
ref_table = read.csv(reftable_path, header = TRUE, sep="\t", row.names = 1)

for (rowid in rownames(ref_table)){
  
  #determine source to use
  ref_selection =  ref_table[rowid,"selection"]
  print(paste0(rowid,"-",ref_selection))
  ### Gencode annotations
  if(ref_selection=="Gencode"){
    return_list = gencodeAnno(rowid)
    df_sub = return_list[[1]]
    content = return_list[[2]]
    
    ### REFSEQ annotations  
  } else if (ref_selection == "RefSeq"){
    print("add refseq")
    
    ### Repeatmasker annotations
  } else if (ref_selection == "Repeatmasker"){
    return_list = rmskAnno(rowid)
    df_sub = return_list[[1]]
    content = return_list[[2]]
    
    ###Soyeong references
  } else if (ref_selection == "SY"){
    return_list = SYAnno(rowid)
    df_sub = return_list[[1]]
    content = return_list[[2]]
    
    ###Invalid entry
  } else{
    df_sub=data.frame()
    content=""
  }
  
  #print into to output
  annotation_output[rowid,"contents"]  = content
  annotation_output[rowid,"count"] = nrow(df_sub)  
  annotation_output[rowid,'source'] = ref_selection
  annotation_output[rowid,'description'] = ref_table[rowid,"description"]
}

#hard code contents
if (ref_species=='mm10') {
  annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12 and various predicted genes'
  annotation_output['rRNA_gencode','contents']='5S, 5.8s, predicted gene'
} else if (ref_species=='hg38') {
  annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12'
  annotation_output['rRNA_gencode','contents']='5S, 5.8s, predicted gene'
}
write.csv(annotation_output,paste0(out_dir,"annotation_in_progress.csv"))

###phil - this is not a listed family, can't add
#unique(rmsk_GRCm38$repFamily)
# SKRNA_rmsk = rmskSelection(rmsk_GRCm38, '7SKRNA')
# annotation_output['7SK RNA','Repeatmasker']=nrow(SKRNA_rmsk)
# annotation_output['7SK RNA','contents']=paste0(unique(SKRNA_rmsk$name),collapse = ', ')


###phil why are these difference from one another?
# if (species=='mm10') {
###phil why are we hardcoding these?
#     annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12 and various predicted genes'
###phil why is this and the next the only ones pulling from soyeong?
#     annotation_output['tRNA','contents']=paste0(unique(tRNA_sy$type),collapse = ', ')
#     annotation_output['sncRNA','contents']=paste0('3 annotations: ',paste0(unique(sncRNA_sy$name),collapse = ', '))
#     annotation_output['miRNA','contents']=paste0(unique(miRNA_ref_gencode$transcript_type),collapse = ', ')
#     annotation_output['rRNA_gencode','contents']='5S, 5.8s, predicted gene'
#     annotation_output['rRNA_rmsk','contents']=paste0(unique(rRNA_rmsk[order(rRNA_rmsk$name),'name']),collapse = ', ')
#     annotation_output['rRNA_DNA','contents']=paste0(unique(rRNA_BK00964$name),collapse = ', ')
#     annotation_output['lncRNA','contents']=paste0('Various ',paste0(unique(lncRNA_ref_gencode[lncRNA_ref_gencode$gene_type_ALL%in%'lncRNA','transcript_type']),collapse = ', '))
#     annotation_output['lincRNA','contents']=paste0('Various ',paste0(unique(lncRNA_ref_gencode[lncRNA_ref_gencode$gene_type_ALL%in%'lncRNA'==F,'transcript_type']),collapse = ', '))
#   }
#   if (species=='hg38') {
#     annotation_output['yRNA','contents']=paste0(unique(YRNA_rmsk$name),collapse = ', ')
#     annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12'
#     annotation_output['tRNA','contents']=paste0(unique(tRNA_rmsk$type),collapse = ', ')
#     annotation_output['scRNA','contents']=paste0(unique(scRNA_rmsk$name),collapse = ', ')
#     annotation_output['miRNA','contents']=paste0(unique(miRNA_ref_gencode$transcript_type),collapse = ', ')
#     annotation_output['rRNA_gencode','contents']='5S, 5.8s, predicted gene'
#     annotation_output['rRNA_rmsk','contents']=paste0(unique(rRNA_rmsk[order(rRNA_rmsk$name),'name']),collapse = ', ')
#     annotation_output['lncRNA','contents']=paste0('Various ',paste0(unique(lncRNA_ref_gencode[lncRNA_ref_gencode$gene_type_ALL%in%'lncRNA','transcript_type']),collapse = ', '))
#     annotation_output['lincRNA','contents']=paste0('Various ',paste0(unique(lncRNA_ref_gencode[lncRNA_ref_gencode$gene_type_ALL%in%'lncRNA'==F,'transcript_type']),collapse = ', '))
#   }  

#write table



# #replaced with code above
# annotation_output['yRNA','Repeatmasker']=nrow(rmskSelection(rmsk_GRCm38[grepl("HY",rmsk_GRCm38$repName)==TRUE,], 
#                                                          'scRNA', 'yRNA'))
# annotation_output['scRNA','Repeatmasker']=nrow(rmskSelection(rmsk_GRCm38[grepl("HY", rmsk_GRCm38$repName)==FALSE,],'scRNA'))  
# annotation_output['snRNA','Repeatmasker']=nrow(rmskSelection(rmsk_GRCm38,'snRNA'))
# 
# annotation_output['srpRNA','Repeatmasker']=nrow(rmskSelection(rmsk_GRCm38,'srpRNA'))
# annotation_output['srpRNA','contents']=paste0(unique(srpRNA_rmsk$name),collapse = ', ')
# 
# annotation_output['tRNA','Repeatmasker']=nrow(rmskSelection(rmsk_GRCm38,'tRNA'))
# annotation_output['rRNA_rmsk','Repeatmasker']=nrow(rmskSelection(rmsk_GRCm38,'rRNA'))
# annotation_output['LISI','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass %in% c('LINE','SINE')))
# annotation_output['LTR','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'LTR'))
# annotation_output['DNA','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'DNA'))
# annotation_output['Satalites','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'Satellite'))
# annotation_output['Simple Repeats','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'Simple_repeat'))
# annotation_output['Low Complexity','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'Low_complexity'))
# annotation_output['Other repeats','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'Other'))  
# annotation_output['Unknown repeats','Repeatmasker']=nrow(subset(rmsk_GRCm38, repClass == 'Unknown'))
# annotation_output$Comb=rowSums(annotation_output[,c('Gencode','Repeatmasker')])

###phil is there a time when we would not want to write these?
#if (WriteClassTable==T) {
###phil is there a reason for the misc folder being created

# if (species=='mm10'){
#   write.table(YRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/yRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(srpRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/srpRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(SKRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/SKRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(scRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/scRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/rRNA_rmsk.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_BK00964[,rmskcol],file=paste0(outdir,"/annotation/rRNA_BK00964.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   
#   write.table(snRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/snRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(snoRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/snoRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(miRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/mirNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/rRNA_gencode.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(lncRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/lncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   write.table(tRNA_sy[,rmskcol],file=paste0(outdir,"/annotation/tRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(sncRNA_sy[,c('chr','start','end','name','V6','strand')],file=paste0(outdir,"/annotation/sncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   rnames_ncRNA=c('yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','sncRNA','miRNA','rRNA_gencode','rRNA_rmsk','rRNA_DNA','lncRNA','lincRNA')
#   cnames_ncRNA=c('source','contents','Description') 
# } else if (species=='hg38'){
#   write.table(YRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/yRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(srpRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/srpRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(tRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/tRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(SKRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/SKRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(scRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/scRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_rmsk[,rmskcol],file=paste0(outdir,"/annotation/rRNA_rmsk.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   write.table(snRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/snRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(snoRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/snoRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(miRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/mirNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/rRNA_gencode.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(lncRNA_ref_gencode[,mm10col],file=paste0(outdir,"/annotation/lncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   rnames_ncRNA=c('yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','miRNA','rRNA_gencode','rRNA_rmsk','lncRNA')
#   cnames_ncRNA=c('source','contents','Description')
# }
  
  
### fix table has NA's 
annotation_output_ncRNA=annotation_output[rnames_ncRNA,cnames_ncRNA]
write.table(annotation_output,file=paste0(annodir,"Annotations.txt"), sep = "\t", row.names = T, col.names = T, append = F, quote= FALSE)
write.table(annotation_output_ncRNA,file=paste0(annodir,"ncRNA_Annotations.txt"), sep = "\t", row.names = T, col.names = T, append = F, quote= FALSE)


#############################################################################
## Identify CLIP peak Gene Location

### Identify CLIP peak Host gene   

# Using GTF file from GENCODE v24  
# Peaks were annotated with overlapping gene  

#### FUNCTIONS
### fix files name should be through params
#################################################################################################
bam_anno2=function(peaksTable,Annotable,ColumnName){
  print ("bam starting")
  p=peaksTable[,c('chr','start','end','ID','ID2','strand')]
  a=Annotable[,c('chr','start','end','chr','chr','strand',ColumnName)]
  a[,4:5]=""
  colnames(a)[colnames(a)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')
  
  write.table(a,file=paste0(outdir,"/",misc,"/annotable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  write.table(p,file=paste0(outdir,"/",misc,"/peakstable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  ### fix bedtools should run be called from snakemake, not within R
  #remove lines with illegal characters
  system(paste0("cat ",outdir,misc,"/peakstable.bed | awk '$2 !~ /e/' > ",outdir,misc,"/peakstable.bed"))
  
  system(paste0('bedtools intersect -a ',gsub(" ","\\\\ ",outdir),misc,'/peakstable.bed -b ',gsub(" ","\\\\ ",outdir),misc,'/annotable.bed -wao -s > ',gsub(" ","\\\\ ",outdir),misc,'/peaks_mm10_OL.txt'))
  ab_OL=fread(paste0(outdir,misc,"/peaks_mm10_OL.txt"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(ab_OL)=c(paste0(colnames(p)),paste0(colnames(a)),'ntOL')
  
  ab_OL$width=ab_OL$end-ab_OL$start
  ab_OL$width_anno=ab_OL$end_anno-ab_OL$start_anno
  
  ab_OL=ab_OL[ab_OL$ntOL>0,]
  ab_OL$width_anno=ab_OL$end_anno-ab_OL$start_anno
  ab_OL$OLper_anno=ab_OL$ntOL/ab_OL$width_anno
  ab_OL$width=ab_OL$end-ab_OL$start
  ab_OL$OLper=ab_OL$ntOL/ab_OL$width
  ab_OL=ab_OL[(ab_OL$OLper_anno>.75 | ab_OL$OLper>.51), ]
  dup=unique(ab_OL[duplicated(ab_OL$ID),'ID'])
  ab_OL_single=ab_OL[!(ab_OL$ID%in%dup),]
  ab_OL_double=ab_OL[(ab_OL$ID%in%dup),]
  
  u=unique(ab_OL_double$ID)
  ab_OL_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(ab_OL_double)));
  colnames(ab_OL_colapsed)=colnames(ab_OL_double)
  
  for(x in 1:length(u)){
    p=u[x]
    pam=ab_OL_double[ab_OL_double$ID%in%p,]
    ab_OL_colapsed[x,colnames(ab_OL_double)%in%c(ColumnName)==F]=(pam[1,colnames(ab_OL_double)%in%c(ColumnName)==F,drop=T])
    
    for (cx in 1:length(ColumnName)) {
      ab_OL_colapsed[x,ColumnName[cx]]=paste(sort(unique((pam[,ColumnName[cx]]))), collapse =",")
    }
  }
  
  ab_OL_colapsed=rbind(ab_OL_colapsed,ab_OL_single)
  ab_OL_colapsed=ab_OL_colapsed[!is.na(ab_OL_colapsed$ID),]
  ab_OL_colapsed[(ab_OL_colapsed$ntOL<=0),ColumnName]=NA
  ab_OL=ab_OL_colapsed
  
  ### fix R shouldn't be removing files
  #system(paste0('rm ',gsub(" ","\\\\ ",outdir),'/',misc,'/annotable.bed'))
  #system(paste0('rm ',gsub(" ","\\\\ ",outdir),'/',misc,'/peakstable.bed'))
  #system(paste0('rm ',gsub(" ","\\\\ ",outdir),'/',misc,'/peaks_mm10_OL.txt'))
  
  return(ab_OL) 
}

rpmsk_anno=function(ColumnName,Annotable,peaksTable){
  s=peaksTable
  s=separate(s,ID,into=c('chr','start'),sep=":",remove=F)
  s=separate(s,start,into=c('start','end'),sep="-",remove=F)
  
  s.GR <- GRanges(seqnames = as.character(s$chr), ranges=IRanges(start = as.numeric(s$start), end = as.numeric(s$end)),strand = s$strand,ID=s$ID )
  
  anno.GR <- GRanges(seqnames = as.character(Annotable$genoName), ranges=IRanges(start = as.numeric(Annotable$genoStart), end = as.numeric(Annotable$genoEnd)),strand = Annotable$strand,repClass=Annotable$repClass,repName=Annotable$repName )
  
  q =s.GR
  s=anno.GR
  xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))
  
  qh=as.data.frame(q[xo$queryHits],row.names = NULL)
  sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_repeat")
  rmskinfo=cbind(qh,sh)
  rmskinfo=rmskinfo[,c('ID','repClass_repeat')];colnames(rmskinfo)[colnames(rmskinfo)%in%'repClass_repeat']= ColumnName      
  
  peaksTable=merge(peaksTable[,!colnames(peaksTable)%in%ColumnName],rmskinfo,by='ID',all.x=T)
  
  dup=unique(peaksTable[duplicated(peaksTable$ID),'ID'])
  peaksTable_single=peaksTable[!(peaksTable$ID%in%dup),]
  peaksTable_double=peaksTable[(peaksTable$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  for(x in 1:length(u)){
    p=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F]=(pam[1,colnames(peaksTable_double)%in%c(ColumnName)==F,drop=T])
    
    peaksTable_colapsed[x,ColumnName]=paste(sort(unique((pam[,ColumnName]))), collapse =" | ")
  }
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  peaksTable_colapsed[is.na(peaksTable_colapsed[,ColumnName]),ColumnName]=NA
  peaksTable_colapsed[peaksTable_colapsed[,ColumnName]%in%"",ColumnName]=NA
  peaksTable=peaksTable_colapsed
  
  return(peaksTable) 
}

### fix this should be a separate script
### fix designations - 1 = "same" 2 = "oppo"? same strand and complementary strand?
##################################################################################################
##################################################################################################
# ANNOTATE
##################################################################################################
for (xopp in 1:2) {
  if (xopp==1) {
    peaks=peaks;nmeprfix='Same_'
  } else if (xopp==2) {
    peaks=peaks_oppo;nmeprfix='Oppo_'
  }
  
  if (species=='hg38'){
    p=bam_anno2(peaks,
                rbind(mm10[,c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')],
                      ref_refseq[,c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')],### added anno for chrm not in gencode
                      lncRNA_ref_gencode[,c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')]
                ), c('ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL'))
  } else if (species=='mm10'){
    p=bam_anno2(peaks,
                rbind(mm10[,c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')]),c('ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL'))
  }
  PeaksdataOut = merge(peaks,
                       p[,c('ID','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')],
                       by='ID',all.x=T)
  
  ### Add Column that commbines all additional annotations
  annocol=c('chr','start','end','strand','type','name')
  if (species=='mm10'){
    Anno_RNA_comb=rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],
                        tRNA_sy[,annocol],sncRNA_sy[,annocol],rRNA_BK00964[,annocol],
                        SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],rRNA_rmsk[,annocol])
  } else if (species=='hg38'){
    Anno_RNA_comb=rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],
                        SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],
                        rRNA_rmsk[,annocol],tRNA_rmsk[,annocol])
  }
  
  # PeaksdataOut=bam_anno('RNA_anno',Anno_RNA_comb,PeaksdataOut)
  ### fix variable naming is confusing - back and forth between p and PeaksdataOut
  p=bam_anno2(PeaksdataOut,
              Anno_RNA_comb,
              c('type','name'))
  PeaksdataOut=merge(PeaksdataOut,p[,c('ID','type','name')],by='ID',all.x=T)
  
  
  ##################################################################################################
  #### CLEAN UP RNA TYPE
  ##################################################################################################
  
  p=PeaksdataOut
  
  ##########################################################################################
  ## change lincRNA,rRNA to rRNA only
  ## do not change the name if there is actually a lincRNA name separate from the rRNA grepl(',',p$name)==F)
  p[(grepl(',',p$name)==F)&(p$type%in%'lincRNA-exon,rRNA'),'type']='rRNA'
  p[(grepl(',',p$name)==F)&(p$type%in%'lincRNA-intron,rRNA'),'type']='rRNA'
  p[(grepl(',',p$name)==F)&(p$type%in%'lncRNA,rRNA'),'type']='rRNA'
  
  ##########################################################################################
  ## change all psueudogene classes to pseudogene
  p[grep('pseudogene',p$gene_type_ALL),'gene_type_ALL']='pseudogene'
  
  ##########################################################################################
  ## change all sc RNA to yRNA names
  p$type=gsub('scRNA,yRNA','yRNA',p$type)
  
  PeaksdataOut=p
  
  #########################
  ## Combine peak RNAtype info
  ########################
  p=PeaksdataOut
  
  p$type_simple_comb=NA
  p$type_comb=NA
  p$gene_name_comb=NA
  
  #######################################################
  dup=p[(is.na(p$gene_type_ALL)==F)|(is.na(p$type)==F),'ID']
  peaksTable_single=p[!(p$ID%in%dup),]
  peaksTable_double=p[(p$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  
  for(x in 1:length(u)){
    p=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    
    
    ###############################################
    ###############################################
    ### Comb Annotation
    
    pam_1=pam$gene_type_ALL ### type from Gencode 
    pam_1=as.data.frame(strsplit(pam_1,','));colnames(pam_1)='a';pam_1$a=as.character(pam_1$a)
    pam_2=pam$type ### Bitype from Annotation
    pam_2=as.data.frame(strsplit(pam_2,','));colnames(pam_2)='a';pam_2$a=as.character(pam_2$a)
    
    ##################################    
    
    ### remove misc_RNA catagory : these are covered by additional annotations (mostly yRNA)
    if (grepl('misc_RNA',pam_1)) {pam_1[pam_1$a%in%'misc_RNA','a']=NA}
    
    ### lincRNA are from an older Gencode version (VM18 or older) so don't use    
    if (grepl('lncRNA',pam_1)&grepl('lincRNA',pam_2)) {pam_1[pam_2$a%in%'lincRNA',]='lncRNA'}
    
    ### change ribozyme (gencode) to RNA type from additional anno
    if (grepl('ribozyme',pam_1)) {pam_1[pam_1$a%in%'ribozyme',]=unique(pam_2$a)}
    
    ##################################    
    ## combine all rna types subtypes into ncRNA
    pam_c=rbind(pam_1,pam_2)
    pam_c=pam_c[!is.na(pam_c$a),,drop=F]
    
    #### Annotation from Gencode : unique(sort(mm10$gene_type_ALL))
    ### fix with loop
    pam_c2=pam_c
    pam_c2$a=gsub('miRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('miscRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('misc_RNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('piRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('rRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('siRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('snRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('snoRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('ribozyme','ncRNA',pam_c2$a)
    pam_c2$a=gsub('scRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('sRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('scaRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('rRNA','ncRNA',pam_c2$a)
    
    #### extra annotation unique(Anno_RNA_comb$type )
    pam_c2$a=gsub('yRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('srpRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('tRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('7SKRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('scRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('sncRNA','ncRNA',pam_c2$a)
    pam_c2$a=gsub('rRNA','ncRNA',pam_c2$a)
    
    if (length(pam_c$a)>1) {pam_c=pam_c[order(pam_c$a),,drop=F]}
    if (length(pam_c2$a)>1) {pam_c2=pam_c2[order(pam_c2$a),,drop=F]}
    
    pam_c=paste(unique(pam_c$a),collapse = ',')
    pam_c2=paste(unique(pam_c2$a),collapse = ',')
    
    
    ########################################################################################################################################
    ########################################################################################################################################
    ### Comb GeneName 
    
    pam_G1=pam$external_gene_name ### gene name from gencode
    pam_G1=as.data.frame(strsplit(pam_G1,','));colnames(pam_G1)='a';pam_G1$a=as.character(pam_G1$a)
    pam_G2=pam$name ### gene name from annotation
    pam_G2=as.data.frame(strsplit(pam_G2,','));colnames(pam_G2)='a';pam_G2$a=as.character(pam_G2$a)
    
    if (grepl('sk',pam_G1)&grepl('7SK',pam_G2)) {
      pam_G2[pam_G2$a%in%'7SK',]=NA
    }
    
    pam_G3=rbind(pam_G1,pam_G2)
    pam_G3=pam_G3[!is.na(pam_G3$a),,drop=F]
    pam_G3=paste(unique(pam_G3$a),collapse = ',')
    
    
    ########################################################################################################################################
    ########################################################################################################################################
    ### Create Table 
    
    peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
    peaksTable_colapsed[x,'type_simple_comb']=pam_c2
    peaksTable_colapsed[x,'type_comb']=pam_c
    peaksTable_colapsed[x,'gene_name_comb']=pam_G3
    
    remove('pam_1','pam_2','pam_c','pam_c2','pam_G1','pam_G2','pam_G3','pam')
  }
  
  p=peaksTable_colapsed
  rnatype=p[,c('ID','ensembl_gene_id','external_gene_name','gene_type',"gene_type_ALL",'type','name','type_simple_comb','type_comb','gene_name_comb')]
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  colnames(peaksTable_colapsed)[colnames(peaksTable_colapsed)%in%c("ensembl_gene_id","external_gene_name","gene_type","gene_type_ALL","type","name","type_simple_comb","type_comb","gene_name_comb")]=
    paste0(nmeprfix,c("ensembl_gene_id","external_gene_name","gene_type","gene_type_ALL","type","name","type_simple_comb","type_comb","gene_name_comb"))
  PeaksdataOut=peaksTable_colapsed
  
  # check lincRNA chr1:156646852-156647108
  ##########################################################################################
  ## change lincRNA,RNA to RNA only
  ## do not change the name if there is actually a lincRNA name separate from the rRNA grepl(',',p$name)==F)
  ##########################################################################################
  
  
  #############################################################################################################
  #############################################################################################################
  ### INTRON EXON ANNOTATION
  
  ### Identify if CLIP peak overlaps with Intron or Exonic region   
  
  #   Using GTF file from GENCODE v23 -mm10
  #   Using GTF file from GENCODE v32 -hg38
  
  # Peaks were annotated by whether they overlap with Host gene intron/exon region
  # Intron coordinates were calculated from GTF file.
  # 
  # A second column was added to idenify if the peak also overlapped with the 5'UTR 3'UTR or CDS (Column: Featrue 2)
  
  #############################################################################################################
  #############################################################################################################
  
  Annotable=intron_exon[grep('protein_coding',intron_exon$gene_type),]
  peaksTable=PeaksdataOut
  ColumnName=paste0(c('ensembl_gene_id','external_gene_name'))
  ColumnName=c("feature","exon_number")
  
  p=peaksTable[,c('chr','start','end','ID','ID2','strand')]
  a=Annotable[,c('chr','start','end','transcript_id','transcript_id','strand',ColumnName)]
  a[,5]=""
  colnames(a)[colnames(a)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')
  
  ### fix - why are we writing the same files over again? 
  write.table(a,file=paste0(outdir,"/",misc,"/annotable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  write.table(p,file=paste0(outdir,"/",misc,"/peakstable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  system(paste0('bedtools intersect -a ',gsub(" ","\\\\ ",outdir),'/',misc,'/peakstable.bed -b ',gsub(" ","\\\\ ",outdir),'/',misc,'/annotable.bed -wao -s  >',gsub(" ","\\\\ ",outdir),'/',misc,'/peaks_mm10_OL.txt'))
  
  ### fix - writing same files again
  exoninof=fread(paste0(outdir,"/",misc,"/peaks_mm10_OL.txt"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(exoninof)=c(paste0(colnames(p)),paste0(colnames(a)),'ntOL')
  exoninof=exoninof[exoninof$ntOL>0,]
  exoninof$width_anno=exoninof$end_anno-exoninof$start_anno
  exoninof$OLper_anno=exoninof$ntOL/exoninof$width_anno
  exoninof$width=exoninof$end-exoninof$start
  exoninof$OLper=exoninof$ntOL/exoninof$width
  
  exoninof=exoninof[(exoninof$OLper_anno>.75 | exoninof$OLper>.51), ]
  
  peaksTable[,paste0(nmeprfix,'feature')]=NA
  peaksTable[,paste0(nmeprfix,'exon_number')]=NA
  peaksTable[,paste0(nmeprfix,'exon_LargeOL')]=NA
  peaksTable[,paste0(nmeprfix,'intron_number')]=NA
  peaksTable[,paste0(nmeprfix,'intron_LargeOL')]=NA
  peaksTable[,paste0(nmeprfix,'intron_5pStart')]=NA
  peaksTable[,paste0(nmeprfix,'intron_length')]=NA
  peaksTable[,paste0(nmeprfix,'exon_length')]=NA
  
  if (xopp==1) {
    ###########
    #### Calc distance of 5' peak to 5' intron/exondist
    exoninof_pos=exoninof[exoninof$strand%in%'+',]
    exoninof_neg=exoninof[exoninof$strand%in%'-',]
    
    exoninof_pos[,paste0(nmeprfix,'feature_Distance')]=((exoninof_pos$start-exoninof_pos$start_anno)/abs(exoninof_pos$end_anno-exoninof_pos$start_anno))*100
    exoninof_pos[,paste0(nmeprfix,'feature_5pStart')]=exoninof_pos$start_anno
    exoninof_pos[,paste0(nmeprfix,'feature_length')]=abs(exoninof_pos$end_anno-exoninof_pos$start_anno)
    
    exoninof_neg[,paste0(nmeprfix,'feature_Distance')]=((exoninof_neg$end_anno-exoninof_neg$end)/abs(exoninof_neg$end_anno-exoninof_neg$start_anno))*100
    exoninof_neg[,paste0(nmeprfix,'feature_5pStart')]=exoninof_neg$end_anno
    exoninof_neg[,paste0(nmeprfix,'feature_length')]=abs(exoninof_neg$end_anno-exoninof_neg$start_anno)
    
    
    exoninof=rbind(exoninof_pos,exoninof_neg);
    
    peaksTable[,paste0(nmeprfix,'Exn_start_dist')]=NA
    peaksTable[,paste0(nmeprfix,'Intron_start_dist')]=NA
    peaksTable[,paste0(nmeprfix,'Intron_5pStart')]=NA
    peaksTable[,paste0(nmeprfix,'Exn_5pStart')]=NA
    
  }
  # xxx1
  ########### 
  # ID="chr19:5490528-5490574"
  for (x in 1:nrow(peaksTable)) {
    l=peaksTable[x,c('ID','ID2')]
    g=exoninof[exoninof$ID%in%l$ID,]
    
    if (nrow(g)>0) {
      g_e=g[g$feature%in%'exon',]
      g_i=g[g$feature%in%'intron',]
      
      
      gname=paste(unique(g$feature),collapse = ",")
      gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
      peaksTable[x,paste0(nmeprfix,'feature')]=gsub('NA,',"",gname)
      
      if (nrow(g_e)>0) {
        gname=paste(unique(g_e$exon_number),collapse = ",")
        gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
        peaksTable[x,paste0(nmeprfix,'exon_number')]=gsub('NA,',"",gname)
        
        peaksTable[x,paste0(nmeprfix,'exon_LargeOL')]=max(g_e$ntOL/(g_e$end-g_e$start))
        
        if (xopp==1) {
          g_e[g_e[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
          peaksTable[x,paste0(nmeprfix,'Exn_start_dist')]=mean(as.numeric(g_e[,paste0(nmeprfix,'feature_Distance')]))
        }
      }
      
      if (nrow(g_i)>0) {
        gname=paste(unique(g_i$exon_number),collapse = ",")
        gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
        peaksTable[x,paste0(nmeprfix,'intron_number')]=gsub('NA,',"",gname)
        
        peaksTable[x,paste0(nmeprfix,'intron_LargeOL')]=max(g_i$ntOL/(g_i$end-g_i$start))
        
        if (xopp==1) {
          g_i[g_i[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
          peaksTable[x,paste0(nmeprfix,'Intron_start_dist')]=mean(as.numeric(g_i[,paste0(nmeprfix,'feature_Distance')]))
        }
      }
    }
    if (nrow(g)==0) {
      if (xopp==1) {peaksTable[x,paste0(nmeprfix,c("feature",'exon_number','intron_number','Exn_start_dist','Intron_start_dist'))]=NA}
      if (xopp==2) {peaksTable[x,paste0(nmeprfix,c("feature",'exon_number','intron_number'))]=NA}
    } 
  }
  
  ############################## to fix this Get introns from refseqNCBI
  peaksTable[is.na(peaksTable[,paste0(nmeprfix,'feature')])&peaksTable[,paste0(nmeprfix,'gene_type')]%in%'protein_coding',paste0(nmeprfix,'feature')]='exon'
  ##############################
  
  peaksTable_inex=peaksTable[grep(',',peaksTable[,paste0(nmeprfix,'feature')]),]
  
  PeaksdataOut=peaksTable
  #############################################################################################################
  #############################################################################################################
  ### IDENTIFY PEAKS IN REPEAT REGIONS   
  
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
  #############################################################################################################
  #############################################################################################################
  
  PeaksdataOut=PeaksdataOut[!is.na(PeaksdataOut$ID),]
  
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LINE_SINE'),rmsk_GRCm38_LISI,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LTR'),rmsk_GRCm38_LTR,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_DNA'),rmsk_GRCm38_DNA,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Satalites'),rmsk_GRCm38_sat,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Simple_Repeats'),rmsk_GRCm38_SR,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Low_Complexity'),rmsk_GRCm38_LowComplx,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Other'),rmsk_GRCm38_Other,PeaksdataOut)
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Unknown'),rmsk_GRCm38_unknown,PeaksdataOut)
  
  
  p=PeaksdataOut
  p[,paste0(nmeprfix,'Repeat_comb')]=NA
  repcol=paste0(nmeprfix,c('Repeat_LINE_SINE','Repeat_LTR','Repeat_DNA','Repeat_Satalites','Repeat_Simple_Repeats','Repeat_Low_Complexity','Repeat_Other','Repeat_Unknown'))
  
  dup=p[rowSums((p[,repcol]>0),na.rm = T)>0,'ID']
  
  peaksTable_single=p[!(p$ID%in%dup),]
  peaksTable_double=p[(p$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  
  for(x in 1:length(u)){
    p=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    pam_c=((pam[1,repcol,drop=F]))
    
    pmat=as.data.frame(matrix(nrow=length(pam_c),ncol=1));colnames(pmat)='a'
    pmat$a=t(pam_c)
    pmat=pmat[is.na(pmat)==F,]
    pmat=paste(unique(pmat),collapse = ',')
    
    ######################     
    
    peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
    peaksTable_colapsed[x,paste0(nmeprfix,'Repeat_comb')]=pmat
    
    remove('p','pam','pam_c')
  }
  
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  p=peaksTable_colapsed
  PeaksdataOut=peaksTable_colapsed
  
  
  #############################################################################################################
  #############################################################################################################
  ### Asigning Clip peak attributes   
  
  # Not all Peaks overlap with a single feature so peak assignments were assigned by priority:  
  # 
  # ncRNA > Protein coding : Exonic > repeats > Pseudogene > Antisense Feature > Protein Coding : Intronic > lncRNA > no Feature  
  # 
  # All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table.   
  
  #############################################################################################################
  #############################################################################################################
  ### fix with loop
  PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]=NA
  PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')] )
  PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')] )
  PeaksdataOut[,paste0(nmeprfix,'type_comb')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_comb')] )
  PeaksdataOut[,paste0(nmeprfix,'type_comb')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_comb')] )
  PeaksdataOut[,paste0(nmeprfix,'gene_type')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type')] )
  PeaksdataOut[,paste0(nmeprfix,'gene_type')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type')] )
  PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')] )
  PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')] )   
  
  # 1. ncRNA        
  comp=( grepl('ncRNA',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) ) 
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='ncRNA'
  
  # 2. protein coding - Exonic
  comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')])&(grepl('exon',PeaksdataOut[,paste0(nmeprfix,'feature')])) )
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='protein_coding: exon'
  
  # 3. repeats
  comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& (is.na(PeaksdataOut[,paste0(nmeprfix,'Repeat_comb')])==F) )
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]="Repeat Element"
  
  # 4. Pseudogene
  comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& grepl('pseudogene',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='pseudogene'
  
  # 5. lncRNA-exon
  comp1=is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])  
  comp2=grepl('linLcRNA-exon',PeaksdataOut[,paste0(nmeprfix,'type_comb')]) | grepl('lnLcRNA',PeaksdataOut[,paste0(nmeprfix,'type_comb')])
  comp=comp1&comp2       
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='lnLcRNA-exon'
  
  # 6. intron
  comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='protein_coding: Intron'
  
  # 7. lncRNA-intron
  comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]) & (grepl('linLcRNA-intron',PeaksdataOut[,paste0(nmeprfix,'type_comb')]) ))
  PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='lnLcRNA-intron'
  PeaksdataOut[is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]),paste0(nmeprfix,'Comb_type_exon')]='no Feature'
  PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]=factor(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lnLcRNA-exon","Antisense Feature","protein_coding: Intron","lnLcRNA-intron","no Feature"))
  
  ##############################
  #### RNA subtypes
  ##############################
  
  p=PeaksdataOut
  p[,paste0(nmeprfix,'Comb_type_ncRNA')]=NA
  p1=p[p[,paste0(nmeprfix,'Comb_type_exon')]%in%'ncRNA',]
  p2=p[!p[,paste0(nmeprfix,'Comb_type_exon')]%in%'ncRNA',]
  
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=p1[,paste0(nmeprfix,'type_comb')]
  
  ### Protein coding + ncRNA annotations -> ncRNA subtype
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('protein_coding,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(',protein_coding',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  
  ### any ncRNA takes priority over pseudogene
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('pseudogene,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(',pseudogene',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  
  ### any double annotations with lncRNA become second annotation only
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('lnLcRNA,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[p1[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-intron,linLcRNA-exon',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
  p1[p1[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-exon,linLcRNA-intron',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
  
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA-exon,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA-intron,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  
  ### fix with loop
  ### any double annotations with yRNA become yRNA only
  p1[grep('yRNA',p1[,paste0(nmeprfix,'Comb_type_ncRNA')]),paste0(nmeprfix,'Comb_type_ncRNA')]='yRNA'
  
  ### tRNA takes priority over miRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,tRNA',"miRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  ### miRNA takes priority over rRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,rRNA',"miRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  ### snRNA takes priority over 7SKRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,snRNA',"snRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('snRNA,7SKRNA',"snRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  ### snRNA takes priority over 7SKRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,snRNA',"snRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  ### tRNA takes priority over rRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('rRNA,tRNA',"tRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('tRNA,rRNA,',"tRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  ### snoRNA takes priority over miRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,snoRNA',"snoRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('snoRNA,miRNA,',"snoRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  ### scaRNA takes priority over miRNA
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,scaRNA',"scaRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('scaRNA,miRNA,',"scaRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  PeaksdataOut=rbind(p1,p2)
  
  if (xopp==1) {PeaksdataOut_same=PeaksdataOut}
  if (xopp==2) {PeaksdataOut=merge(PeaksdataOut_same,PeaksdataOut[,colnames(PeaksdataOut)[!colnames(PeaksdataOut)%in%c("chr","start","end","strand","ID2" )]],by='ID')}
  
}### for pos or neg anno


PeaksdataOut$Comb_type_exon_Oppo=NA
# PeaksdataOut$type_simple_comb

# 1. ncRNA        
comp=(grepl('ncRNA',PeaksdataOut[,paste0('Same_','type_simple_comb')],fixed = T) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']='ncRNA'

# 2. protein coding - Exonic
comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])&grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')])&(grepl('exon',PeaksdataOut[,paste0('Same_','feature')])) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']=paste0('protein_coding: exon')

# 3. repeats
comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])&(is.na(PeaksdataOut[,paste0('Same_','Repeat_comb')])==F) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']="Repeat Element"

# 4. Pseudogene
comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& grepl('pseudogene',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
PeaksdataOut[comp,'Comb_type_exon_Oppo']='pseudogene'

# 5. Antisense - non LncRNA
comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'

# 6. lncRNA-exon
comp1=is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])
comp2=((PeaksdataOut[,paste0('Same_','type_comb')]%in%c('linLcRNA-exon')) | (PeaksdataOut[,paste0('Same_','type_comb')]%in%c('lnLcRNA')) )
comp=comp1&comp2       
PeaksdataOut[comp,'Comb_type_exon_Oppo']='lnLcRNA-exon'

# 7. Antisense - LncRNA - exon
comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA-intron',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'

# 8. intron
comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
PeaksdataOut[comp,'Comb_type_exon_Oppo']=paste0('protein_coding: Intron')

# 9. Antisense - LncRNA - intron
comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA-exon',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'

# 10. lncRNA -intron
comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo']) & ((PeaksdataOut[,paste0('Same_','type_comb')]%in%'linLcRNA-intron')==T) )
PeaksdataOut[comp,'Comb_type_exon_Oppo']='lnLcRNA-intron'
PeaksdataOut[is.na(PeaksdataOut[,'Comb_type_exon_Oppo']),'Comb_type_exon_Oppo']='no Feature'


### only distances for final annotation_output
## distances to any protein coading annotation
PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Exn_start_dist")]=NA
PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Intron_start_dist")]=NA

################################################  
### Revert lnLcRNA to lncRNA  

PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
PeaksdataOut[,'Comb_type_exon_Oppo']=factor(PeaksdataOut[,'Comb_type_exon_Oppo'], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))

unlink(misc, recursive = TRUE)  
return(PeaksdataOut)
}


flag = "run"

if (flag=="R"){
  Peaksdata2=read.csv("/Volumes/sevillas2/git/iCLIP/workflow/scripts/peakstest.csv")
  peaks = Peaksdata2[,c('chr','start','end','strand')]
  WriteClassTable = T
  species = "hg38"
  outdir = "/Volumes/data/iCLIP/marco/14_annotation/"
  Ref = "/Volumes/iCLIP/ref/CLIP_Anno"
  CLIPannotation(peaks,WriteClassTable,species,outdir,Ref)
} else if (flag == "local") {
  Peaksdata2=read.csv("/home/sevillas2/git/iCLIP/workflow/scripts/peakstest.csv")
  peaks = Peaksdata2[,c('chr','start','end','strand')]
  WriteClassTable = T
  species = "hg38"
  outdir = "/data/sevillas2/iCLIP/marco/14_annotation/"
  Ref = "/data/RBL_NCI/iCLIP/ref/CLIP_Anno"
  CLIPannotation(peaks,WriteClassTable,species,outdir,Ref)
}

