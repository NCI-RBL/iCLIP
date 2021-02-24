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
ref_species = args[1]
soyeong_flag = args[2]
soyeong_path = args[3]
alias_path = args[4]
gencode_path = args[5]
refseq_path = args[6]
canonical_path = args[7]
intron_path = args[8]
rmsk_path = args[9]
out_dir = args[10]
reftable_path = args[11]

#testing information
testing="N"
if(testing=="Y"){
  ref_dir = "/Volumes/RBL_NCI/iCLIP/ref/annotations/"
  ref_species = "hg38" ### need better name, match to snakemake

  ### fix need to figure out how to keep all this info - maybe a config? dict?
  alias_path = paste0(ref_dir,ref_species,"/",ref_species,".chromAlias.txt")
  if(ref_species == "mm10"){
    out_dir = "/Volumes/data/iCLIP/mm10/15_annotation/"
    gencode_path = paste0(ref_dir, "mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "/mm10/NCBI_RefSeq/GCF_000001635.26_GRCm38.p6_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"/mm10/Gencode_VM23/fromUCSC/KnownCanonical/KnownCanonical_GencodeM23_GRCm38.txt")
    intron_path = paste0(ref_dir, "/mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
    rmsk_path = paste0(ref_dir,"/mm10/repeatmasker/rmsk_GRCm38.txt")
    soyeong_flag = "Y" #Y or N
    reftable_path = "/Volumes/sevillas2/git/iCLIP/config/annotation_config_mm10.txt"
  } else if (ref_species == "hg38"){
    out_dir = "/Volumes/data/iCLIP/marco/15_annotation/"
    gencode_path = paste0(ref_dir,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"/hg38/Gencode_V32/fromUCSC/KnownCanonical/KnownCanonical_GencodeM32_GRCh38.txt")
    intron_path = paste0(ref_dir,"/hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
    rmsk_path = paste0(ref_dir,"/hg38/repeatmasker/rmsk_GRCh38.txt")
    soyeong_flag = "N" #always now (for now)
    reftable_path = "/Volumes/sevillas2/git/iCLIP/config/annotation_config_hg.txt"
  } 
}

##########################################################################################
############### Annotation info
##########################################################################################
#read in from reference files
alias_anno=fread(alias_path, header=T, sep="\t",
                 stringsAsFactors = F,data.table=F,fill=TRUE)

###phil why are we naming NCBI (CM###) when it's GENCODE
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

write.csv(alias_anno,paste0(out_dir,"ref_alias.csv"))
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

write.csv(ref_gencode_t,paste0(out_dir,"ref_gencode.csv"))

##########################################################################################
############### REFSEQ ANNOTATION
##########################################################################################
###phil the ref being read has NA in all cols past protein_in so data just gets removd with 
#filtering 
#head /data/RBL_NCI/iCLIP/ref/CLIP_Anno/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt 
ref_refseq = fread(refseq_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

###phil why do we have dup col names (source)
#handle duplicate col name in hg38
if(ref_species=="hg38"){
  colnames(ref_refseq)[34] = "source2"
}

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
###phil - there is  no "transcript_type col in refseq
#ref_refseq=ref_refseq[(ref_refseq$gene_type%in%'TEC')==F,]  
ref_refseq = subset(ref_refseq, gene_type == "TEC")   
#ref_refseq = subset(ref_refseq, transcript_type == "TEC") 

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
write.csv(ref_refseq,paste0(out_dir,"ref_refseq.csv"))


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
###phil this file doesn't exist
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
##########################################################################################
############### Annotation functions
##########################################################################################
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
  write.csv(lncRNA_ref_gencode,paste0(out_dir,"ref_lncRNA.csv"))
  
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
  
  write.table(df_sub, file=paste0(out_dir,rowid,".bed"), 
              sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
  
  
  return(list(df_sub,content))
}

rmskAnno<-function(rowid){
  
  #define family and class lists
  family_list = c("rRNA","snRNA", "srpRNA","tRNA")
  
  #remove other from hg38 - not found in rmsk_GRCm38
  if(ref_species=="hg38"){
    class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','Unknown')
  }else{
    class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','Unknown','Other')
  }
  
  #if annotation id (rowid) is in the defined annotaiton list, subset, and outuput bed
  annotation_list = c(family_list,class_list,"7SK RNA","LINE SINE","rRNA_DNA",
                      "scRNA", "yRNA")
  
  if(rowid %in% annotation_list){
    #define cols to select and new col names
    rmskcol=c('genoName','genoStart','genoEnd','repName','swScore','strand')
    newcolnames=c('chr','start','end','transcript_name','score','strand')
    
    #if in the family or class list
    if(rowid %in% family_list){
      df_sub = subset(rmsk_GRCm38, repFamily == rowid) %>%
        select(rmskcol) %>%
        setnames(newcolnames)
    } else if (rowid %in% class_list){
      df_sub = subset(rmsk_GRCm38, repClass == rowid)
      #Special handling  
    } else if(rowid== "7SK RNA"){
      df_sub = subset(rmsk_GRCm38, repFamily == 'RNA') %>%
        select(rmskcol) %>%
        setnames(newcolnames)
      rowid="7SKRNA"
    } else if (rowid == "LINE SINE"){
      df_sub = subset(rmsk_GRCm38, repClass == c("LINE","SINE"))
    } else if (rowid == "rRNA_DNA"){
      df_sub = subset(rmsk_GRCm38, repClass == "DNA")
    } else if (rowid == "scRNA"){
      df_sub = subset(rmsk_GRCm38, repFamily == 'scRNA' & grepl("HY",rmsk_GRCm38$repName)==FALSE) %>%
        select(rmskcol) %>%
        setnames(newcolnames)
    } else if (rowid == "yRNA"){
      df_sub = subset(rmsk_GRCm38, repFamily == 'scRNA' & grepl("HY",rmsk_GRCm38$repName)==TRUE) %>%
        select(rmskcol) %>%
        setnames(newcolnames)
    } 
    df_sub$type=rowid
    write.table(df_sub,file=paste0(out_dir,rowid,".bed"), 
                sep = "\t", row.names = F, append = F, quote= FALSE)
    
    content = paste0(unique(df_sub$name),collapse = ", ")
  } else{
    df_sub=data.frame()
    content=""
  }
  
  return(list(df_sub,content))
}
   
SYAnno<-function(rowid,ref_species){
  #read bedfile
  file_name = ref_table[rowid,"SY_1"]
  df_sub = read.table(paste0(soyeong_path,file_name))
  
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
              file=paste0(out_dir,rowid,".bed"), 
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
###phil need to fix
if (ref_species=='mm10') {
  annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12 and various predicted genes'
  annotation_output['rRNA_gencode','contents']='5S, 5.8s, predicted gene'
} else if (ref_species=='hg38') {
  annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12'
  annotation_output['rRNA_gencode','contents']='5S, 5.8s, predicted gene'
}

#write out annotation files  
rnames_ncRNA = c('yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','miRNA','rRNA','lncRNA')
cnames_ncRNA = c('source','contents','description')
write.table(annotation_output,
            file=paste0(out_dir,"annotations.txt"), 
            sep = "\t", row.names = T, col.names = T, append = F, quote= FALSE)
write.table(annotation_output[rnames_ncRNA,cnames_ncRNA],
            file=paste0(out_dir,"ncRNA_annotations.txt"), 
            sep = "\t", row.names = T, col.names = T, append = F, quote= FALSE)

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
#   write.table(YRNA_rmsk[,rmskcol],file=paste0(out_dir,"yRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(srpRNA_rmsk[,rmskcol],file=paste0(out_dir,"srpRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(SKRNA_rmsk[,rmskcol],file=paste0(out_dir,"SKRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(scRNA_rmsk[,rmskcol],file=paste0(out_dir,"scRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_rmsk[,rmskcol],file=paste0(out_dir,"rRNA_rmsk.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_BK00964[,rmskcol],file=paste0(out_dir,"rRNA_BK00964.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   
#   write.table(snRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"snRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(snoRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"snoRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(miRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"mirNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"rRNA_gencode.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(lncRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"lncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   write.table(tRNA_sy[,rmskcol],file=paste0(out_dir,"tRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(sncRNA_sy[,c('chr','start','end','name','V6','strand')],file=paste0(out_dir,"sncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   rnames_ncRNA=c('yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','sncRNA','miRNA','rRNA_gencode','rRNA_rmsk','rRNA_DNA','lncRNA','lincRNA')
#   cnames_ncRNA=c('source','contents','Description') 
# } else if (species=='hg38'){
#   write.table(YRNA_rmsk[,rmskcol],file=paste0(out_dir,"yRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(srpRNA_rmsk[,rmskcol],file=paste0(out_dir,"srpRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(tRNA_rmsk[,rmskcol],file=paste0(out_dir,"tRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(SKRNA_rmsk[,rmskcol],file=paste0(out_dir,"SKRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(scRNA_rmsk[,rmskcol],file=paste0(out_dir,"scRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_rmsk[,rmskcol],file=paste0(out_dir,"rRNA_rmsk.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   write.table(snRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"snRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(snoRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"snoRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(miRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"mirNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(rRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"rRNA_gencode.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   write.table(lncRNA_ref_gencode[,newcolnames],file=paste0(out_dir,"lncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
#   
#   rnames_ncRNA=c('yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','miRNA','rRNA_gencode','rRNA_rmsk','lncRNA')
#   cnames_ncRNA=c('source','contents','Description')
# }
