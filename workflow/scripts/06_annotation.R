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
soyeong_path = args[2]
alias_path = args[3]
gencode_path = args[4]
refseq_path = args[5]
canonical_path = args[6]
intron_path = args[7]
rmsk_path = args[8]
out_dir = args[9]
reftable_path = args[10]

setwd("../../../../../")

#testing information
if(length(args)==0){
  ref_dir = "Volumes/RBL_NCI/iCLIP/ref/annotations/"
  ref_species = "hg38" ### need better name, match to snakemake
  reftable_path = "Volumes/sevillas2/git/iCLIP/config/annotation_config.txt"

  ### fix need to figure out how to keep all this info - maybe a config? dict?
  alias_path = paste0(ref_dir,ref_species,"/",ref_species,".chromAlias.txt")
  if(ref_species == "mm10"){
    out_dir = "Volumes/data/CCBR_Projects/iCLIP/testing/"
    gencode_path = paste0(ref_dir, "mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "mm10/NCBI_RefSeq/GCF_000001635.26_GRCm38.p6_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"mm10/Gencode_VM23/fromUCSC/KnownCanonical/KnownCanonical_GencodeM23_GRCm38.txt")
    intron_path = paste0(ref_dir, "mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
    rmsk_path = paste0(ref_dir,"mm10/repeatmasker/rmsk_GRCm38.txt")
    soyeong_path= paste0(ref_dir,'mm10/AdditionalAnno/')
    
  } else if (ref_species == "hg38"){
    out_dir = "Volumes/data/CCBR_Projects/iCLIP/testing/"
    gencode_path = paste0(ref_dir,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"hg38/Gencode_V32/fromUCSC/KnownCanonical/KnownCanonical_GencodeM32_GRCh38.txt")
    intron_path = paste0(ref_dir,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
    rmsk_path = paste0(ref_dir,"hg38/repeatmasker/rmsk_GRCh38.txt")
  } 
}

##########################################################################################
############### GENCODE ANNOTATION
##########################################################################################
ref_gencode = fread(gencode_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

#dplyr::rename cols
ref_gencode = ref_gencode %>% 
  dplyr::rename(
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

# #remove rows with missing data
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

# create a column $gene_type_ALL to keep track of ncRNA sub catagories. 
# when we plot CLIP peaks we plot both ncRNA as a single group and ncRNA as subgroups
ReplaceType <-function(df_in,col_in,list_in,replace_id){
  for (id in list_in){
    df_in[,col_in]=gsub(id,replace_id,df_in[,col_in])
  } 
  return(df_in)
}

type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA","tRNA","ribozyme")
ref_gencode = ReplaceType(ref_gencode,"gene_type", type_list,"ncRNA")

#create subsets used to annotate clip peaks separately 
ref_gencode_t = subset(ref_gencode, feature == 'transcript')
ref_gencode_e = subset( ref_gencode, feature == "exon")
write.csv(ref_gencode_t,paste0(out_dir,"ref_gencode.csv"))


##########################################################################################
############### canonical paths
##########################################################################################
canonical=fread(canonical_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

#remove version
canonical$transcript=removeVersion(canonical$transcript)
canonical$protein=removeVersion(canonical$protein)

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

# intron table is 0 based while exon starts at 1
# same for exon number the exon # starts at 0
introns$start=introns$start+1
introns$exon_number=as.numeric(introns$exon_number)+1

## add gene type type
introns=merge(introns,
              (unique(ref_gencode_t[,c('transcript_id','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL','transcript_type','transcript_name','score')])),by='transcript_id',all.x=T)

##########################################################################################
############### binding exon, intron
##########################################################################################
intron_exon=rbind(ref_gencode_e[,c('chr','feature','start','end','strand','transcript_id','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL','transcript_type','transcript_name','score','exon_number')],
                  introns[,c('chr','feature','start','end','strand','transcript_id','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL','transcript_type','transcript_name','score','exon_number')])
intron_exon$ID=paste0(intron_exon$chr,':',intron_exon$start,'-',intron_exon$end)

##########################################################################################
############### REPEATMASKER
##########################################################################################
# get repeat regions and extra small RNA locations
rmsk_GRCm38=fread(rmsk_path, header=T, sep="\t",stringsAsFactors = F,
                  data.table=F)


##########################################################################################
############### Annotation functions
##########################################################################################
gencodeLNCRNA<-function(){
  #identfy a subset of lncRNA that have intronic sequences.
  # we are then combining the Intronic and Exonic locations and changing the gene_type_ALL from lnc to linc
  # finally we will combine the linc + lnc (no introns) together into one lnc table
  
  #Get all lncRNA transcripts
  lncRNA_ref_gencode= subset(ref_gencode_t, gene_type_ALL == "lncRNA")
  
  #lncRNA with introns
  lincRNA_ref_gencode=lncRNA_ref_gencode[lncRNA_ref_gencode$transcript_id%in%unique(introns$transcript_id),]
  lincRNA_ref_gencode$gene_type_ALL="lincRNA"
  
  #lncRNA with out introns
  lncRNA_ref_gencode_noInt=lncRNA_ref_gencode[!lncRNA_ref_gencode$transcript_id%in%introns$transcript_id,]
  
  ##lincRNA exon coordinates
  lincRNA_ref_gencode_IntronExon=intron_exon[intron_exon$transcript_id%in%lincRNA_ref_gencode$transcript_id,]
  lincRNA_ref_gencode_IntronExon$gene_type_ALL="lincRNA"
  lincRNA_ref_gencode_IntronExon$gene_type_ALL=paste0(lincRNA_ref_gencode_IntronExon$gene_type_ALL,'-',lincRNA_ref_gencode_IntronExon$feature)
  
  
  lncRNA_ref_gencode = rbind(lncRNA_ref_gencode_noInt[,c('chr','start','end','strand',
                                                         'ensembl_gene_id','external_gene_name',
                                                         'gene_type','gene_type_ALL','feature',
                                                         'transcript_id','transcript_type',
                                                         'transcript_name','score')],
                             lincRNA_ref_gencode_IntronExon[,c('chr','start','end','strand',
                                                               'ensembl_gene_id','external_gene_name',
                                                               'gene_type','gene_type_ALL','feature',
                                                               'transcript_id','transcript_type',
                                                               'transcript_name','score')])
  
  # somewhat redefined lncRNA with linc and included exon and intron sequences - 
  # treating lncRNA as separate annotation table and removing these transcripts from ref_gencode_t
  ref_gencode_t=ref_gencode_t[!ref_gencode_t$transcript_id%in%unique(lncRNA_ref_gencode$transcript_id),]
  write.csv(lncRNA_ref_gencode,paste0(out_dir,"ref_lncRNA.csv"))
  
  return(lncRNA_ref_gencode)
}

gencodeAnno<-function(rowid){
  #introns
  df_sub = subset(ref_gencode_t, transcript_type == rowid)
  content = paste0(unique(df_sub$transcript_type),collapse = ', ')
  
  if(rowid=="intron"){
    df_sub = introns
    content = paste0(unique(df_sub$feature),collapse = ', ')
  } 
  if (rowid %in% c("lncRNA")){
    df_sub = gencodeLNCRNA() %>%
      subset(gene_type_ALL == "lncRNA")
    content = paste0(unique(df_sub$transcript_type),collapse = ', ')
  }
  if (rowid %in% c("lincRNA")){
    df_sub = gencodeLNCRNA() 
    df_sub = subset(df_sub,grepl("lincRNA",df_sub$gene_type_ALL))
    content = paste0(unique(df_sub$transcript_type),collapse = ', ')
  }
  
  write.table(df_sub, file=paste0(out_dir,rowid,".bed"), 
              sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
  
  return(list(df_sub,content))
}

rmskAnno<-function(rowid){
  #define family and class lists
  if(ref_species=="hg38"){
    class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','DNA','Unknown')
  }else{
    class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','DNA','Unknown','Other')
  }
  
  family_list=unique(rmsk_GRCm38$repFamily)
  class_list=unique(rmsk_GRCm38$repClass)
  
  # if annotation id (rowid) is in the defined annotation list, subset, and outuput bed
  annotation_list=unique(c(class_list,'7SK RNA','yRNA'))
  
  ## PH -elseif statements weren't being recognized for some reason
  if(rowid %in% annotation_list){
    #define cols to select and new col names
    rmskcol=c('genoName','genoStart','genoEnd','repName','swScore','strand')
    newcolnames=c('chr','start','end','transcript_name','score','strand')
    #if in the family or class list
    if (rowid %in% class_list){
      df_sub = subset(rmsk_GRCm38, repClass == rowid) %>%
        select(rmskcol) %>%
        setnames(newcolnames)
    }
    if(rowid== "7SK RNA"){
      df_sub = subset(rmsk_GRCm38, repFamily == 'RNA') %>%
        select(rmskcol) %>%
        setnames(newcolnames)
      rowid="7SKRNA"
    } 
    if (rowid == "LINE SINE"){
      df_sub = subset(rmsk_GRCm38, repClass == c("LINE","SINE")) %>%
        select(rmskcol) %>%
        setnames(newcolnames)
    } 
    if (rowid == "yRNA"){
      df_sub = subset(rmsk_GRCm38, repFamily == 'scRNA' & grepl("HY",rmsk_GRCm38$repName)==TRUE) %>%
        select(rmskcol) %>%
        setnames(newcolnames)
    } 
    if (rowid == "scRNA"){
      df_sub = subset(rmsk_GRCm38, repFamily == 'scRNA' & grepl("HY",rmsk_GRCm38$repName)==FALSE) %>%
        select(rmskcol) %>% setnames(newcolnames)
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
  file_name = ref_table[rowid,paste0(ref_species,"_SY")]
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


##########################################################################################
############### Create annotation_output table 
##########################################################################################

annotation_output=data.frame()

#input data from annotation df
ref_table = read.csv(reftable_path, header = TRUE, sep="\t", row.names = 1)
#rowid="lincRNA"
for (rowid in rownames(ref_table)){
  #determine source to use
  ref_selection =  ref_table[rowid,paste0(ref_species,"_selection")]
  ref_selection=unlist(strsplit(ref_selection,split ="," ))  
  
  
  for (x in 1:length(ref_selection) ) {
    ref_selectionx=ref_selection[x]
    print(paste0(rowid,"-",ref_selectionx))
    
    df_sub=data.frame()
    content=""    
    ### Gencode annotations
    if(ref_selectionx=="Gencode"){
      return_list = gencodeAnno(rowid)
      df_sub = return_list[[1]]
      content = return_list[[2]]
      
      ### REFSEQ annotations  
    } 
    if (ref_selectionx == "RefSeq"){
      print("add refseq")
      
      ### Repeatmasker annotations
    } 
    if (ref_selectionx == "Repeatmasker"){
      return_list = rmskAnno(rowid)
      df_sub = return_list[[1]]
      content = return_list[[2]]
      
      ###Soyeong references
    } 
    if (ref_selectionx%in%c("Gencode","Repeatmasker")==F){
      return_list = SYAnno(rowid,ref_species)
      df_sub = return_list[[1]]
      content = return_list[[2]]
      
      ###Invalid entry
    } 
    
    #print into to output
    if (length(ref_selection)>1) {
      annotation_output[paste0(rowid,"_",ref_selectionx),"contents"]  = content
      annotation_output[paste0(rowid,"_",ref_selectionx),"count"] = nrow(df_sub)  
      annotation_output[paste0(rowid,"_",ref_selectionx),'source'] = ref_selection[x]
      annotation_output[paste0(rowid,"_",ref_selectionx),'description'] = ref_table[rowid,"description"]
      
      if (rowid=='lincRNA' & ref_selectionx =='Gencode') {annotation_output[paste0(rowid,"_",ref_selectionx),"count"] = length(unique(df_sub$transcript_id))}
      
    } 
    if (length(ref_selection)==1){
      annotation_output[rowid,"contents"]  = content
      annotation_output[rowid,"count"] = nrow(df_sub)  
      annotation_output[rowid,'source'] = ref_selection
      annotation_output[rowid,'description'] = ref_table[rowid,"description"]
      
      if (rowid=='lincRNA'&ref_selectionx=='Gencode') {annotation_output[rowid,"count"] = length(unique(df_sub$transcript_id))}
    }
  }
}

#hard code contents
if (ref_species=='mm10') {
  annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12 and various predicted genes'
  annotation_output['rRNA_Gencode','contents']='5S, 5.8s, predicted gene'
} else if (ref_species=='hg38') {
  annotation_output['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12'
  annotation_output['rRNA_Gencode','contents']='5S, 5.8s, predicted gene'
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