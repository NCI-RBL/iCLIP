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
testing="Y"
if(testing=="Y"){
  ref_dir = "/Volumes/RBL_NCI/iCLIP/ref/annotations/"
  ref_dir = "/Users/homanpj/Documents/Resources/ref/"
  ref_species = "hg38" ### need better name, match to snakemake
  reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/iCLIP/config/annotation_config.txt"
  
  ### fix need to figure out how to keep all this info - maybe a config? dict?
  alias_path = paste0(ref_dir,ref_species,"/",ref_species,".chromAlias.txt")
  if(ref_species == "mm10"){
    out_dir = "/Users/homanpj/Desktop/CLIPtest/annotation/mm10/"
    gencode_path = paste0(ref_dir, "mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "/mm10/NCBI_RefSeq/GCF_000001635.26_GRCm38.p6_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"/mm10/Gencode_VM23/fromUCSC/KnownCanonical/KnownCanonical_GencodeM23_GRCm38.txt")
    intron_path = paste0(ref_dir, "/mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
    rmsk_path = paste0(ref_dir,"/mm10/repeatmasker/rmsk_GRCm38.txt")
    soyeong_flag = "Y" #Y or N
    soyeong_path= paste0(ref_dir,'/mm10/AdditionalAnno/')
  } else if (ref_species == "hg38"){
    out_dir = "/Users/homanpj/Desktop/CLIPtest/annotation/hg38/"
    gencode_path = paste0(ref_dir,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
    refseq_path = paste0(ref_dir, "/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt")
    canonical_path = paste0(ref_dir,"/hg38/Gencode_V32/fromUCSC/KnownCanonical/KnownCanonical_GencodeM32_GRCh38.txt")
    intron_path = paste0(ref_dir,"/hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
    rmsk_path = paste0(ref_dir,"/hg38/repeatmasker/rmsk_GRCh38.txt")
    soyeong_flag = "N" #always now (for now)
  } 
}

##########################################################################################
############### Annotation info
##########################################################################################
#### Lets not include this as part of the script. 
#### I think if they want to bring in refseq annotations that can be an idditional imput


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

###phil why are we merging these all into ncRNA? Can we leave as separate types?
## combine all ncRNA
#### we created a column $gene_type_ALL to keep track of ncRNA sub catagories. 
#### when we plot CLIP peaks we plot both ncRNA as a single group and ncRNA as subgroups
ReplaceType <-function(df_in,col_in,list_in,replace_id){
  for (id in list_in){
    df_in[,col_in]=gsub(id,replace_id,df_in[,col_in])
  } 
  return(df_in)
}

type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA","tRNA","ribozyme")
ref_gencode = ReplaceType(ref_gencode,"gene_type", type_list,"ncRNA")

#create subsets
### PH - Subsets used to annotate clip peaks seperatly 
###phil FTR is never used PH-remove FTR
ref_gencode_t = subset(ref_gencode, feature == 'transcript')
ref_gencode_e = subset( ref_gencode, feature == "exon")

# ref_gencode_te = subset( ref_gencode, feature == c("transcript","exon"))


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

###phil why are we adding 1 to start and exon? 
### PH - intron table is 0 based while exon starts at 1
### same for exon number the exon # starts at 0
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

###phil why is this flag set? none of the following code is used with flag
###phil this file doesn't exist
###PH- I created this early on to Idneify Intron locations from GTF file (only includes exon coordinates)
### we can extract this section of code and include it as a tool when updating reference GTF
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

## PH - to save time i have added Gene type to intron table then added to exons
# intron_exon=merge(introns,
#                   ref_gencode_t[,c('transcript_id','gene_type')],by='transcript_id',all.x=T)



##########################################################################################
############### REPEATMASKER
##########################################################################################
# get  repeat regions and extra small RNA locations

rmsk_GRCm38=fread(rmsk_path, header=T, sep="\t",stringsAsFactors = F,
                  data.table=F)


##########################################################################################
############### Annotation functions
##########################################################################################
gencodeLNCRNA<-function(){
  ## PH- Simplified function
  ## What i am trying to do is identfy a subset of lncRNA that have intronic sequences.
  ## we are then combining the Intronic and Exonic locations and changing the gene_type_ALL from lnc to linc
  ## finally we will combibne the linc + lnc (no introns) together into one lnc table
 
   #####
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

  ##PH- because I have somewhat redefined lncRNA with linc and included exon and intron sequences I am 
  ## treating lncRNA as separate annotation table and removing these transcripts from ref_gencode_t
  
  ref_gencode_t=ref_gencode_t[!ref_gencode_t$transcript_id%in%unique(lncRNA_ref_gencode$transcript_id),]
  
  write.csv(lncRNA_ref_gencode,paste0(out_dir,"ref_lncRNA.csv"))
  
  return(lncRNA_ref_gencode)
}



rowid='lincRNA'
# if(ref_table[ref_table%in%rowid,paste0(ref_species,'_mm10_options')])
  
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
  # family_list = c("rRNA","snRNA", "srpRNA","tRNA")
  #remove other from hg38 - not found in rmsk_GRCm38
  if(ref_species=="hg38"){
    class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','DNA','Unknown')
  }else{
    class_list =c("Satellite","Low_complexity","LTR",'Simple_repeat','DNA','Unknown','Other')
  }
  family_list=unique(rmsk_GRCm38$repFamily)
  class_list=unique(rmsk_GRCm38$repClass)
  #if annotation id (rowid) is in the defined annotaiton list, subset, and outuput bed
  # annotation_list = c(family_list,class_list,"7SK RNA","LINE SINE","rRNA_DNA","scRNA", "yRNA")
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
      #Special handling  
    }
    # if(rowid== "tRNA"){rep family contains tRNA labeling that includes SINE/LINE in repCLASS
    #   df_sub = subset(rmsk_GRCm38, repFamily == 'tRNA') %>%
    #     select(rmskcol) %>%
    #     setnames(newcolnames)
    # } 
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

###phil - review the reftable to fill in missing data
#input data from annotation df
ref_table = read.csv(reftable_path, header = TRUE, sep="\t", row.names = 1)
rowid="lincRNA"
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
###phil need to fix
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
