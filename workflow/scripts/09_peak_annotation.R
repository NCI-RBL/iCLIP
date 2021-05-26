#library
library(tidyr)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(data.table)

#TO DO
#fix variables
#remove all paths and link directly in snakemake
# fix arg list

args <- commandArgs(trailingOnly = TRUE)
peak_type = args[1]
peak_unique = args[2]
peak_all = args[3]
join_junction = args[4]
read_depth = args[5]
DEmethod = args [6]
sample_id = args[7]
nt_merge = args[8]
out_dir = args[9]
anno_dir = args[10]
ref_species = args[11]
reftable_path = args[12]
ref_dir = args[13]


if(length(args)==0){
  wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/iCLIP_Git/"
  setwd(wd)
  wd="."
    
  peak_type= "ALL"
  peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/NextSeq_fCLIP_0218/13_counts/allreadpeaks/WT1_fCLIP_50nt_uniqueCounts.txt"
  peak_all = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/NextSeq_fCLIP_0218/13_counts/allreadpeaks/WT1_fCLIP_50nt_allFracMMCounts.txt"
  reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/iClip_Git/config/annotation_config.txt"
  
  #output
  out_dir = paste0(wd,"/workflow/scripts/AnnoTestPIPE/")
  out_file = paste0(wd,"/workflow/scripts/AnnoTestPIPE/15_peaks/WT1_fCLIP_50nt_peakannotation_complete.txt")
  
  #project annotation files
  rmsk_path = "/Users/homanpj/Documents/Resources/ref/hg38/repeatmasker/rmsk_GRCh38.txt"
  anno_dir = paste0(wd,"/workflow/scripts/AnnoTestPIPE/15_project_annotation/")
    ref_dir = "/Users/homanpj/Documents/Resources/ref/"
  
  #feature information
  join_junction = "TRUE"
  read_depth = 5
  DEmethod = "MANORM"
  ref_species="hg38"
  sample_id = "WT1"
  nt_merge = "50nt"
}

removeVersion <- function(ids){
  return(unlist(lapply(stringr::str_split(ids, "[.]"), "[[",1)))
}

varname <- function(x) {
  deparse(substitute(x))}

#annotation paths
gencode_transc_path = paste0(anno_dir,"ref_gencode.txt")
lncra_path = paste0(anno_dir,"lncRNA_gencode.txt")
alias_path =  paste0("/Users/homanpj/Documents/Resources/ref/mm10/mm10.chromAlias.txt") 
YRNA_path = paste0(anno_dir, "yRNA_repeatmasker.bed")
srpRNA_path = paste0(anno_dir, "srpRNA_repeatmasker.bed")
SKRNA_path = paste0(anno_dir, "7SKRNA_repeatmasker.bed")
scRNA_path = paste0(anno_dir, "scRNA_repeatmasker.bed")
tRNA_path = paste0(anno_dir, "tRNA_repeatmasker.bed")
sncRNA_path = paste0(anno_dir, "sncRNA_gencode.bed")
rRNA_BK00964_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
rRNA_rmsk_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
tRNA_rmsk_path = paste0(anno_dir, "tRNA_repeatmasker.bed")

if(ref_species == "mm10"){
  gencode_path = paste0(ref_dir, "mm10/Gencode_VM23/fromGencode/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.txt")
  intron_path = paste0(ref_dir, "mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
  rmsk_path = paste0(ref_dir,"mm10/repeatmasker/rmsk_GRCm38.txt")

} else if (ref_species == "hg38"){
  gencode_path = paste0(ref_dir,"hg38/Gencode_V32/fromGencode/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.txt")
  intron_path = paste0(ref_dir,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
  rmsk_path = paste0(ref_dir,"hg38/repeatmasker/rmsk_GRCh38.txt")
} 

#set id for files
file_id = paste0(sample_id,"_",nt_merge,"_")

##########################################################################################
############### unique read Count
##########################################################################################
FtrCount_uniq=read.delim(peak_unique, header=T, sep="\t",
                         stringsAsFactors = F,comment.char = '#')
colnames(FtrCount_uniq)[7]='Counts'; #FtrCount_uniq[is.na(FtrCount_uniq$Start),]
FtrCount_uniq$Start=as.numeric(FtrCount_uniq$Start)#+1
FtrCount_uniq$End=as.numeric(FtrCount_uniq$End)
FtrCount_uniq=FtrCount_uniq[!is.na(FtrCount_uniq$Start),]
FtrCount_uniq$ID=paste0(FtrCount_uniq$Chr,":",FtrCount_uniq$Start,"-",FtrCount_uniq$End)
FtrCount_uniq$ID2=paste0(FtrCount_uniq$Chr,":",FtrCount_uniq$Start,"-",FtrCount_uniq$End,"_",FtrCount_uniq$Strand)

##########################################################################################
############### all reads
##########################################################################################
FtrCount_frac=read.delim(peak_all, header=T, sep="\t",
                         stringsAsFactors = F,comment.char = '#')
colnames(FtrCount_frac)[7]='Counts'
FtrCount_frac$Start=as.numeric(FtrCount_frac$Start)#+1
FtrCount_frac$End=as.numeric(FtrCount_frac$End)
FtrCount_frac=FtrCount_frac[!is.na(FtrCount_frac$Start),]
FtrCount_frac$ID=paste0(FtrCount_frac$Chr,":",FtrCount_frac$Start,"-",FtrCount_frac$End)

##########################################################################################
############### merge
##########################################################################################
FtrCount_merged=merge(FtrCount_uniq[,(colnames(FtrCount_uniq) %in% "Geneid")==F],
                      FtrCount_frac[,c('ID','Counts')],
                      by='ID',suffixes=c("_Unique","_fracMM"),all.x=T)
colnames(FtrCount_merged)[which(colnames(FtrCount_merged)%in%c('Chr','Start','End','Strand'))]=c('chr','start','end','strand')
FtrCount_merged=FtrCount_merged[duplicated(FtrCount_merged)==F,]

##########################################################################################
##############################  splice junctions
##########################################################################################
if (join_junction==TRUE) {
  
  print("Running Join Junction")
  
  #read in jcounts file
  if (peak_type=="UNIQUE") {
    FtrCount_fracJCount = read.delim(paste0(peak_unique,".jcounts"),
                                     header=TRUE,sep="\t",stringsAsFactors = FALSE,comment.char = '#') 
    paste(nrow(FtrCount_fracJCount),"Unique")
  }  else {
    FtrCount_fracJCount = read.delim(paste0(peak_all,".jcounts"), 
                                     header=TRUE,sep="\t",stringsAsFactors = FALSE,comment.char = '#')
    paste(nrow(FtrCount_fracJCount),"all")
  } 
  
  #rename last col
  colnames(FtrCount_fracJCount)[ncol(FtrCount_fracJCount)]='counts'
  
  #clear na's
  FtrCount_fracJCount = FtrCount_fracJCount[is.na(FtrCount_fracJCount$PrimaryGene)==FALSE,]
  
  ###phil need to determine what to do if there are no splice junctions identified
  #if there are junctions identifitied merge
  if(nrow(FtrCount_fracJCount)>0){
    print("Junctions were identitified")
    JoinJunc=T
    
    FtrCount_fracJCount_merged = merge(FtrCount_fracJCount, 
                                       FtrCount_merged[,c("ID","strand")],by.x="PrimaryGene",by.y="ID",all.x=T)
    #split primary gene col
    FtrCount_fracJCount_merged = separate(FtrCount_fracJCount_merged,
                                          PrimaryGene,into=c('chr','start','end'),sep = ":|-",remove = F)
    
    ### remove junctions where splicing is within peak
    site_1 = FtrCount_fracJCount_merged$Site1_location >=FtrCount_fracJCount_merged$start&FtrCount_fracJCount_merged$Site1_location<=FtrCount_fracJCount_merged$end
    site_2 = FtrCount_fracJCount_merged$Site2_location>=FtrCount_fracJCount_merged$start&FtrCount_fracJCount_merged$Site2_location<=FtrCount_fracJCount_merged$end
    FtrCount_fracJCount_merged=FtrCount_fracJCount_merged[site_1!=site_2,]
    FtrCount_fracJCount_merged$JunctionID = paste0(FtrCount_fracJCount_merged$Site1_chr,':',
                                                   FtrCount_fracJCount_merged$Site1_location,'-', 
                                                   FtrCount_fracJCount_merged$Site2_location)
    
    #create Unique Row/junction ID
    FtrCount_fracJCount_merged$rID=seq(1,nrow(FtrCount_fracJCount_merged))
    
    ###phil there is a row counts error here; changed the rID from FtrCount_fracJCount$rID
    #to FtrCount_fracJCount_merged$rID
    jcount1.GR = GRanges(seqnames = as.character(FtrCount_fracJCount_merged$Site1_chr), 
                         ranges=IRanges(start = as.numeric(FtrCount_fracJCount_merged$Site1_location), 
                                        end = as.numeric(FtrCount_fracJCount_merged$Site1_location)),
                         strand = FtrCount_fracJCount_merged$strand,
                         PrimaryGene = FtrCount_fracJCount_merged$PrimaryGene,
                         count = FtrCount_fracJCount_merged$counts,
                         JunctionID = FtrCount_fracJCount_merged$JunctionID,
                         rID = FtrCount_fracJCount_merged$rID)
    jcount2.GR = GRanges(seqnames = as.character(FtrCount_fracJCount_merged$Site2_chr),
                         ranges=IRanges(start = as.numeric(FtrCount_fracJCount_merged$Site2_location),
                                        end = as.numeric(FtrCount_fracJCount_merged$Site2_location)),
                         strand = FtrCount_fracJCount_merged$strand,
                         PrimaryGene=FtrCount_fracJCount_merged$PrimaryGene,
                         count=FtrCount_fracJCount_merged$counts,
                         JunctionID=FtrCount_fracJCount_merged$JunctionID,
                         rID=FtrCount_fracJCount_merged$rID)
    
    FtrCount.GR=GRanges(seqnames = as.character(FtrCount_merged$chr), 
                        ranges=IRanges(start = as.numeric(FtrCount_merged$start), 
                                       end = as.numeric(FtrCount_merged$end)),
                        strand = FtrCount_merged$strand,
                        ID=FtrCount_merged$ID)
    
    #Add Junction ID
    junctionID <- function(Ftr_in, j_in, j_name, p_name){
      xo=as.data.frame(GenomicRanges::findOverlaps(j_in, Ftr_in, 
                                                   type = "any",ignore.strand=F))
      qh=as.data.frame(j_in[xo$queryHits],row.names = NULL)
      colnames(qh)=paste0(colnames(qh),j_name)
      
      sh=as.data.frame(Ftr_in[xo$subjectHits],row.names = NULL)
      colnames(sh)=paste0(colnames(sh),p_name)
      
      return(cbind(qh,sh))
    }
    
    Junc_PL1 = junctionID(FtrCount.GR, jcount1.GR,"_Junction1","_Peaks1")
    Junc_PL2 = junctionID(FtrCount.GR, jcount2.GR,"_Junction2","_Peaks2")
    Junc_PL_merged = merge(Junc_PL1,Junc_PL2,by.x='rID_Junction1',by.y='rID_Junction2')
    Junc_PL_merged = Junc_PL_merged[(is.na(Junc_PL_merged$ID_Peaks1)|is.na(Junc_PL_merged$ID_Peaks2))==F,]
    
    # identify PEAK ID for Site1 and Site2 location of splice junction
    Junc_PLOut = Junc_PL_merged[,c('PrimaryGene_Junction1','JunctionID_Junction1',
                                   'strand_Peaks1','ID_Peaks1','ID_Peaks2')]
    
    # Identify peaks that are linked between splice junctions
    # Remove Peaks where PrimaryGene - Peak1, Peak2 are duplicated
    # Junction ID can be off by a few NT because different splice start/end in peak
    Junc_PLOut=Junc_PLOut[duplicated(Junc_PLOut[,c(1,4,5)])==F,]
    
    # Combine by primaryGene ID - Junction ID acts as ID to distinguish all possible junctions
    # get all unique PrimaryGene_Junction
    prim_junct=unique(Junc_PLOut$PrimaryGene_Junction1)
    PGene_TBL = setNames(data.frame(matrix(nrow=length(prim_junct), ncol=3)),
                         c('PrimaryGene_Junction1','JunctionID_Junction1','ID_comb'))
    
    for (junctions in 1:length(prim_junct)) {
      PrimGene=prim_junct[junctions]
      g=Junc_PL_merged[Junc_PL_merged$PrimaryGene_Junction1%in%PrimGene,]
      PGene_TBL[junctions,'PrimaryGene_Junction1']=unique(g$PrimaryGene_Junction1)
      PGene_TBL[junctions,'JunctionID_Junction1']=paste(unique(sort(g$JunctionID_Junction1)),collapse = ",")
      PGene_TBL[junctions,'ID_comb']=paste(unique(c(g$ID_Peaks1,g$ID_Peaks2)),collapse = ",")
    }
    #PGene_TBL=as.data.frame(PGene_TBL)
    #remove duplicated junctions
    PGene_TBL_unique=PGene_TBL[duplicated(PGene_TBL)==FALSE,]
    
    ##### Combine all connected peaks
    # get all unique PrimaryGene_Junction
    prim_junct = unique(PGene_TBL_unique$PrimaryGene_Junction1)
    PGene_TBL2=(matrix(nrow=nrow(PGene_TBL_unique),ncol=ncol(PGene_TBL_unique)+4))
    colnames(PGene_TBL2)=c(colnames(PGene_TBL_unique),'chr','start','end','linkedID')
    
    for (junction in 1:length(prim_junct)) {
      a = PGene_TBL[junction,]
      id = unlist(str_split (a$ID_comb,pattern = ","))
      idcomb = PGene_TBL[grep(paste(id,collapse = "|"),PGene_TBL$ID_comb),]
      
      ## look again to see if more locations with found connected peaks
      id2 = sort(unique(unlist(str_split (idcomb$ID_comb,pattern = ","))))
      
      ## If new peaks come up check until no new peaks  
      while (FALSE%in%(id2%in%id)) {
        id=id2
        
        # Look Again
        idcomb=PGene_TBL[grep(paste(id2,collapse = "|"),PGene_TBL$ID_comb),]
        
        # xxx3 Get Peak IDS to end while loop
        id2=sort(unique(unlist(str_split (idcomb$ID_comb,pattern = ","))))
      }
      
      id2_coord = separate(as.data.frame(id2),1,sep = ":|-",
                           into = c('chr','start','end'))
      id2_coord$start=as.numeric(id2_coord$star)
      id2_coord$end=as.numeric(id2_coord$end)
      PGene_TBL2[junction,'PrimaryGene_Junction1']=paste(sort(unique(idcomb$PrimaryGene_Junction1)),collapse = ',')
      PGene_TBL2[junction,'JunctionID_Junction1']=paste(sort(unique(idcomb$JunctionID_Junction1)),collapse = ',')
      PGene_TBL2[junction,'ID_comb']=paste(sort(unique(idcomb$ID_comb)),collapse = ',')
      PGene_TBL2[junction,'chr']=paste(unique(id2_coord$chr),collapse=",") ###phil added paste here
      PGene_TBL2[junction,'start']=min((id2_coord[,c('start','end')]))
      PGene_TBL2[junction,'end']=max((id2_coord[,c('start','end')]))
    }
    
    PGene_TBL2=as.data.frame(PGene_TBL2)
    PGene_TBL2=PGene_TBL2[duplicated(PGene_TBL2)==F,]
    PGene_TBL2$linkedID=paste0(PGene_TBL2$chr,":",PGene_TBL2$start,"-",PGene_TBL2$end)
    
    ##### check no reeated peaks 
    ##### get all unique PrimaryGene_Junction
    d3=unique(PGene_TBL2$PrimaryGene_Junction1)
    PGene_TBL3 = (matrix(nrow=nrow(PGene_TBL2),
                         ncol=ncol(PGene_TBL2)+1));colnames(PGene_TBL3)=c(colnames(PGene_TBL2),"nrows")
    
    for (x in 1:length(d3)) {
      a2=PGene_TBL2[x,]
      id2=unlist(str_split (a2$ID_comb,pattern = ","))
      idcomb2=PGene_TBL2[grep(paste(id2,collapse = "|"),PGene_TBL2$ID_comb),]
      
      PGene_TBL3[x,1:3]=as.matrix(a2[,1:3])
      PGene_TBL3[x,'nrows']=nrow(idcomb2)
    }
    PGene_TBL3=as.data.frame(PGene_TBL3)
    
    # ##### Get all peaks with junctions
    Junc_peaks=unique(c(Junc_PL_merged$ID_Peaks1,Junc_PL_merged$ID_Peaks2))
    FtrCount_Junc_peaks=FtrCount_merged[FtrCount_merged$ID%in%Junc_peaks,]
    FtrCount_merged_junction=FtrCount_merged[FtrCount_merged$ID%in%Junc_peaks==F,]
    
    ### Trim peaks without Splicing
    FtrCount_nosplice=FtrCount_merged_junction[FtrCount_merged_junction$Counts_fracMM>=read_depth,]
    FtrCount_splice_junc=rbind(FtrCount_nosplice,FtrCount_Junc_peaks)
    FtrCount_splice_junc=FtrCount_splice_junc[duplicated(FtrCount_splice_junc)==F,]
    FtrCount_splice_junc=FtrCount_splice_junc[!FtrCount_splice_junc$chr%in%c('chrM'),]
    
    if (DEmethod=='MANORM') {
      print("Running MANORM")
      FtrCount_splice_junc$Counts_fracMM >
        write.table(FtrCount_splice_junc[, c('chr','start','end','ID','ID','strand')],
                    file=paste0(out_dir, sample_id,"_", nt_merge,"_peakDepth.bed"), 
                    sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
      write.table(FtrCount_splice_junc[FtrCount_splice_junc$strand%in%"+", c('chr','start','end','ID','ID','strand')],
                  file=paste0(out_dir,sample_id,"_",nt_merge,"_peakDepth_P.bed"), 
                  sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
      write.table(FtrCount_splice_junc[FtrCount_splice_junc$strand%in%"-", c('chr','start','end','ID','ID','strand')],
                  file=paste0(out_dir,sample_id,"_",nt_merge,"_peakDepth_N.bed"), 
                  sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
    } else{
      print("MANORM skipped")
    } 
    
    #format final table
    colnames(FtrCount_splice_junc)[which(colnames(FtrCount_splice_junc) %in%
                                           c("V1","V2","V3","V6"))]=c('chr','start','end','strand')
    
  } else{
    print("No junctions were identified")
    ###phil if there are no junctoins, then just use the merge?
    FtrCount_splice_junc=FtrCount_merged[c("chr","start","end","strand")]
    JoinJunc=F
  }
}  else {
  print("Junctions joining not selected")
  FtrCount_splice_junc=FtrCount_merged[c("chr","start","end","strand")]
  JoinJunc=F
}

#write final junction output
write.table(FtrCount_splice_junc,paste0(out_dir,file_id,"peakjunction.txt"),sep = "\t")
## end join_junc
##########################################################################################
##############################  Peak info - starts old 07_peak_generation
##########################################################################################
#read in peaks and alias file generated from 06_peak_junction.R
peaks = FtrCount_splice_junc

peaks_Oppo=peaks 
peaks_Oppo$strand=gsub("\\+","pos",peaks_Oppo$strand) 
peaks_Oppo$strand=gsub("\\-","+",peaks_Oppo$strand)
peaks_Oppo$strand=gsub("pos","-",peaks_Oppo$strand)

##########################################################################################
##############################  reference, annotation table
##########################################################################################
#read in references
# ref_lncra = read.table(lncra_path,header=TRUE,sep="\t")
ref_gencode_trans = read.table(gencode_transc_path,header=TRUE,sep="\t")

### Add Column that combines all additional annotations
# YRNA_rmsk = read.table(YRNA_path,header=TRUE,sep="\t")
# srpRNA_rmsk = read.table(srpRNA_path,header=TRUE,sep="\t")
# SKRNA_rmsk = read.table(SKRNA_path,header=TRUE,sep="\t")
# scRNA_rmsk = read.table(scRNA_path,header=TRUE,sep="\t")

#phil we should think about how to organize this - we should pull in from anno_config? 
#can also call the annotation script into this if we dont need to create BEDS
###phil this assumes that SY flag - is there no option here? if so need to update previous code

ref_table = read.table(reftable_path, header = TRUE, sep="\t", row.names = 1)


rpmk_Anno_RNA_comb=data.frame()
gencode_Anno_RNA_comb=data.frame()
  gencodeCol=c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')
  gencode_Anno_RNA_comb=ref_gencode_trans[,gencodeCol]  
custom_Anno_RNA_comb=data.frame()
refseq_Anno_RNA_comb=data.frame()

for (rowid in rownames(ref_table)){
  ref_type = ref_table[rowid,'AnnoType']
  ref_selection =  ref_table[rowid,paste0(ref_species,"_selection")]
  ref_selection=unlist(strsplit(ref_selection,split ="," ))  
  rowid=gsub(" ","",rowid)
  
  if (ref_type=='transcript_anno') {
    for (selec in ref_selection) {
    if (selec =='Gencode') {
      print(paste0('Combine ',paste0(rowid,'_',selec,'.bed')))
      Gtable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=T,sep="\t")
      gencode_Anno_RNA_comb=rbind(gencode_Anno_RNA_comb,Gtable[,gencodeCol])
    } else if (selec =='Repeatmasker') {
      print(paste0('Combine ',paste0(rowid,'_',selec,'.bed')))
      Rtable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=TRUE,sep="\t")
      rpmk_Anno_RNA_comb=rbind(rpmk_Anno_RNA_comb,Rtable)
    } else if (selec =='Custom') {
      print(paste0('Combine ',paste0(rowid,'_',selec,'.bed')))
      Ctable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=TRUE,sep="\t")
      custom_Anno_RNA_comb=rbind(custom_Anno_RNA_comb,Ctable)
    } else if (selec =='refseq') {
      print(paste0('Combine ',paste0(rowid,'_',selec,'.bed')))
      Reftable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=TRUE,sep="\t")
      refseq_Anno_RNA_comb=rbind(refseq_Anno_RNA_comb,Reftable)
    }

    }
  }else { next }
}  

gencode_Anno_RNA_comb=gencode_Anno_RNA_comb[duplicated(gencode_Anno_RNA_comb)==F,]
rpmk_Anno_RNA_comb=rpmk_Anno_RNA_comb[duplicated(rpmk_Anno_RNA_comb)==F,]
custom_Anno_RNA_comb=custom_Anno_RNA_comb[duplicated(custom_Anno_RNA_comb)==F,]
refseq_Anno_RNA_comb=refseq_Anno_RNA_comb[duplicated(refseq_Anno_RNA_comb)==F,]

annocol=c("chr","start","end","strand","type","transcript_name")

# if (dim(gencode_Anno_RNA_comb)[1] == 0) {
#   gencode_Anno_RNA_comb=as.data.frame(matrix(ncol=length(annocol),nrow=0))
#   colnames(gencode_Anno_RNA_comb)=annocol
# } else {gencode_Anno_RNA_comb=gencode_Anno_RNA_comb[,c("chr","start","end","strand","gene_type","ensembl_gene_id")]
# colnames(gencode_Anno_RNA_comb)=annocol}


if (dim(custom_Anno_RNA_comb)[1] == 0) {
  custom_Anno_RNA_comb=as.data.frame(matrix(ncol=length(annocol),nrow=0))
  colnames(custom_Anno_RNA_comb)=annocol
} else {custom_Anno_RNA_comb=custom_Anno_RNA_comb[,c("chr","start","end","strand","gene_type","ensembl_gene_id",Gene)]
        colnames(custom_Anno_RNA_comb)=annocol}

if (dim(rpmk_Anno_RNA_comb)[1] == 0) {
  rpmk_Anno_RNA_comb=as.data.frame(matrix(ncol=length(annocol),nrow=0))
  colnames(rpmk_Anno_RNA_comb)=annocol
} else {rpmk_Anno_RNA_comb=rpmk_Anno_RNA_comb[,c('chr','start','end','strand','type','name')]
        colnames(rpmk_Anno_RNA_comb)=annocol}

if (dim(refseq_Anno_RNA_comb)[1] == 0) {
  refseq_Anno_RNA_comb=as.data.frame(matrix(ncol=length(annocol),nrow=0))
  colnames(refseq_Anno_RNA_comb)=annocol
} else {refseq_Anno_RNA_comb=refseq_Anno_RNA_comb[,c("chr","start","end","strand","gene_type","ensembl_gene_id")]
        colnames(refseq_Anno_RNA_comb)=annocol}


Anno_RNA_comb=rbind(
      # gencode_Anno_RNA_comb,
      rpmk_Anno_RNA_comb,
      custom_Anno_RNA_comb,
      refseq_Anno_RNA_comb
)

Anno_RNA_comb=Anno_RNA_comb[duplicated(Anno_RNA_comb)==F,]

### above should create anno tables based on config table
# if(ref_species == "mm10"){
#   #pull in references
#   tRNA_sy = read.table(tRNA_path,header=TRUE,sep="\t")
#   sncRNA_sy = read.table(sncRNA_path,header=TRUE,sep="\t")
#   rRNA_BK00964 = read.table(rRNA_BK00964_path,header=TRUE,sep="\t")
#   
#   Anno_RNA_comb=rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],
#                       tRNA_sy[,annocol],sncRNA_sy[,annocol],rRNA_BK00964[,annocol],
#                       SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],rRNA_rmsk[,annocol])
# } else if (ref_species == "hg38"){
#   #pull in references
#   rRNA_rmsk = read.table(rRNA_rmsk_path,header=TRUE,sep="\t")
#   tRNA_rmsk = read.table(tRNA_rmsk_path,header=TRUE,sep="\t")
#   
#   Anno_RNA_comb = rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],
#                         SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],
#                         rRNA_rmsk[,annocol],tRNA_rmsk[,annocol])
# }

##########################################################################################
##############################  call peaks
##########################################################################################

#merge references depending on species, run annotation
bam_anno2=function(peaksTable,Annotable,ColumnName,pass_n){
  #testing
  if(length(args)==0){
    peaksTable=peaks_Oppo
    Annotable=Anno_RNA_comb
    ColumnName=c('type','transcript_name')
    pass_n='pass2'
    # 
    # peaksTable =peaks_OppoAnno[,c('chr','start','end','ID','ID2','strand')]
    # Annotable=anno_tmp
    # ColumnName=c('feature','exon_number')
    # pass_n="pass3"
  }

  toString(anno_dir)
  ###phil why are we adding multiple chr cols? 
  ## the multiple chr is a space filler was to make a bed6 file for bedtools 
  anno_output = Annotable[,c('chr','start','end','chr','chr','strand',ColumnName)] %>%
  # anno_output = Annotable[,c('chr','start','end','strand',ColumnName)] %>%
    dplyr::rename(
      chr_anno = chr,
      start_anno = start,
      end_anno = end,
      strand_anno = strand
    )
  ###phil why are we clearing these? we could only input the ones we want
  ## this was a simple formatting issue for creating a bed 6 file
  #clear chr cols 4,5
  anno_output[,4:5]="."
  
  #set file prefix
  file_prefix = paste0(file_id,pass_n,"_")
  
  ###phil there are non chrm found in the peaksTable - need to remove for merge
  #output bed files
  write.table(anno_output,
              file = paste0(out_dir,file_prefix,"annotable.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  
  peaksTable_output=peaksTable[is.na(peaksTable$start)==F,]
  peaksTable_output=subset(peaksTable_output,chr %in% unique(anno_output$chr_anno))
  peaksTable_output=peaksTable_output[,c('chr','start','end','ID','ID2','strand')]
  write.table(peaksTable_output,
              file = paste0(out_dir,file_prefix,"peakstable.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  #remove lines with illegal characters, remove NAs
  ## this was not working for me so did similar command above
  # cmd = dQuote("NA")
  # system(paste0("cat ",anno_dir,file_prefix,"peakstable.bed | awk '$2 !~ /e/' | awk '$2 != ", cmd, "' > ",
                # anno_dir,file_prefix,"peakstable_clean.bed"))
  
  ###phil i think we can do this with GRRanges? 
  ### i liked using bedtools because i have more control over the overlap parameters. I have an older version that uses GRanges
  # merge bedfiles into output text
  system(paste0('bedtools intersect -a ', 
                # anno_dir, file_prefix,'peakstable_clean.bed -b ', 
                out_dir, file_prefix,'peakstable.bed -b ', 
                out_dir, file_prefix,'annotable.bed -wao -s > ',
                out_dir, file_prefix,'peaks_OL.txt'))
  
  
  #read in merged file
  ab_OL = read.table(paste0(out_dir,file_prefix,"peaks_OL.txt"), 
                   header=F, sep="\t",stringsAsFactors = F)
  colnames(ab_OL) = c(paste0(colnames(peaksTable_output)),
                      paste0(colnames(anno_output)),'ntOL')
  head(ab_OL)
  
  #add width
  ab_OL$width = ab_OL$end-ab_OL$start
  ab_OL$width_anno = ab_OL$end_anno-ab_OL$start_anno
  
  ###phil what do we do if everything is 0? 
  ## That should not haappen if all 0 then the file would be empty and there would be no overlap features.
  #subset 
  ab_OL=ab_OL[ab_OL$ntOL>0,]
  if (nrow(ab_OL)>0){
    #add nt / width in annotation and read
    ab_OL$OLper_anno = ab_OL$ntOL/ab_OL$width_anno
    ab_OL$OLper = ab_OL$ntOL/ab_OL$width
    
    #subset
    # do not use bedtools -F or -f cut off because will not catch peaks that are larger than annotation (because stitched togheter to form wider peak) 
    # eg is rRNA subunit that was missed by -f flag
    ab_OL=ab_OL[(ab_OL$OLper_anno>.75 | ab_OL$OLper>.51), ]
    
    #ID duplicates and unique
    dup=unique(ab_OL[duplicated(ab_OL$ID),'ID'])
    ab_OL_single=ab_OL[!(ab_OL$ID%in%dup),]
    ab_OL_double=ab_OL[(ab_OL$ID%in%dup),]
    
    unique_list=unique(ab_OL_double$ID)
    ab_OL_colapsed=as.data.frame(matrix(nrow=length(unique_list),ncol=ncol(ab_OL_double)));
    colnames(ab_OL_colapsed)=colnames(ab_OL_double)
    
    #for every unique id
    for(x in 1:length(unique_list)){
      peaksTableU = unique_list[x]
      pam = ab_OL_double[ab_OL_double$ID%in%peaksTableU,]
      ab_OL_colapsed[x,colnames(ab_OL_double)%in%c(ColumnName)==F]=(pam[1,colnames(ab_OL_double)%in%c(ColumnName)==F,drop=T])
      
      for (cx in 1:length(ColumnName)) {
        ab_OL_colapsed[x,ColumnName[cx]]=paste(sort(unique((pam[,ColumnName[cx]]))), collapse =",")
      }
    }
    
    ab_OL_colapsed=rbind(ab_OL_colapsed,ab_OL_single)
    ab_OL_colapsed=ab_OL_colapsed[!is.na(ab_OL_colapsed$ID),]
    ab_OL_colapsed[(ab_OL_colapsed$ntOL<=0),ColumnName]=NA
    
    ###fix rename variable
    ab_OL=ab_OL_colapsed
    print(paste0('Overlaps with ',nrow(ab_OL),' Regions from '))[1]
  } else{
    print("does not overlap with ")[1]
    # print("DO SOMETHING")
  }
  
  print((match.call())$Annotable)
  
  ab_OL = merge(peaksTable,
                       ab_OL[,c('ID',ColumnName)],
                       by='ID',all.x=T)
  
  
  return(ab_OL) 
}

#covers 1123 through 1342, 1397 through 13404
peak_calling<-function(peak_in,nmeprfix){
  if(length(args)==0){
    peak_in=peaks
    nmeprfix='Same_'
    }
  
  print((match.call())$peak_in)
  
  colMerge = c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')
  colSelect=c('ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')
  
  # if (ref_species=='hg38'){
  #   PeaksdataOut=bam_anno2(peak_in,
  #               rbind(ref_gencode_trans[,colMerge],
  #                     # ref_refseq[,colMerge],
  #                     ref_lncra[,colMerge]), 
  #               colSelect,
  #               "pass1")
  # } else if (ref_species=='mm10'){
  # 
  #   
  #   PeaksdataOut=bam_anno2(peak_in,
  #               ref_gencode_trans[,colMerge],
  #               colSelect,
  #               "pass1")
  # }
  
  PeaksdataOut=bam_anno2(peak_in,
                         gencode_Anno_RNA_comb[,colMerge],
                         colSelect,
                         "pass1")
  
  #merge peaks_jucntion with p file just created
  ###phil isn't this is redundant? 
  ## fixed in bam_anno2 function
  # PeaksdataOut = merge(peaks,
  #                      p[,c('ID','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')],
  #                      by='ID',all.x=T)
  
  # PeaksdataOut=bam_anno('RNA_anno',Anno_RNA_comb,PeaksdataOut)
  ### fix variable naming is confusing - back and forth between p and PeaksdataOut
  ###phil Same as above... seems redundant especially when merging back?
  ### I kept seperate so that we could potentially make this section more flexible if there are additional ncRNA annotation files
  ### the idea is first use the general gencode annnotations then use specific ncRNA annotations
  PeaksdataOut=bam_anno2(PeaksdataOut,
              Anno_RNA_comb,
              c('type','transcript_name'),
              "pass2")
 
  
  ##########################################################################################
  ############### Clean up
  ##########################################################################################
  # change lincRNA,rRNA to rRNA only
  # do not change the name if there is actually a lincRNA name separate from the rRNA grepl(',',p$name)==F)
  PeaksdataOut[(grepl(',',PeaksdataOut$name)==F)&(PeaksdataOut$type%in%'lincRNA-exon,rRNA'),'type']='rRNA'
  PeaksdataOut[(grepl(',',PeaksdataOut$name)==F)&(PeaksdataOut$type%in%'lincRNA-intron,rRNA'),'type']='rRNA'
  PeaksdataOut[(grepl(',',PeaksdataOut$name)==F)&(PeaksdataOut$type%in%'lncRNA,rRNA'),'type']='rRNA'
  
  # change all psueudogene classes to pseudogene
  PeaksdataOut[grep('pseudogene',PeaksdataOut$gene_type_ALL),'gene_type_ALL']='pseudogene'
  
  # change all sc RNA to yRNA names
  PeaksdataOut$type=gsub('scRNA,yRNA','yRNA',PeaksdataOut$type)
  
  
  
  ##########################################################################################
  ############### duplicates / unique 
  ##########################################################################################
  #handle duplicates
  dup=PeaksdataOut[(is.na(PeaksdataOut$gene_type_ALL)==F)|(is.na(PeaksdataOut$type)==F),'ID']
  peaksTable_single=PeaksdataOut[!(PeaksdataOut$ID%in%dup),]
  peaksTable_double=PeaksdataOut[(PeaksdataOut$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  #for each unique id create...?
  ###phil help explain what is being created here
  ## Checiking for conflicting annotations from gencode and additional supplied annotations
  ## If conflict prioritize Supplied annotations 
  ## Additionally Peak may overlap multiple features and so I am setting hierarchy
  
  for(x in 1:length(u)){
    PeaksdataOutU=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%PeaksdataOutU,]
    
    #Comb Annotation
    pam_1=pam$gene_type_ALL ### type from Gencode 
    pam_1=as.data.frame(strsplit(pam_1,','));colnames(pam_1)='a';pam_1$a=as.character(pam_1$a)
    
    pam_2=pam$type ### type from Annotation
    pam_2=as.data.frame(strsplit(pam_2,','));colnames(pam_2)='a';pam_2$a=as.character(pam_2$a)
    
    ###phil - if we are doing this for all cases, we should remove them in the anno script
    #remove misc_RNA catagory : these are covered by additional annotations (mostly yRNA)
    ## We should set this up to only display misc_RNA if no additional annotations
    if (grepl('misc_RNA',pam_1)&is.na(pam_2)==F) {pam_1[pam_1$a%in%'misc_RNA','a']=NA}
    
    ###phil - Same here - this should all be handled in anno script
    #lincRNA are from an older Gencode version (VM18 or older) so don't use    
    if (grepl('lncRNA',pam_1)&grepl('lincRNA',pam_2)) {pam_1[pam_2$a%in%'lincRNA',]='lncRNA'}
    
    ###phil Same as above
    #change ribozyme (gencode) to RNA type from additional anno
    if (grepl('ribozyme',pam_1)) {pam_1[pam_1$a%in%'ribozyme',]=unique(pam_2$a)}
    
    ###phil - my hope would be all of this would be unnnecessary; we should have the
    #annotations set up to be exactly as we want to avoid having to do so much manipulation
    #to the actual peak data before adding on the annotation
    #combine all rna types subtypes into ncRNA
    pam_c=rbind(pam_1,pam_2)
    pam_c=pam_c[!is.na(pam_c$a),,drop=F]
    
    #### Annotation from Gencode : unique(sort(mm10$gene_type_ALL))
    ReplaceType <-function(df_in,col_in,list_in,replace_id){
      for (id in list_in){
        df_in[,col_in]=gsub(id,replace_id,df_in[,col_in])
      } 
      return(df_in)
    }
    
    type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA",
                  "ribozyme","scRNA", "sRNA","scaRNA", "yRNA", "srpRNA", "tRNA", "7SKRNA",
                  "sncRNA")
    pam_c2 = ReplaceType(pam_c,"a", type_list,"ncRNA")
    
    ###phil what if the length isn't greater than 1? keep the full file? 
    ## yes keep full file
    if (length(pam_c$a)>1) {pam_c=pam_c[order(pam_c$a),,drop=F]}
    if (length(pam_c2$a)>1) {pam_c2=pam_c2[order(pam_c2$a),,drop=F]}
    
    pam_c=paste(unique(pam_c$a),collapse = ',')
    pam_c2=paste(unique(pam_c2$a),collapse = ',')
    
    #Comb GeneName 
    pam_G1=pam$external_gene_name ### gene name from gencode
    pam_G1=as.data.frame(strsplit(pam_G1,','));colnames(pam_G1)='a';pam_G1$a=as.character(pam_G1$a)
    
    ###phil this is confusing - what do annotation are we using? isngt anno a combo from gencode?
    ### G1 is gencode annotation and G2 is additional annotation
    pam_G2=pam$transcript_name ### gene name from annotation
    pam_G2=as.data.frame(strsplit(pam_G2,','));colnames(pam_G2)='a';pam_G2$a=as.character(pam_G2$a)
    
    ###phil is this for the Same reason as the mirna above?
    ## there are non uniform naming for 7SK rRNA so i am making a uniform naming by selecting only from gencode
    if (grepl('sk',pam_G1)&grepl('7SK',pam_G2)) {xxxx
      pam_G2[pam_G2$a%in%'7SK',]=NA
    }
    
    #merge gene names from "gencode" and "annotation"
    pam_gmerge=rbind(pam_G1,pam_G2)
    pam_gmerge=pam_gmerge[!is.na(pam_gmerge$a),,drop=F]
    pam_gmerge=paste(unique(pam_gmerge$a),collapse = ',')
    
    ### Create Table 
    # ###phil can we just remove these cols? 
    # ## we use them later on in final table
    # PeaksdataOut$type_simple_comb=NA
    # PeaksdataOut$type_comb=NA
    # PeaksdataOut$gene_name_comb=NA
    
    peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
    peaksTable_colapsed[x,'type_simple_comb']=pam_c2
    peaksTable_colapsed[x,'type_comb']=pam_c
    peaksTable_colapsed[x,'gene_name_comb']=pam_gmerge
    
    #remove('pam_1','pam_2','pam_c','pam_c2','pam_G1','pam_G2','pam_gmerge','pam')
  }
  
  ###fix p variable swap
  rnatype=peaksTable_colapsed[,c('ID','ensembl_gene_id','external_gene_name','gene_type',"gene_type_ALL",
               'type','transcript_name','type_simple_comb','type_comb','gene_name_comb')]
  
  
  peaksTable_single[,c('type_simple_comb','type_comb','gene_name_comb')]=NA
 
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  rnames=c("ensembl_gene_id","external_gene_name","gene_type","gene_type_ALL","type","transcript_name","type_simple_comb","type_comb","gene_name_comb")
  
  colnames(peaksTable_colapsed)[colnames(peaksTable_colapsed)%in%rnames]=
    paste0(nmeprfix,rnames)
  

  
  type_start = c('lincRNA-exon,lncRNA','lincRNA-intron,lncRNA')
  type_end = c('lincRNA-exon','lincRNA-intron')
  col_list= c('type_comb','type_simple_comb','gene_type_ALL')
    #change the type (type_list) to new id (change_list)

  for (cols_to_change in col_list){
    
    #change the type (type_list) to new id (change_list)
    for (i in range(1:length(type_start))){
      peaksTable_colapsed[,paste0(nmeprfix,cols_to_change)] = gsub(type_start[i],
                                                        type_end[i],
                                                        peaksTable_colapsed[,paste0(nmeprfix,cols_to_change)] )
    }
  }
  ###fix variable switch
  ###phil this is never used again
  # PeaksdataOut=peaksTable_colapsed
  
  ###phil why are we creating NA cols?
  ## no longer applies
  # col_list = c("feature","exon_number","exon_LargeOL","intron_number","intron_LargeOL",
  #              "intron_5pStart","intron_length","exon_length")
  # for (col_names in col_list){
  #   peaksTable[,paste0(nmeprfix,col_names)]=NA
  # }
  
  # return(peaksTable)
  return(peaksTable_colapsed)
}

### fix designations - 1 = "Same" 2 = "Oppo"? Same strand and complementary strand?
#can prob  be more descriptive
peaks_SameAnno = peak_calling(peaks,"Same_")
peaks_OppoAnno = peak_calling(peaks_Oppo,"Oppo_")

##########################################################################################
##############################  INTRON EXON ANNOTATION
##########################################################################################
#lines 1368 through 1395, 1406 through 1490
### Identify if CLIP peak overlaps with Intron or Exonic region   

#   Using GTF file from GENCODE v23 -mm10
#   Using GTF file from GENCODE v32 -hg38

# Peaks were annotated by whether they overlap with Host gene intron/exon region
# Intron coordinates were calculated from GTF file.
# 
# A second column was added to idenify if the peak also overlapped with the 5'UTR 3'UTR or CDS (Column: Featrue 2)
#########################################################################################3
#### Sam - I think we should create Exon_intron file here as well and not just in 06. This is repetitive but it will be a really big file and i do not think we should save it for each dataset. 
ref_gencode = fread(gencode_path, header=T, sep="\t",stringsAsFactors = F,data.table=F)

#dplyr::rename cols
ref_gencode = ref_gencode %>% 
  dplyr::rename(
    chr = seqname,
    ensembl_gene_id = gene_id,
    external_gene_name = gene_name
  )

ref_gencode_e = subset( ref_gencode, feature == "exon")

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
# Same for exon number the exon # starts at 0
introns$start=introns$start+1
introns$exon_number=as.numeric(introns$exon_number)+1

## add gene type type
introns=merge(introns,
              (unique(ref_gencode_trans[,c('transcript_id','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL','transcript_type','transcript_name','score')])),by='transcript_id',all.x=T)

#############################3
############### binding exon, intron

intron_exon=rbind(ref_gencode_e[,c('chr','feature','start','end','strand','transcript_id','ensembl_gene_id','external_gene_name','gene_type','transcript_type','transcript_name','score','exon_number')],
                  introns[,c('chr','feature','start','end','strand','transcript_id','ensembl_gene_id','external_gene_name','gene_type','transcript_type','transcript_name','score','exon_number')])
intron_exon$ID=paste0(intron_exon$chr,':',intron_exon$start,'-',intron_exon$end)


##########################################################################################
############### Intron Exon Function
##########################################################################################

IE_calling <- function(peak_in,Peak_Strand_ref,nmeprfix){
  ###phil designations - 1 = "Same" 2 = "Oppo"? Same strand and complementary strand?
  #xopp = 1 #Same
  #xopp = 2 #Oppo
  
  if(length(args)==0){
    peak_in=peaks_OppoAnno
    Peak_Strand_ref='Oppo'
    nmeprfix='Oppo_'
  }
  
  # print((match.call())$peak_in)
  
  
  #create annotation table
  ColumnName = c("feature","exon_number")
  anno_IntExn = intron_exon[grep('protein_coding',intron_exon$gene_type),] #%>%
    # select(c('chr','start','end','transcript_id','transcript_id','strand',ColumnName))
  anno_IntExn=anno_IntExn[,c('chr','start','end','transcript_id','transcript_id','strand',ColumnName)]
  # colnames(anno_tmp)[colnames(anno_tmp)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')

  ################################################################################################################################  
  #
  # #run bam_anno 
  # exoninof = bam_anno2(peak_in[,c('chr','start','end','ID','ID2','strand')],
  #                      anno_IntExn,
  #                      c('feature','exon_number'),
  #                      "pass3")
  # 
################################################################################################################################  
  ##phil was there a reason for not using bam_anno2?
  ## this is a little different from bam_anno2 because trying to calculate distance to feature and then collapsing to each feature
  

  ##phil col name is being re-written
  ColumnName=paste0(c('ensembl_gene_id','external_gene_name'))
  ColumnName=c("feature","exon_number")

  ###fix variable swap
  p=peak_in[,c('chr','start','end','ID','ID2','strand')]
  a=anno_IntExn[,c('chr','start','end','transcript_id','transcript_id','strand',ColumnName)]

  ###phil why are we clearing this column
  a[,5]=""
  colnames(a)[colnames(a)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')

  ###phil - why are we writing the Same files over again?
  write.table(a,file=paste0(out_dir,"/","/annotable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  write.table(p,file=paste0(out_dir,"/","/peakstable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)

  system(paste0('bedtools intersect -a ',gsub(" ","\\\\ ",out_dir),'/','/peakstable.bed -b ',gsub(" ","\\\\ ",out_dir),'/','/annotable.bed -wao -s  >',gsub(" ","\\\\ ",out_dir),'/','/peaks_OL.txt'))

  ### fix - writing Same files again
  ##phil are there differences between this and the bam_anno2 function? If not we can run together
  exoninof=fread(paste0(out_dir,"/peaks_OL.txt"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(exoninof)=c(paste0(colnames(p)),paste0(colnames(a)),'ntOL')
  exoninof=exoninof[exoninof$ntOL>0,]
  exoninof$width_anno=exoninof$end_anno-exoninof$start_anno
  exoninof$OLper_anno=exoninof$ntOL/exoninof$width_anno
  exoninof$width=exoninof$end-exoninof$start
  exoninof$OLper=exoninof$ntOL/exoninof$width

  exoninof=exoninof[(exoninof$OLper_anno>.75 | exoninof$OLper>.51), ]
  ################################################################################################################################  
  
  ###phil what if this is opp? we ignore?
  ## yes only doing distance calculation for Same strand features
  if (Peak_Strand_ref=='Same') {
    #Calc distance of 5' peak to 5' intron/exondist
    exoninof_pos=exoninof[exoninof$strand%in%'+',]
    exoninof_neg=exoninof[exoninof$strand%in%'-',]
    
    #run positive
    exoninof_pos[,paste0(nmeprfix,'feature_Distance')]=((exoninof_pos$start-exoninof_pos$start_anno)/abs(exoninof_pos$end_anno-exoninof_pos$start_anno))*100
    exoninof_pos[,paste0(nmeprfix,'feature_5pStart')]=exoninof_pos$start_anno
    exoninof_pos[,paste0(nmeprfix,'feature_length')]=abs(exoninof_pos$end_anno-exoninof_pos$start_anno)
    
    #run negative
    exoninof_neg[,paste0(nmeprfix,'feature_Distance')]=((exoninof_neg$end_anno-exoninof_neg$end)/abs(exoninof_neg$end_anno-exoninof_neg$start_anno))*100
    exoninof_neg[,paste0(nmeprfix,'feature_5pStart')]=exoninof_neg$end_anno
    exoninof_neg[,paste0(nmeprfix,'feature_length')]=abs(exoninof_neg$end_anno-exoninof_neg$start_anno)
    
    #bind positive and negative
    exoninof=rbind(exoninof_pos,exoninof_neg);
    
    ###phil why are we pushing NA cols?
    peak_in[,paste0(nmeprfix,c('Exn_start_dist','Intron_start_dist','Intron_5pStart','Exn_5pStart'))]=NA
  }
  
  #for each peak in table
  for (x in 1:nrow(peak_in)) {
    
    #subset ids
    l=peak_in[x,c('ID','ID2')]
    g=exoninof[exoninof$ID%in%l$ID,]
    
    #if there are ids in both lists
    ###fix we can handle this better
    if (nrow(g)>0) {
      
      #split into exon and intron
      g_e=g[g$feature%in%'exon',]
      g_i=g[g$feature%in%'intron',]
      
      
      gname=paste(unique(g$feature),collapse = ",")
      gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
      peak_in[x,paste0(nmeprfix,'feature')]=gsub('NA,',"",gname)
      
      ###phil we need to handle all the else statements
      #if there are exons
      if (nrow(g_e)>0) {
        gname=paste(unique(g_e$exon_number),collapse = ",")
        gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
        peak_in[x,paste0(nmeprfix,'exon_number')]=gsub('NA,',"",gname)
        
        peak_in[x,paste0(nmeprfix,'exon_LargeOL')]=max(g_e$ntOL/(g_e$end-g_e$start))
        
        #if the Same
        if (Peak_Strand_ref=='Same') {## Do not need feature distance for Opposite strand features
          g_e[g_e[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
          peak_in[x,paste0(nmeprfix,'Exn_start_dist')]=mean(as.numeric(g_e[,paste0(nmeprfix,'feature_Distance')]))
        } 
      } #else{
      #  print("Not Exonic Region")
      #}
      ###fix create function to loop both sets
      if (nrow(g_i)>0) {
        gname=paste(unique(g_i$exon_number),collapse = ",")
        gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
        peak_in[x,paste0(nmeprfix,'intron_number')]=gsub('NA,',"",gname)
        
        peak_in[x,paste0(nmeprfix,'intron_LargeOL')]=max(g_i$ntOL/(g_i$end-g_i$start))
        
        if (Peak_Strand_ref=='Same') {
          g_i[g_i[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
          peak_in[x,paste0(nmeprfix,'Intron_start_dist')]=mean(as.numeric(g_i[,paste0(nmeprfix,'feature_Distance')]))
        }
      } #else{
      #  print ("Not Intronic region")
      #}
      ###phil why do we have different cols depending on the strand - can we make these the Same?
      ### We are not calculating distance to feature on Opposite strand
    } else {
      if (Peak_Strand_ref=='Same') {peak_in[x,paste0(nmeprfix,c("feature",'exon_number','intron_number','Exn_start_dist','Intron_start_dist'))]=NA}
      if (Peak_Strand_ref=='Oppo') {peak_in[x,paste0(nmeprfix,c("feature",'exon_number','intron_number'))]=NA}
    } 
  }
  
  

  ###phil what is this fixing?
  ### protien coding peaks with not Correctly assigned exon overlap
  peak_in[is.na(peak_in[,paste0(nmeprfix,'feature')]) & 
               peak_in[,paste0(nmeprfix,'gene_type')]%in%'protein_coding',paste0(nmeprfix,'feature')]='exon'
  peaksTable_inex=peak_in[grep(',',peak_in[,paste0(nmeprfix,'feature')]),]
  
  return(peak_in)
}

intronexon_Same = IE_calling(peaks_SameAnno,'Same',"Same_")
intronexon_Opposite = IE_calling(peaks_OppoAnno,'Oppo',"Oppo_")

##########################################################################################
##############################  IDENTIFY PEAKS IN REPEAT REGIONS   
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
#read in rmsk
rmsk_GRCm38=fread(rmsk_path, header=T, sep="\t",stringsAsFactors = F,
                  data.table=F)

rpmsk_anno <- function(ColumnName,Annotable,peaksTable){
  
  if(length(args)==0){
    nmeprfix="Same_"
    peaksTable=intronexon_Same
    Annotable=rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('LINE','SINE'),]
    ColumnName=paste0(nmeprfix,'Repeat_LINE_SINE')
  }
  
  ###fix variable

  #create df of ids and start location
  peaksTable=separate(peaksTable,ID,into=c('chr','start'),sep=":",remove=F)
  peaksTable=separate(peaksTable,start,into=c('start','end'),sep="-",remove=F)
  
  #create granges obj for peaks
  peaksTable.GR <- GRanges(seqnames = as.character(peaksTable$chr), 
                  ranges=IRanges(start = as.numeric(peaksTable$start),
                                 end = as.numeric(peaksTable$end)),
                  strand = peaksTable$strand,ID=peaksTable$ID )
  
  #create granges obj for annotation
  anno.GR <- GRanges(seqnames = as.character(Annotable$genoName), 
                     ranges=IRanges(start = as.numeric(Annotable$genoStart), 
                                    end = as.numeric(Annotable$genoEnd)),
                     strand = Annotable$strand,repClass=Annotable$repClass,repName=Annotable$repName )
  
  ###fix variable renames
  q =peaksTable.GR
  s=anno.GR
  
  #find annotation overlaps
  xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))
  
  #run query
  qh=as.data.frame(q[xo$queryHits],row.names = NULL)
  sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_repeat")
  rmskinfo=cbind(qh,sh)
  rmskinfo=rmskinfo[,c('ID','repClass_repeat')];colnames(rmskinfo)[colnames(rmskinfo)%in%'repClass_repeat']= ColumnName      
  
  ###fix variable
  #merge rmskinfo with peakstable
  peaksTable=merge(peaksTable[,!colnames(peaksTable)%in%ColumnName],rmskinfo,by='ID',all.x=T)
  
  #check dups and unique peaks
  dup=unique(peaksTable[duplicated(peaksTable$ID),'ID'])
  peaksTable_single=peaksTable[!(peaksTable$ID%in%dup),]
  
  peaksTable_double=peaksTable[(peaksTable$ID%in%dup),]
  u=unique(peaksTable_double$ID)
  
  #collapse duplicate peaks
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)))
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  #for each unique peak
  for(x in 1:length(u)){
    ###fix variable
    ###fix code similarly to above - function?
    peaksTable = u[x]
    pam = peaksTable_double[peaksTable_double$ID%in%peaksTable,]
    
    #if the cols of double are not in ColumnName list
    ###phil what is the purpose of this? 
    peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F] = (pam[1,colnames(peaksTable_double)%in%c(ColumnName)==F,drop=T])
    
    #past the uqique values of the col to collapsed table
    peaksTable_colapsed[x,ColumnName]=paste(sort(unique((pam[,ColumnName]))), collapse =" | ")
  }
  
  #bind collapsed and unqiue
  peaksTable_colapsed = rbind(peaksTable_colapsed,peaksTable_single)
  
  #if there is no ID remove
  peaksTable_colapsed = peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  ###phil why are we adding NA columns
  peaksTable_colapsed[is.na(peaksTable_colapsed[,ColumnName]),ColumnName]=NA
  peaksTable_colapsed[peaksTable_colapsed[,ColumnName]%in%"",ColumnName]=NA
  
  
  return(peaksTable_colapsed) 
}

rpmsk_calling <- function(peaks_in,nmeprfix){
  ###fix variable names
  if(length(args)==0){
    nmeprfix="Same_"
    peaks_in=intronexon_Same
  }
  
  PeaksdataOut=peaks_in[!is.na(peaks_in$ID),]
  
  rpmsk_list = c("LTR","DNA","Satellite","Simple_repeat","Low_complexity","Other","Unknown")
  rpmsk_names = c("Repeat_LTR","Repeat_DNA","Repeat_Satalites","Repeat_Simple_Repeats",
                  "Repeat_Low_Complexity","Repeat_Other","Repeat_Unknown")
  i = 1
  for(find_type in rpmsk_list){
    rmsk_subset = subset(rmsk_GRCm38, repClass == find_type)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,rpmsk_names[i]),rmsk_subset,PeaksdataOut)
    i = i+1
  }
  
  #run with final list
  PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LINE_SINE'),
                          rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('LINE','SINE'),],
                          PeaksdataOut)
  
  ###fix variable naming  
  p=PeaksdataOut
  p[,paste0(nmeprfix,'Repeat_comb')]=NA
  repcol=paste0(nmeprfix,c('Repeat_LINE_SINE','Repeat_LTR','Repeat_DNA',
                           'Repeat_Satalites','Repeat_Simple_Repeats','Repeat_Low_Complexity',
                           'Repeat_Other','Repeat_Unknown'))
  
  dup=p[rowSums((p[,repcol]>0),na.rm = T)>0,'ID']
  
  #run duplicated and unique
  peaksTable_single=p[!(p$ID%in%dup),]
  peaksTable_double=p[(p$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  #for each unique id
  for(x in 1:length(u)){
    p=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    pam_c=((pam[1,repcol,drop=F]))
    
    pmat=as.data.frame(matrix(nrow=length(pam_c),ncol=1))
    
    ###fix col name 
    colnames(pmat)='a'
    pmat$a=t(pam_c)
    
    #remove na's
    pmat=pmat[is.na(pmat)==F,]
    
    #paste uniques
    pmat=paste(unique(pmat),collapse = ',')
    
    ###phil what is this doing?
    peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
    peaksTable_colapsed[x,paste0(nmeprfix,'Repeat_comb')]=pmat
    
  }
  
  #bind collapsed with single
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  
  #remove na's
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  ###fix variable naming
  p=peaksTable_colapsed
  PeaksdataOut=peaksTable_colapsed
  
  return(PeaksdataOut)
}

rpmsk_Same = rpmsk_calling(intronexon_Same,"Same_")
rpmsk_Opposite = rpmsk_calling(intronexon_Opposite,"Oppo_")

##########################################################################################
##############################  Assigning Clip peak attributes by strand
##########################################################################################
# Not all Peaks overlap with a single feature so peak assignments were assigned by priority:
# ncRNA > Protein coding : Exonic > repeats > Pseudogene > Antisense Feature > 
#Protein Coding : Intronic > lncRNA > no Feature  
# All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table

peak_attributes <- function(peaks_in,nmeprfix){
  if(length(args)==0){
    nmeprfix="Oppo_"
    peaks_in=rpmsk_Opposite
    # colnames(peaks_in)
  }  
  #convert lincRNA and lncRNA types
  col_list = c("type_simple_comb","type_comb","gene_type","gene_type_ALL")
  type_start = c("lincRNA","lncRNA")
  type_end = c("linLcRNA","lnLcRNA")
  
  #for each col selected -
  for (cols_to_change in col_list){
    
    #change the type (type_list) to new id (change_list)
    for (i in range(1:length(type_start))){
      peaks_in[,paste0(nmeprfix,cols_to_change)] = gsub(type_start[i],
                                                            type_end[i],
                                                            peaks_in[,paste0(nmeprfix,cols_to_change)] )
    }
  }
  
  #replaced above
  #PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')] )
  #PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')] )
  
  #PeaksdataOut[,paste0(nmeprfix,'type_comb')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_comb')] )
  #PeaksdataOut[,paste0(nmeprfix,'type_comb')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_comb')] )
  
  #PeaksdataOut[,paste0(nmeprfix,'gene_type')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type')] )
  #PeaksdataOut[,paste0(nmeprfix,'gene_type')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type')] )
  
  #PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')] )
  #PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')] )   
  
  AddTypetoComb <- function(df_in,subset_df,type_id){
    df_in[subset_df,paste0(nmeprfix,'Comb_type_exon')]=type_id
    return(df_in)
  }
  
  #create column
  peaks_in[,paste0(nmeprfix,'Comb_type_exon')]="TBD"
  
  # 1. ncRNA
  peaks_in = AddTypetoComb(peaks_in,
                               grepl('ncRNA',peaks_in[,paste0(nmeprfix,'type_simple_comb')]),
                               "ncRNA")
  
  # 2. protein coding - Exonic
  peaks_in = AddTypetoComb(peaks_in,
                                        peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD" & 
                                        grepl('protein_coding',peaks_in[,paste0(nmeprfix,'type_simple_comb')]) & 
                                        grepl('exon',peaks_in[,paste0(nmeprfix,'feature')]),
                               "protein_coding: exon")

  # 3. repeats
  ###phil check the repeat_comb is correct
  peaks_in = AddTypetoComb(peaks_in,
                               peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD" & 
                                 is.na(peaks_in[,paste0(nmeprfix,'Repeat_comb')])==F,
                               "Repeat Element")

    # 4. Pseudogene
  peaks_in = AddTypetoComb(peaks_in,
                           peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD" &
                             grepl('pseudogene',peaks_in[,paste0(nmeprfix,'type_simple_comb')]),
                               "pseudogene")
  
  # 5. lncRNA-exon
  peaks_in = AddTypetoComb(peaks_in,
                           peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD" &
                           (grepl('linLcRNA-exon',peaks_in[,paste0(nmeprfix,'type_comb')]) |
                               grepl('lnLcRNA',peaks_in[,paste0(nmeprfix,'type_comb')]) ),
                               "lnLcRNA-exon")

  # 6. intron
  peaks_in = AddTypetoComb(peaks_in,
                           peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD" &
                             grepl('protein_coding',peaks_in[,paste0(nmeprfix,'type_simple_comb')]),
                               "protein_coding: Intron")

  # 7. lncRNA-intron
  peaks_in = AddTypetoComb(peaks_in,
                           peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD" &
                             grepl('linLcRNA-intron',peaks_in[,paste0(nmeprfix,'type_comb')]),
                               "lnLcRNA-intron")

  #all other features
  peaks_in = AddTypetoComb(peaks_in,
                           peaks_in[,paste0(nmeprfix,'Comb_type_exon')] == "TBD",
                               "no Feature")

  #create factor of col
  peaks_in[,paste0(nmeprfix,'Comb_type_exon')] = 
    factor(peaks_in[,paste0(nmeprfix,'Comb_type_exon')], 
           levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene",
                      "lnLcRNA-exon","Antisense Feature","protein_coding: Intron","lnLcRNA-intron","no Feature"))
  
  #replaced above
  #####  
  # # 1. ncRNA        
  # comp = ( grepl('ncRNA',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) ) 
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='ncRNA'
  # 
  # # 2. protein coding - Exonic
  # PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]=NA
  # comp=(is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]) &
  #         grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) &
  #         (grepl('exon',PeaksdataOut[,paste0(nmeprfix,'feature')])))
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='protein_coding: exon'
  
  # 
  # # 3. repeats
  # comp =(is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]) &
  #          (is.na(PeaksdataOut[,paste0(nmeprfix,'Repeat_comb')])==F) )
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]="Repeat Element"
   
  # 
  # # 4. Pseudogene
  # comp = (is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])
  #         & grepl('pseudogene',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]))
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='pseudogene'
  # 
  # # 5. lncRNA-exon
  # comp1 = is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])  
  # comp2 = grepl('linLcRNA-exon',PeaksdataOut[,paste0(nmeprfix,'type_comb')]) |
  #   grepl('lnLcRNA',PeaksdataOut[,paste0(nmeprfix,'type_comb')])
  # comp = comp1&comp2       
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='lnLcRNA-exon'
  # 
  # # 6. intron
  # comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]) &
  #          grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='protein_coding: Intron'

  # 
  # # 7. lncRNA-intron
  # comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]) &
  #          (grepl('linLcRNA-intron',PeaksdataOut[,paste0(nmeprfix,'type_comb')]) ))
  # PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='lnLcRNA-intron'
  
  # 
  # 
  # #all other features
  # PeaksdataOut[is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]),paste0(nmeprfix,'Comb_type_exon')]='no Feature'
  # PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')] = 
  #   factor(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')], 
  #          levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene",
  #                     "lnLcRNA-exon","Antisense Feature","protein_coding: Intron","lnLcRNA-intron","no Feature"))
  #####
  
  #create col
  peaks_in[,paste0(nmeprfix,'Comb_type_ncRNA')]=NA
  
  #seleect ncRNA
  ###fix variable
  ## subset was not working correctly
  p1=peaks_in[peaks_in[,paste0(nmeprfix,'Comb_type_exon')]=='ncRNA',]
  
  #make col type_comb with prefix  
  p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=p1[,paste0(nmeprfix,'type_comb')]
  
  ###phil why are we removing all of theses types? 
  #remove types
  ### Protein coding + ncRNA annotations -> ncRNA subtype
  ### any ncRNA takes priority over pseudogene
  ### any double annotations with lncRNA become second annotation only
  type_list = c("protein_coding,",",protein_coding","pseudogene,",",pseudogene",'lnLcRNA,')
  for (type_id in type_list){
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(type_id,"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  }  
  
  #merge linLcRNA  
  p1[p1[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-intron,linLcRNA-exon',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
  p1[p1[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-exon,linLcRNA-intron',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
  
  type_list = c("lnLcRNA,",'linLcRNA,','linLcRNA-exon,',"linLcRNA-intron,")
  for (type_id in type_list){
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(type_id,"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
  } 
  
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
  
  #select not ncRNA
  p2=peaks_in[!peaks_in[,paste0(nmeprfix,'Comb_type_exon')]%in%'ncRNA',]
  
  ###fix variable
  PeaksdataOut=rbind(p1,p2)
  
  return(PeaksdataOut)
}

peak_attrib_Same = peak_attributes(rpmsk_Same,"Same_")
peak_attrib_Oppo = peak_attributes(rpmsk_Opposite,"Oppo_")

###fix variable swap
#merge attributes
PeaksdataOut = merge(peak_attrib_Same,
                     peak_attrib_Oppo[,colnames(peak_attrib_Oppo)[!colnames(peak_attrib_Oppo) %in% 
                                                                    c("chr","start","end","strand","ID2" )]],
                     by='ID')
##########################################################################################
##############################  Assigning Clip peak attributes - merged strands
##########################################################################################  
#create column
PeaksdataOut$Comb_type_exon_Oppo="no Feature"

#select subset of df, input specific type
AddTypetoComb <- function(df_in,subset_df,type_id){
  df_in[subset_df,paste0('Comb_type_exon_Oppo')]=type_id
  return(df_in)
}

###phil - i tried to make this more readable but not sure if it is; we can keep the old
#way if you prefer
###fix variable naming
# 1. ncRNA        
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             grepl('ncRNA',PeaksdataOut[,paste0('Same_','type_simple_comb')],fixed = T), 
                             'ncRNA')

# 2. protein coding - Exonic
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')]) & 
                               grepl('exon',PeaksdataOut[,paste0('Same_','feature')]),
                             'protein_coding: exon')

# 3. repeats
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" &
                               (is.na(PeaksdataOut[,paste0('Same_','Repeat_comb')])==F),
                             "Repeat Element")

# 4. Pseudogene
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               grepl('pseudogene',PeaksdataOut[,paste0('Same_','type_simple_comb')]),
                             "pseudogene")

# 5. Antisense - non LncRNA
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) &
                               grepl('lnLcRNA',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F,
                             "Antisense Feature")

# 6. lncRNA-exon
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               (grepl('linLcRNA-exon',PeaksdataOut[,paste0('Same_','type_comb')]) | 
                                grepl('lnLcRNA',PeaksdataOut[,paste0('Same_','type_comb')]) ),
                             'lnLcRNA-exon')
#####
# 7. Antisense - LncRNA - exon
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F & 
                               grepl('intron',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F ),
                             'Antisense Feature')

# 8. intron
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')]),
                             'protein_coding: Intron')

# 9. Antisense - LncRNA - intron
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & 
                               (grepl('linLcRNA-exon',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F),
                             'Antisense Feature')



# 10. lncRNA -intron
PeaksdataOut = AddTypetoComb(PeaksdataOut, 
                             PeaksdataOut[,'Comb_type_exon_Oppo'] == "no Feature" & 
                               grepl('linLcRNA-intron',PeaksdataOut[,paste0('Same_','type_comb')]),
                             'lnLcRNA-intron')

#replaced above

# # 1. ncRNA        
# comp=(grepl('ncRNA',PeaksdataOut[,paste0('Same_','type_simple_comb')],fixed = T) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='ncRNA'
# 
# # 2. protein coding - Exonic
# comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])&grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')])&(grepl('exon',PeaksdataOut[,paste0('Same_','feature')])) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']=paste0('protein_coding: exon')
# 
# # 3. repeats
# comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])&(is.na(PeaksdataOut[,paste0('Same_','Repeat_comb')])==F) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']="Repeat Element"
# 
# # 4. Pseudogene
# comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& grepl('pseudogene',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='pseudogene'
# 
# # 5. Antisense - non LncRNA
# comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'
###### 
# # 6. lncRNA-exon
# comp1=is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])
# comp2=((PeaksdataOut[,paste0('Same_','type_comb')]%in%c('linLcRNA-exon')) | (PeaksdataOut[,paste0('Same_','type_comb')]%in%c('lnLcRNA')) )
# comp=comp1&comp2       
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='lnLcRNA-exon'
#
# # 7. Antisense - LncRNA - exon
# comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA-intron',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'
#  
# # 8. intron
# comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
# PeaksdataOut[comp,'Comb_type_exon_Oppo']=paste0('protein_coding: Intron')
# 
# # 9. Antisense - LncRNA - intron
# comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA-exon',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'
# 
# # 10. lncRNA -intron
# comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo']) & ((PeaksdataOut[,paste0('Same_','type_comb')]%in%'linLcRNA-intron')==T) )
# PeaksdataOut[comp,'Comb_type_exon_Oppo']='lnLcRNA-intron'
# PeaksdataOut[is.na(PeaksdataOut[,'Comb_type_exon_Oppo']),'Comb_type_exon_Oppo']='no Feature'
########  
# ### only distances for final annotation_output
# ## distances to any protein coading annotation
# PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Exn_start_dist")]=NA
# PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Intron_start_dist")]=NA
# 
# ### Revert lnLcRNA to lncRNA  
# 
PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
PeaksdataOut[,'Comb_type_exon_Oppo']=factor(PeaksdataOut[,'Comb_type_exon_Oppo'], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))


# write.table(PeaksdataOut,paste0(file_id,"peakannotation.txt"),sep = "\t")
###phil - we dont need to factor this since it's the end of the script
#PeaksdataOut[,'Comb_type_exon_Oppo']=factor(PeaksdataOut[,'Comb_type_exon_Oppo'], 
# levels = c("ncRNA", "protein_coding: exon", 
#            "Repeat Element","pseudogene","lncRNA-exon",
#            "Antisense Feature","protein_coding: Intron",
#            "lncRNA-intron","lncRNA","no Feature"))

##########################################################################################
##############################  Merge peak annotations with peak junctions
##########################################################################################

Peaksdata2_anno = merge(FtrCount_splice_junc,
                        PeaksdataOut[,!colnames(PeaksdataOut)%in%c('chr','start','end','strand')],
                        by='ID')
##########################################################################################
##############################  Add annotations to junctions
##########################################################################################  

if (JoinJunc==TRUE) {
  
  #select peaks with junction
  Peaksdata2_anno_trns_exon=(Peaksdata2_anno[Peaksdata2_anno$ID%in%Junc_peaks,])
  
  #create new df
  CollapsedOut=as.data.frame(matrix(nrow = (1),ncol = ncol(Peaksdata2_anno_trns_exon)))
  colnames(CollapsedOut)=colnames(Peaksdata2_anno_trns_exon)
  CollapsedOut$IDmerge=NA
  
  #for each unique peak
  for (x in 1:nrow(PGene_TBL2)) {
    
    #select row
    Trnsc=PGene_TBL2[x,]
    
    #split ID
    TrnscUL=unlist(unique(str_split(Trnsc$ID_comb,pattern = ",")))
    
    #create df of id annotation information
    d2=Peaksdata2_anno_trns_exon[Peaksdata2_anno_trns_exon$ID%in%TrnscUL,]
    
    #if the line isn't in the annotation file
    ###phil - we should prob do some error handling here; we shouldn't have instances
    #where this is true, unless there are no peaks? and then this whole thing should stop, right? 
    if (nrow(d2)==0) {
      next
    } else if (unique(d2$strand)=="+") { if (length(unique(d2$strand)=="+")>1){ 'Join junc not strand specific'}
      d5=d2[d2$start%in%min(d2$start),] 
      dmax=d2[d2$Counts_fracMM%in%max(d2$Counts_fracMM),]
      if (nrow(dmax)>1) {dmax=dmax[dmax$start%in%min(dmax$start),]}
      #negative strand
    } else if (unique(d2$strand)=="-") {
      d5=d2[d2$start%in%max(d2$start),]
      dmax=d2[d2$Counts_fracMM%in%max(d2$Counts_fracMM),]
      if (nrow(dmax)>1) {dmax=dmax[dmax$start%in%max(dmax$start),]}
    }
    
    ### collapse all columns 
    d3 = (apply(d2 ,2, function(x){paste(unique(x[!is.na(x)]),collapse = ',')}))
    d3[d3==""]<-NA
    
    cols=c('Counts_Unique','Counts_fracMM','Length')
    d3=as.data.frame(t(d3))
    d3[1,cols] = t(as.data.frame(colSums(d2[,cols])))   
    d3[,cols] = as.numeric(d3[,cols])
    d3$IDmerge=d3$ID
    
    ## Select summary annotation based on Max(dmax) expression peak or 5'(d5) peak   
    danno=dmax # d5 dmax
    cols_AnnotSelect=c('Same_Comb_type_exon', 'Same_Comb_type_ncRNA',
                       'Oppo_Comb_type_exon','Oppo_Comb_type_ncRNA','Comb_type_exon_Oppo')
    d3[,cols_AnnotSelect]= (danno[,cols_AnnotSelect])
    
    # Select ID for 5' most read or max counts
    d3$start=d5$start
    d3$end=d5$end
    d3$ID=paste0(d5$chr,":",d5$start,"-",d5$end)   
    
    ## if count < readdepth skip
    if (d3$Counts_fracMM<read_depth) {
      next
    }
    
    CollapsedOut=rbind(CollapsedOut[,colnames(CollapsedOut)],d3[1,colnames(CollapsedOut)])
  }
  
  #remove rows without id
  CollapsedOut = CollapsedOut[complete.cases(CollapsedOut[,"ID"]),]
  
  ###fix variable we dont want to overwrite the original
  #subset for all non-junction peaks
  Peaksdata2_anno$IDmerge=NA
  Peaksdata2_anno=subset(Peaksdata2_anno, !(ID%in%Junc_peaks))
  
  #merge non-junction peaks (Peaksdata2_anno) with junction peaks (CollapsedOut)
  Peaksdata2_anno=(rbind(Peaksdata2_anno,CollapsedOut))
  
  #### Re filter reads after combing counts of spliced reads
  Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$Counts_fracMM>=read_depth,]
  
} else if (JoinJunc==FALSE) {
  
  ###phil need to fill else statement
  # collapse exon
  if (params$Condense==TRUE) {
    Peaksdata2_anno_trns_exon = subset(Peaksdata2_anno, !(Same_feature == "") &
                                         Comb_type_exon_Oppo == 'protein_coding: exon')
    #create list
    trns=Peaksdata2_anno_trns_exon$Same_ensembl_gene_id
    
    ###fix variable naming
    trns=trns[grep(",",trns)]
    trns=c(Peaksdata2_anno_trns_exon$Same_ensembl_gene_id,unlist(strsplit(trns,",")))
    
    #remove dups
    dupGeneName=unique(trns[duplicated(trns)])
    
    #create df
    CollapsedOut = as.data.frame(matrix( nrow = 1,
                                         ncol = ncol(Peaksdata2_anno_trns_exon)+1))
    colnames(CollapsedOut)=c(colnames(Peaksdata2_anno_trns_exon),'IDmerge')
    
    #for each unique val
    for (x in 1:length(dupGeneName)) {
      d1=dupGeneName[x]
      d2=Peaksdata2_anno_trns_exon[grep(d1,Peaksdata2_anno_trns_exon$Same_ensembl_gene_id),]
      d2ID=paste0(d2$chr[1],":",min(d2$start),"-",max(d2$end))
      
      #positive strand
      if (unique(d2$strand)=="+") {
        d3=d2[d2$start%in%min(d2$start),]
        #negative strand
      } else if (unique(d2$strand)=="-") {
        d3=d2[d2$start%in%max(d2$start),]
      }
      
      d3[,'start']=min(d2$start)
      d3[,'end']=max(d2$end)
      d3$IDmerge=paste(d2$ID,collapse = ",")
      
      trns=trns[grep(",",d2$Same_ensembl_gene_id)]
      d3$Same_ensembl_gene_id = paste0(unique(c(d2$Same_ensembl_gene_id,unlist(strsplit(trns,",")))),
                                       collapse = ",")
      trns=trns[grep(",",d2$Same_gene_name_comb)]
      d3$Same_gene_name_comb=paste0(unique(c(d2$Same_gene_name_comb,unlist(strsplit(trns,",")))),
                                    collapse = ",")
      trns=trns[grep(",",d2$`Opposite Strand: Host_gene_ensembl_id`)]
      d3$Oppo_ensembl_gene_id=paste0(unique(c(d2$Oppo_ensembl_gene_id,unlist(strsplit(trns,",")))),
                                     collapse = ",")
      trns=trns[grep(",",d2$Oppo_gene_name_comb)]
      d3$Oppo_gene_name_comb=paste0(unique(c(d2$Oppo_gene_name_comb,unlist(strsplit(trns,",")))),
                                    collapse = ",")
      CollapsedOut=rbind(CollapsedOut,d3)
    }
    
    #remove rownames, col
    rownames(CollapsedOut)=NULL
    CollapsedOut=CollapsedOut[-1,]
    CollapsedOutID= unique(unlist(strsplit(CollapsedOut$IDmerge,",")))
    
    ###fix we dont want to overwrite original df
    Peaksdata2_anno$IDmerge=NA
    Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$ID%in%CollapsedOutID==F,]
    Peaksdata2_anno=(rbind(Peaksdata2_anno,CollapsedOut))
  }
}

#write out for junction annotation 
write.table(Peaksdata2_anno,paste0(out_file),sep = "\t")

#write out for mapq
write.table(Peaksdata2_anno[,c('chr','start','end','strand','ID')],
            file=paste0(out_dir,file_id,"peakannotation_mapq_IN.txt"), 
            sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")

