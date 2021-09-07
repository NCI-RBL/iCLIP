#library
library(tidyr)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(data.table)
library(argparse)

#set args
parser <- ArgumentParser()
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

##testing
testing="Y"
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
  
  #feature information
  join_junction = "TRUE"
  condense_exon="TRUE"
  read_depth = 3
  DEmethod = "None"
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
FtrCount_input<-function(count.file){
  FTR_df=read.delim(count.file, header=T, sep="\t",
                    stringsAsFactors = F,comment.char = '#')
  colnames(FTR_df)[7]='Counts'
  
  #remove rows without start,end
  FTR_df=FTR_df[!is.na(FTR_df$Start) | !is.na(FTR_df$End),]
  
  #create unique ids
  FTR_df$ID=paste0(FTR_df$Chr,":",FTR_df$Start,"-",FTR_df$End,"_",FTR_df$Strand)
  
  return(FTR_df)
}
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
      
      #merge JCount with unique/all merge
      FtrCount_fracJCount=FtrCount_fracJCount %>%
        merge(
          FtrCount[,c("ID","strand")],
          by.x="PrimaryGene",
          by.y="ID",
          all.x=T) %>%
        separate(
          PrimaryGene,
          into=c('ID','strand'),
          sep = "_",
          remove = F) %>%
        separate(
          ID,
          into=c('chr','start','end'),
          sep = ":|-",
          remove = T)
      print(paste0('check rows: merge jcounts - ',nrow(FtrCount_fracJCount)))
      
      #remove junctions where splicing is within peak
      site_1=FtrCount_fracJCount$Site1_location>=FtrCount_fracJCount$start&
        FtrCount_fracJCount$Site1_location<=FtrCount_fracJCount$end
      site_2=FtrCount_fracJCount$Site2_location>=FtrCount_fracJCount$start&
        FtrCount_fracJCount$Site2_location<=FtrCount_fracJCount$end
      FtrCount_fracJCount=FtrCount_fracJCount[site_1!=site_2,]
      FtrCount_fracJCount$JunctionID=paste0(FtrCount_fracJCount$Site1_chr,':',
                                            FtrCount_fracJCount$Site1_location,'-',
                                            FtrCount_fracJCount$Site2_location,'_',
                                            FtrCount_fracJCount$strand)
      #renumber rownames
      FtrCount_fracJCount$rID = 1:nrow(FtrCount_fracJCount)
      
      #create GR objects
      jcount1.GR=GRanges(seqnames = as.character(FtrCount_fracJCount$Site1_chr), 
                         ranges=IRanges(start = as.numeric(FtrCount_fracJCount$Site1_location), 
                                        end = as.numeric(FtrCount_fracJCount$Site1_location)),
                         strand = FtrCount_fracJCount$strand,
                         PrimaryGene=FtrCount_fracJCount$PrimaryGene,
                         count=FtrCount_fracJCount$counts,
                         JunctionID=FtrCount_fracJCount$JunctionID,
                         rID=FtrCount_fracJCount$rID)
      jcount2.GR=GRanges(seqnames = as.character(FtrCount_fracJCount$Site2_chr), 
                         ranges=IRanges(start = as.numeric(FtrCount_fracJCount$Site2_location), 
                                        end = as.numeric(FtrCount_fracJCount$Site2_location)),
                         strand = FtrCount_fracJCount$strand,
                         PrimaryGene=FtrCount_fracJCount$PrimaryGene,
                         count=FtrCount_fracJCount$counts,
                         JunctionID=FtrCount_fracJCount$JunctionID,
                         rID=FtrCount_fracJCount$rID)
      FtrCount.GR=GRanges(seqnames = as.character(FtrCount$chr), 
                          ranges=IRanges(start = as.numeric(FtrCount$start), 
                                         end = as.numeric(FtrCount$end)),
                          strand = FtrCount$strand,
                          ID=FtrCount$ID)
      
      ############### Add Junction ID
      junctionID <- function(Ftr_in, j_in, j_name, p_name){
        xo=as.data.frame(GenomicRanges::findOverlaps(j_in, 
                                                     Ftr_in, 
                                                     type = "any",
                                                     ignore.strand=F))
        qh=as.data.frame(j_in[xo$queryHits],row.names = NULL)
        colnames(qh)=paste0(colnames(qh),j_name)
        
        sh=as.data.frame(Ftr_in[xo$subjectHits],row.names = NULL)
        colnames(sh)=paste0(colnames(sh),p_name)
        
        return(cbind(qh,sh))
      }
      
      ### Junction 1 ID
      Junc_PL1 = junctionID(FtrCount.GR, jcount1.GR,"_Junction1","_Peaks1")
      ### Junction 2 ID
      Junc_PL2 = junctionID(FtrCount.GR, jcount2.GR,"_Junction2","_Peaks2")
      
      #merge and filter
      Junc_PL=merge(Junc_PL1,Junc_PL2,by.x='rID_Junction1',by.y='rID_Junction2')[,c('PrimaryGene_Junction1',
                                                                                    'JunctionID_Junction1',
                                                                                    'strand_Peaks1',
                                                                                    'ID_Peaks1',
                                                                                    'ID_Peaks2')]
      Junc_PL= Junc_PL %>%
        base::subset(!is.na(ID_Peaks1) | !is.na(ID_Peaks2)) %>%
        dplyr::distinct(ID_Peaks1,ID_Peaks2, .keep_all= TRUE)
      
      #create join id col
      Junc_PL$JoinID=paste0(Junc_PL$ID_Peaks1,",",Junc_PL$ID_Peaks2)
      
      ############### Combine all connected peaks
      #create output df
      Junc_PL$matchid=paste(Junc_PL$ID_Peaks1,Junc_PL$ID_Peaks2,sep="|")
      Junc_PL = Junc_PL %>%
        separate(ID_Peaks1,
                 sep = "_",
                 into = c('ID','strand'),remove = F) %>%
        separate(ID, 
                 sep = ":|-",
                 into = c('chr','start','end'),
                 remove = T)
      
      Junc_PL[c("start","end")] <- sapply(Junc_PL[c("start","end")],as.numeric)
      
      PGene_TBL2 = setNames(data.frame(matrix(ncol = 9, nrow = length(unique(Junc_PL$JoinID))),row.names =unique(Junc_PL$matchid)), 
                            c("PrimaryGene_Junction1", "JunctionID_Junction1", "JoinID","chr",
                              "start","end","strand","ID_Peaks1","ID_Peaks2"))   
      
      Junc_PL_short=Junc_PL
      for (rowid in unique(Junc_PL$matchid)){
        #generate search id
        #IE "chr1:100007068-100007156_\\+|chr1:100082000-100082149_\\+"
        #determine which rows in df partially match any part of search id
        
        match_df=Junc_PL_short[grepl(gsub("_","_\\\\",rowid),
                                     Junc_PL_short$JoinID),]
        if (nrow(match_df)==0) {next}
        
        ## search only peaks in same chrm of interest
        Junc_PL_chrslct=Junc_PL_short[Junc_PL_short$chr%in%unique(match_df$chr),]
        
        #until we find all rows matching all partial matches
        match_original=str_split(rowid,"\\|")%>%unlist()%>%unique()
        match_current=str_split(match_df$matchid,"\\|")%>%unlist()%>%unique()
        
        
        while(length(setdiff(match_current,match_original))>0){ 
          for (rowid2 in setdiff(match_current,match_original)){
            #determine which rows in df partially match any part of search id
            rowid2=paste(rowid2,collapse = '|')
            match_df2=Junc_PL_chrslct[grepl(gsub("_","_\\\\",rowid2),
                                            Junc_PL_chrslct$matchid),]
            
            #original list is current, and current list is new matched df
            match_original=match_current
            match_add=str_split(match_df2$matchid,"\\|")%>%unlist()%>%unique()
            
            match_current=unique(c(match_add,match_current))
            match_df=rbind(match_df,match_df2) %>%
              distinct
          }
        }
        #update df
        PGene_TBL2[rowid,] = list(paste(sort(unique(match_df$PrimaryGene_Junction1)),collapse = ','),
                                  paste(sort(unique(match_df$JunctionID_Junction1)),collapse = ','),
                                  paste(sort(unique(match_df$JoinID)),collapse = ','),
                                  unique(match_df$chr),
                                  min((match_df[,c('start','end')])),
                                  max((match_df[,c('start','end')])),
                                  unique(as.character(match_df$strand_Peaks1)),
                                  paste(sort(unique(match_df$ID_Peaks1)),collapse = ','),
                                  paste(sort(unique(match_df$ID_Peaks2)),collapse = ','))
        
        #### Remove Joined ID from master junction table Junc_PLOut
        # Junc_PLOut=Junc_PLOut[idcomb$ID_Peaks2%in%Junc_PLOut==F,]
        Junc_PL_short=Junc_PL_short[grepl(gsub("_","_\\\\",paste(unique(c(match_df$ID_Peaks1,match_df$ID_Peaks2)),collapse = "|")),Junc_PL_short$JoinID)==F,]
      }
      
      filter
      PGene_TBL2=PGene_TBL2 %>%
        subset(!is.na(PrimaryGene_Junction1)) %>%
        distinct(.keep_all=TRUE)
      
      #create linked id
      PGene_TBL2$linkedID=paste0(PGene_TBL2$chr,":",
                                 PGene_TBL2$start,"-",
                                 PGene_TBL2$end,"_",
                                 PGene_TBL2$strand)
      
      PGene_TBL2=PGene_TBL2[duplicated(PGene_TBL2)==F,]
      
      
      ############### Get all peaks with junctions
      Junc_peaks=unique(c(Junc_PL$ID_Peaks1,Junc_PL$ID_Peaks2))
      
      #Trim peaks without Splicing
      # FtrCount=FtrCount[FtrCount$ID%in%Junc_peaks==F & FtrCount$Counts_fracMM>=read_depth,]
      FtrCount_trimmed=rbind(FtrCount[FtrCount$ID%in%Junc_peaks==F & FtrCount$Counts_fracMM>=read_depth,],## peaks not spliced and count > readDepth cut
                             FtrCount[FtrCount$ID%in%Junc_peaks,]) ## all peaks with spliced reads 
      
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
if (DEmethod!='none') {
  print(paste0("Running ",DEmethod," Peaks file prep"))
  
  manorm_bed=FtrCount_trimmed[,c('chr','start','end','ID','ID','strand')]
  manorm_bed[,5]=1
  
  #write peaks bed file
  write.table(manorm_bed,
              file=paste0(out_dir_DEP, sample_id, "_", "Peaksfor",DEmethod,".bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
  #write peaks + bed file
  write.table(manorm_bed[manorm_bed$strand%in%"+",],
              file=paste0(out_dir_DEP, sample_id, "_", "Peaksfor",DEmethod,"_P.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
  
  #write peaks - bed file
  write.table(manorm_bed[manorm_bed$strand%in%"-",],
              file=paste0(out_dir_DEP, sample_id, "_", "Peaksfor",DEmethod,"_N.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
} else{
  print("Differential Peaks skipped")
} 

##########################################################################################
############### reference, annotation table
##########################################################################################
print("Prep References")

#read in references
ref_table = read.table(reftable_path, header = TRUE, sep="\t", row.names = 1)

gencodeCol=c('chr','start','end','strand','ensembl_gene_id',
             'transcript_id','external_gene_name','gene_type','gene_type_ALL')
ref_gencode_trans = read.table(gencode_transc_path,header=TRUE,sep="\t")
gencode_Anno_RNA_comb=ref_gencode_trans[,gencodeCol]  

#create annotation dfs
rpmk_Anno_RNA_comb=data.frame()
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
        Gtable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=T,sep="\t")
        gencode_Anno_RNA_comb=rbind(gencode_Anno_RNA_comb,Gtable[,gencodeCol])
      } else if (selec =='Repeatmasker') {
        Rtable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=TRUE,sep="\t")
        rpmk_Anno_RNA_comb=rbind(rpmk_Anno_RNA_comb,Rtable)
      } else if (selec =='Custom') {
        Ctable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=TRUE,sep="\t")
        custom_Anno_RNA_comb=rbind(custom_Anno_RNA_comb,Ctable)
      } else if (selec =='refseq') {
        Reftable=read.table(paste0(anno_dir,rowid,'_',selec,'.bed'),header=TRUE,sep="\t")
        refseq_Anno_RNA_comb=rbind(refseq_Anno_RNA_comb,Reftable)
      }
      
    }
  }else { next }
}  

gencode_Anno_RNA_comb=gencode_Anno_RNA_comb[duplicated(gencode_Anno_RNA_comb)==F,]

#prepare annotations
prep_annotation<-function(df.in,col.in){
  
  df_clean = df.in %>%
    distinct()
  
  if (dim(df_clean)[1] == 0) {
    df_out=as.data.frame(matrix(ncol=length(c("chr","start","end","strand","type","transcript_name")),nrow=0))
    colnames(df_out)=c("chr","start","end","strand","type","transcript_name")
  } else {
    df_out=df_clean[,col.in]
    colnames(df_out)=c("chr","start","end","strand","type","transcript_name")
  }
  return(df_out)
}

Anno_RNA_comb = rbind(
  prep_annotation(rpmk_Anno_RNA_comb,c("chr","start","end","strand","type","name")),
  prep_annotation(custom_Anno_RNA_comb,c("chr","start","end","strand","type","name")),
  prep_annotation(refseq_Anno_RNA_comb,c("chr","start","end","strand","gene_type","ensembl_gene_id"))) %>%
  distinct()
Anno_RNA_comb=Anno_RNA_comb[is.na(Anno_RNA_comb$chr)==F,]
################################################################################## SSC cleaned above

##########################################################################################
############### call peaks
##########################################################################################
print("Call peaks")

#merge references depending on species, run annotation
bam_anno2<-function(peaksTable,Annotable,ColumnName,pass_n){
  
  #set file prefix
  file_prefix = paste0(file_id,pass_n,"_")
  
  ## the multiple chr is a space filler was to make a bed6 file for bedtools 
  anno_output = Annotable[,c('chr','start','end','chr','chr','strand',ColumnName)] %>%
    dplyr::rename(
      chr_anno = chr,
      start_anno = start,
      end_anno = end,
      strand_anno = strand
    )
  
  #clear chr cols 4,5
  anno_output[,4:5]="."
  
  #output bed files
  write.table(anno_output,
              file = paste0(tmp_dir,file_prefix,"annotable.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  #write out peakstable
  peaksTable_output=subset(peaksTable,!is.na(peaksTable$start) & 
                             chr %in% unique(anno_output$chr_anno))[,c('chr','start','end','ID','ID','strand')]
  write.table(peaksTable_output,
              file = paste0(tmp_dir,file_prefix,"peakstable.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  # merge bedfiles into output text
  system(paste0('bedtools intersect -a ', 
                tmp_dir, file_prefix,'peakstable.bed -b ', 
                tmp_dir, file_prefix,'annotable.bed -wao -s > ',
                tmp_dir, file_prefix,'peaks_OL.txt'))
  
  #read in merged file
  ab_OL = read.table(paste0(tmp_dir,file_prefix,"peaks_OL.txt"), 
                     header=F, sep="\t",stringsAsFactors = F)
  colnames(ab_OL) = c(paste0(colnames(peaksTable_output)),
                      paste0(colnames(anno_output)),'ntOL')
  
  #add width
  ab_OL$width = ab_OL$end-ab_OL$start
  ab_OL$width_anno = ab_OL$end_anno-ab_OL$start_anno
  
  #subset 
  ab_OL=ab_OL[ab_OL$ntOL>0,]
  if (nrow(ab_OL)>0){
    #add nt / width in annotation and read
    ab_OL$OLper_anno = ab_OL$ntOL/ab_OL$width_anno
    ab_OL$OLper = ab_OL$ntOL/ab_OL$width
    
    #subset
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
    print(paste0('Overlaps with ',nrow(ab_OL_colapsed),' Regions from '))[1]
    
    #merge tables
    ab_OL_colapsed = merge(peaksTable,
                           ab_OL_colapsed[,c('ID',ColumnName)],
                           by='ID',all.x=T)
    
  } else{ 
    print("does not overlap with ")[1]
    # next
    ab_OL_colapsed=peaksTable
    ab_OL_colapsed[,ColumnName]=NA
  }
  
  print((match.call())$Annotable)
  
  return(ab_OL_colapsed) 
}

peak_calling<-function(peak_in,nmeprfix){
  # peak_in=peaks_Oppo
  # nmeprfix="Oppo_"
  
  colMerge = c('chr','start','end','strand','ensembl_gene_id','transcript_id',
               'external_gene_name','gene_type','gene_type_ALL')
  colSelect=c('ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')
  
  PeaksdataOut=bam_anno2(peak_in,
                         gencode_Anno_RNA_comb[,colMerge],
                         colSelect,
                         "pass1")
  
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
  ############### Summarize Peak annotation
  ##########################################################################################
  #separate annotated and non annotated peaks 
  anno=PeaksdataOut[(is.na(PeaksdataOut$gene_type_ALL)==F)|(is.na(PeaksdataOut$type)==F),'ID']
  
  #Only subset of peaks are annotated (most likely)
  if (length(anno)<nrow(PeaksdataOut)&length(anno)>0) { 
    peaksTable_Anno=PeaksdataOut[(PeaksdataOut$ID%in%anno),,drop=F]
    peaksTable_noAnno=PeaksdataOut[!(PeaksdataOut$ID%in%anno),,drop=F]
  } else if (length(anno)==nrow(PeaksdataOut)) { #"All peaks are annotated"
    peaksTable_Anno=PeaksdataOut[(PeaksdataOut$ID%in%anno),,drop=F]
    peaksTable_noAnno=as.data.frame(matrix(nrow=1,ncol=ncol(PeaksdataOut)))
    colnames(peaksTable_noAnno)=colnames(PeaksdataOut)
  }
  
  if(length(anno)>0) { 
    
    u=unique(peaksTable_Anno$ID)
    peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_Anno)));
    colnames(peaksTable_colapsed)=colnames(peaksTable_Anno)
    
    ## Checking for conflicting annotations from gencode and additional supplied annotations
    ## If conflict prioritize Supplied annotations 
    ## Additionally Peak may overlap multiple features and so I am setting hierarchy
    for(x in 1:length(u)){
      PeaksdataOutU=u[x]
      pam=peaksTable_Anno[peaksTable_Anno$ID%in%PeaksdataOutU,]
      
      #Comb Annotation
      pam_1=as.character(pam$gene_type_ALL) ### type from Gencode 
      pam_1=as.data.frame(strsplit(pam_1,','));colnames(pam_1)='a';pam_1$a=as.character(pam_1$a)
      
      pam_2=as.character(pam$type) ### type from Annotation
      pam_2=as.data.frame(strsplit(pam_2,','));colnames(pam_2)='a';pam_2$a=as.character(pam_2$a)
      
      type_list = c("miRNA","piRNA","rRNA","siRNA","snRNA","snoRNA",
                    "ribozyme","scRNA", "sRNA","scaRNA", "yRNA", "srpRNA", "7SKRNA",
                    "sncRNA","vaultRNA", "tRNA")
      #remove misc_RNA catagory : these are covered by additional annotations (mostly yRNA)
      ## We should set this up to only display misc_RNA if no additional annotations
      if (nrow(pam_1)>1&grepl('misc_RNA',pam_1)&grepl(paste(type_list,collapse = "|"),pam_1)){
        pam_1[pam_1$a%in%'misc_RNA','a']=NA
      }
      if (grepl('misc_RNA',pam_1)&'FALSE'%in%is.na(pam_2)){
        pam_1[pam_1$a%in%'misc_RNA','a']=NA
      }
      
      #lincRNA are from an older Gencode version (VM18 or older) so don't use    
      if (grepl('lncRNA',pam_1)&grepl('lincRNA',pam_2)) {
        pam_1[pam_2$a%in%'lincRNA',]='lncRNA'
      }
      
      #change ribozyme (gencode) to RNA type from additional anno
      if (grepl('ribozyme',pam_1)) {
        pam_1[pam_1$a%in%'ribozyme',]=unique(pam_2$a)
      }
      
      #combine all rna types subtypes into ncRNA
      pam_c=rbind(pam_1,pam_2)
      pam_c=pam_c[!is.na(pam_c$a),,drop=F]
      
      #### Annotation from Gencode : unique(sort(mm10$gene_type_ALL))
      ReplaceType <-function(df_in,col_in,list_in,replace_id){
        for (id in list_in){
          df_in[,col_in]=gsub(paste0("\\<",id,"\\>"),replace_id,df_in[,col_in])
        } 
        return(df_in)
      }
      
      type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA",
                    "ribozyme","scRNA", "sRNA","scaRNA", "yRNA", "srpRNA", "7SKRNA",
                    "sncRNA",'vaultRNA', "tRNA")
      pam_c2 = ReplaceType(pam_c,"a", type_list,"ncRNA")
      
      if (length(pam_c$a)>1) {pam_c=pam_c[order(pam_c$a),,drop=F]}
      if (length(pam_c2$a)>1) {pam_c2=pam_c2[order(pam_c2$a),,drop=F]}
      
      pam_c=paste(unique(pam_c$a),collapse = ',')
      pam_c2=paste(unique(pam_c2$a),collapse = ',')
      
      #Comb GeneName 
      pam_G1=as.character(pam$external_gene_name) ### gene name from gencode
      pam_G1=as.data.frame(strsplit(pam_G1,','));colnames(pam_G1)='a';pam_G1$a=as.character(pam_G1$a)
      
      ### G1 is gencode annotation and G2 is additional annotation
      pam_G2=as.character(pam$transcript_name) ### gene name from annotation
      pam_G2=as.data.frame(strsplit(pam_G2,','));colnames(pam_G2)='a';pam_G2$a=as.character(pam_G2$a)
      
      ## there are non uniform naming for 7SK rRNA so i am making a uniform naming by selecting only from gencode
      if (grepl('sk',pam_G1)&grepl('7SK',pam_G2)) {
        pam_G2[pam_G2$a%in%'7SK',]=NA
      }
      
      #merge gene names from "gencode" and "annotation"
      pam_gmerge=rbind(pam_G1,pam_G2)
      pam_gmerge=pam_gmerge[!is.na(pam_gmerge$a),,drop=F]
      pam_gmerge=paste(unique(pam_gmerge$a),collapse = ',')
      
      peaksTable_colapsed[x,colnames(peaksTable_Anno)]=(pam[1,colnames(peaksTable_Anno),drop=T])
      peaksTable_colapsed[x,'type_simple_comb']=pam_c2
      peaksTable_colapsed[x,'type_comb']=pam_c
      peaksTable_colapsed[x,'gene_name_comb']=pam_gmerge
    }
    
    
    rnatype=peaksTable_colapsed[,c('ID','ensembl_gene_id','external_gene_name','gene_type',"gene_type_ALL",
                                   'type','transcript_name','type_simple_comb','type_comb','gene_name_comb')]
    
    peaksTable_noAnno[,c('type_simple_comb','type_comb','gene_name_comb')]=NA
    
    peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_noAnno)
    peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
    
    rnames=c("ensembl_gene_id","external_gene_name","gene_type","gene_type_ALL","type","transcript_name","type_simple_comb","type_comb","gene_name_comb")
    
    colnames(peaksTable_colapsed)[colnames(peaksTable_colapsed)%in%rnames]=paste0(nmeprfix,rnames)
    
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
    
  } else { 
    print("No Peaks Annotated")
    peaksTable_colapsed=PeaksdataOut
    rnames=c("ensembl_gene_id","external_gene_name","gene_type","gene_type_ALL","type","transcript_name","type_simple_comb","type_comb","gene_name_comb")
    peaksTable_colapsed[,rnames]=NA
    colnames(peaksTable_colapsed)[colnames(peaksTable_colapsed)%in%rnames]=paste0(nmeprfix,rnames)
  }  
  
  return(peaksTable_colapsed)
}


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

#gencode
ref_gencode = fread(gencode_path, header=T, sep="\t",stringsAsFactors = F,data.table=F) %>%
  subset(transcript_type != "TEC" & gene_type != "TEC") %>% 
  dplyr::rename(
    chr = seqname,
    ensembl_gene_id = gene_id,
    external_gene_name = gene_name
  )

#exons
ref_gencode_e = base::subset( ref_gencode, feature == "exon")

#introns
introns=fread(intron_path, 
              header=F, sep="\t",stringsAsFactors = F,data.table=F, 
              col.names = c('chr','start','end','attribute','V5','strand')) %>%
  separate(attribute,
           into=c('transcript_id','feature','exon_number','level','chr2','intronnumber','dir'),
           remove = T,sep = "_")

#remove Non-Chromosome contigs
introns$chr=(gsub('chr[0-9]+_|chr[X-Y]_|chrUn_|_alt|_random|_fix','',introns$chr))
introns$chr=(gsub('v','.',introns$chr))

#split out transcript id
introns$transcript_id = unlist(lapply(stringr::str_split(introns$transcript_id, "[.]"), "[[",1))

# intron table is 0 based while exon starts at 1
# Same for exon number the exon # starts at 0
introns$start=introns$start+1
introns$exon_number=as.numeric(introns$exon_number)+1

# add gene type type
introns=merge(introns,
              unique(ref_gencode_trans[,c('transcript_id',
                                          'ensembl_gene_id','external_gene_name',
                                          'gene_type','gene_type_ALL','transcript_type',
                                          'transcript_name','score')]),
              by='transcript_id',all.x=T)

#binding exon, intron
intron_exon=rbind(ref_gencode_e[,c('chr','feature','start','end','strand','transcript_id',
                                   'ensembl_gene_id','external_gene_name','gene_type',
                                   'transcript_type','transcript_name','score','exon_number')],
                  introns[,c('chr','feature','start','end','strand','transcript_id',
                             'ensembl_gene_id','external_gene_name','gene_type','transcript_type',
                             'transcript_name','score','exon_number')])
intron_exon$ID=paste0(intron_exon$chr,':',intron_exon$start,'-',intron_exon$end,"_",intron_exon$strand)

##########################################################################################
############### Intron Exon Calling
##########################################################################################
print("IE Calling")
IE_calling <- function(peak_in,Peak_Strand_ref,nmeprfix){
  #create annotation table
  ColumnName = c("feature","exon_number")
  anno_IntExn = intron_exon[grep('protein_coding',intron_exon$gene_type),][,c('chr','start',
                                                                              'end','transcript_id',
                                                                              'transcript_id','strand',ColumnName)]
  
  ###fix variable swap
  p=peak_in[,c('chr','start','end','ID','ID','strand')]
  
  # Setting up bed 6 file but need column to be clear for bedtools intersect
  anno_IntExn[,5]=""
  
  #add _anno to col names
  colnames(anno_IntExn)[colnames(anno_IntExn)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')
  
  ## setting up bedtools intersect
  write.table(anno_IntExn,file=paste0(out_dir,file_id,"annotable.bed"), sep = "\t", 
              row.names = FALSE, col.names = F, append = F, quote= FALSE)
  write.table(p,file=paste0(tmp_dir,file_id,"peakstable.bed"), sep = "\t", 
              row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  system(paste0('bedtools intersect -a ', tmp_dir, file_id, 'peakstable.bed -b ', 
                out_dir, file_id, 'annotable.bed -wao -s  >', tmp_dir, file_id, 'peaks_OL.txt'))
  
  exoninof=fread(paste0(tmp_dir,file_id,"peaks_OL.txt"), 
                 header=F, sep="\t",stringsAsFactors = F,data.table=F)
  
  #set up df
  colnames(exoninof)=c(paste0(colnames(p)),paste0(colnames(anno_IntExn)),'ntOL')
  exoninof=exoninof[exoninof$ntOL>0,]
  exoninof$width_anno=exoninof$end_anno-exoninof$start_anno
  exoninof$OLper_anno=exoninof$ntOL/exoninof$width_anno
  exoninof$width=exoninof$end-exoninof$start
  exoninof$OLper=exoninof$ntOL/exoninof$width
  exoninof=exoninof[(exoninof$OLper_anno>.75 | exoninof$OLper>.51), ]
  
  #run same strand
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
    peak_in[,paste0(nmeprfix,c('Exn_start_dist','Intron_start_dist','Intron_5pStart','Exn_5pStart'))]=NA
  }
  
  #for each peak in table
  for (x in 1:nrow(peak_in)) {
    
    #subset ids
    l=peak_in[x,c('ID','ID')]
    g=exoninof[exoninof$ID%in%l$ID,]
    
    #if there are ids in both lists
    if (nrow(g)>0) {
      
      #split into exon and intron
      g_e=g[g$feature%in%'exon',]
      g_i=g[g$feature%in%'intron',]
      
      gname=paste(unique(g$feature),collapse = ",")
      gname=gsub('NA,',"",gname)
      gname=gsub(',NA',"",gname)
      peak_in[x,paste0(nmeprfix,'feature')]=gsub('NA,',"",gname)
      
      #if there are exons
      if (nrow(g_e)>0) {
        gname=paste(unique(g_e$exon_number),collapse = ",")
        gname=gsub('NA,',"",gname)
        gname=gsub(',NA',"",gname)
        peak_in[x,paste0(nmeprfix,'exon_number')]=gsub('NA,',"",gname)
        peak_in[x,paste0(nmeprfix,'exon_LargeOL')]=max(g_e$ntOL/(g_e$end-g_e$start))
        
        #if the Same
        if (Peak_Strand_ref=='Same') {
          g_e[g_e[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
          peak_in[x,paste0(nmeprfix,'Exn_start_dist')]=mean(as.numeric(g_e[,paste0(nmeprfix,'feature_Distance')]))
        } 
      } 
      
      #if there are introns
      if (nrow(g_i)>0) {
        gname=paste(unique(g_i$exon_number),collapse = ",")
        gname=gsub('NA,',"",gname)
        gname=gsub(',NA',"",gname)
        peak_in[x,paste0(nmeprfix,'intron_number')]=gsub('NA,',"",gname)
        peak_in[x,paste0(nmeprfix,'intron_LargeOL')]=max(g_i$ntOL/(g_i$end-g_i$start))
        
        if (Peak_Strand_ref=='Same') {
          g_i[g_i[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
          peak_in[x,paste0(nmeprfix,'Intron_start_dist')]=mean(as.numeric(g_i[,paste0(nmeprfix,'feature_Distance')]))
        }
      } 
      
    } else {
      if (Peak_Strand_ref=='Same') {
        peak_in[x,paste0(nmeprfix,c("feature",'exon_number','intron_number','Exn_start_dist','Intron_start_dist'))]=NA
      } else if (Peak_Strand_ref=='Oppo') {
        peak_in[x,paste0(nmeprfix,c("feature",'exon_number','intron_number'))]=NA
      }
    }
  }
  
  ### protein coding peaks with not Correctly assigned exon overlap
  peak_in[is.na(peak_in[,paste0(nmeprfix,'feature')]) & 
            peak_in[,paste0(nmeprfix,'gene_type')]%in%'protein_coding',paste0(nmeprfix,'feature')]='exon'
  return(peak_in)
}

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
rmsk_GRCm38=fread(rmsk_path, header=T, sep="\t",stringsAsFactors = F,
                  data.table=F)
rmsk_GRCm38$genoName=gsub('chr[0-9]+_|chr[X-Y]_|chrUn_|_alt|_random|_fix',"",rmsk_GRCm38$genoName)
rmsk_GRCm38$genoName=gsub("v",".",rmsk_GRCm38$genoName)

rpmsk_anno <- function(ColumnName,Annotable,peaksTable){
  
  if(length(args)==0){
    nmeprfix="Same_"
    peaksTable=intronexon_Same
    Annotable=rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('LINE','SINE'),]
    ColumnName=paste0(nmeprfix,'Repeat_LINE_SINE')
  }
  
  #create df of ids and start location
  peaksTable=separate(peaksTable,ID,into=c('start','strand'),sep="_",remove=F)
  peaksTable=separate(peaksTable,start,into=c('chr','start','end'),sep=":|-",remove=F)
  
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
    peaksTable = u[x]
    pam = peaksTable_double[peaksTable_double$ID%in%peaksTable,]
    
    #if the cols of double are not in ColumnName list
    peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F] = (pam[1,colnames(peaksTable_double)%in%c(ColumnName)==F,drop=T])
    
    #past the uqique values of the col to collapsed table
    peaksTable_colapsed[x,ColumnName]=paste(sort(unique((pam[,ColumnName]))), collapse =" | ")
  }
  
  #bind collapsed and unqiue
  peaksTable_colapsed = rbind(peaksTable_colapsed,peaksTable_single)
  
  #if there is no ID remove
  peaksTable_colapsed = peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  ###phil why are we adding NA columns ## This is a check to make sure rows with empty annotations get an NA
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
    rmsk_subset = base::subset(rmsk_GRCm38, repClass == find_type)
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
############### Assigning Clip peak attributes by strand
##########################################################################################
# Not all Peaks overlap with a single feature so peak assignments were assigned by priority:
# ncRNA > Protein coding : Exonic > repeats > Pseudogene > Antisense Feature > 
#Protein Coding : Intronic > lncRNA > no Feature  
# All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table
print("Add peak attributes")
peak_attributes <- function(peaks_in,nmeprfix){
  if(length(args)==0){
    nmeprfix="Same_"
    peaks_in=rpmsk_Same
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
  
  #create col
  peaks_in[,paste0(nmeprfix,'Comb_type_ncRNA')]=NA
  
  #select ncRNA
  p_ncrna=  p_ncrna=peaks_in[peaks_in[,paste0(nmeprfix,'Comb_type_exon')]=='ncRNA',]
  
  #make col type_comb with prefix  
  p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=p_ncrna[,paste0(nmeprfix,'type_comb')]
  
  #remove types
  ### Protein coding + ncRNA annotations -> ncRNA subtype
  ### any ncRNA takes priority over pseudogene
  ### any double annotations with lncRNA become second annotation only
  type_list = c("protein_coding,",",protein_coding","pseudogene,",",pseudogene",'lnLcRNA,')
  for (type_id in type_list){
    p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(type_id,"",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  }  
  
  #merge linLcRNA  
  p_ncrna[p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-intron,linLcRNA-exon',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
  p_ncrna[p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-exon,linLcRNA-intron',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
  
  type_list = c("lnLcRNA,",'linLcRNA,','linLcRNA-exon,',"linLcRNA-intron,")
  for (type_id in type_list){
    p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(type_id,"",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  } 
  
  
  #### This section will give priority to ncRNA assignment based on order listed in ncRNAPriority
  #### This replaces previous approach to replace on a case by case baseis
  ncRNAPriority=c(
  'yRNA',
  'snRNA',
  'snoRNA',
  'srpRNA',
  'tRNA',
  '7SKRNA',
  'scRNA',
  'sncRNA',
  'scaRNA',
  'miRNA',
  'vaultRNA',
  'miscRNA',
  'misc_RNA',
  'rRNA',
  'linLcRNA-intron',
  'linLcRNA-exon',
  'lnLcRNA') 
  
  for (x in 1:length(ncRNAPriority)) {
    nc=paste0(",",ncRNAPriority[x],"|",ncRNAPriority[x],",")
    p_ncrna[grep(nc,p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]),paste0(nmeprfix,'Comb_type_ncRNA')]=ncRNAPriority[x]
  }
  
  p_ncrna[grep('miscRNA',p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]),paste0(nmeprfix,'Comb_type_ncRNA')]='misc_RNA'
  
  
  # ### any double annotations with yRNA become yRNA only
  # p_ncrna[grep('yRNA',p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]),paste0(nmeprfix,'Comb_type_ncRNA')]='yRNA'
  # 
  # ### tRNA takes priority over miRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,tRNA',"miRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # 
  # ### miRNA takes priority over rRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,rRNA',"miRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # 
  # ### snRNA takes priority over 7SKRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,snRNA',"snRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('snRNA,7SKRNA',"snRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # 
  # ### 7SKRNA takes priority over lncRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,linLcRNA-intron',"7SKRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(',linLcRNA-intron,7SKRNA',"7SKRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,linLcRNA-exon',"7SKRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA-exon,7SKRNA',"7SKRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,lnLcRNA',"7SKRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('lnLcRNA,7SKRNA',"7SKRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # 
  # ### tRNA takes priority over rRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('rRNA,tRNA',"tRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('tRNA,rRNA,',"tRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # 
  # ### snoRNA takes priority over miRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,snoRNA',"snoRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('snoRNA,miRNA,',"snoRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # 
  # ### scaRNA takes priority over miRNA
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,scaRNA',"scaRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  # p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('scaRNA,miRNA,',"scaRNA",p_ncrna[,paste0(nmeprfix,'Comb_type_ncRNA')])
  
  #select not ncRNA
  p2=peaks_in[!peaks_in[,paste0(nmeprfix,'Comb_type_exon')]%in%'ncRNA',]
  
  ###fix variable
  PeaksdataOut=rbind(p_ncrna,p2)
  
  return(PeaksdataOut)
}

#merge attributes
peak_attrib_Oppo = peak_attributes(rpmsk_Opposite,"Oppo_")
PeaksdataOut = merge(peak_attributes(rpmsk_Same,"Same_"),
                     peak_attrib_Oppo[,colnames(peak_attrib_Oppo)[!colnames(peak_attrib_Oppo) %in% 
                                                                    c("chr","start","end","strand" )]],
                     by='ID')
##########################################################################################
############### Assigning Clip peak attributes - merged strands
##########################################################################################  
#create column
PeaksdataOut$Comb_type_exon_Oppo="no Feature"

#select subset of df, input specific type
AddTypetoComb <- function(df_in,subset_df,type_id){
  df_in[subset_df,paste0('Comb_type_exon_Oppo')]=type_id
  return(df_in)
}

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

PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
PeaksdataOut[,'Comb_type_exon_Oppo']=factor(PeaksdataOut[,'Comb_type_exon_Oppo'], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))

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
  Peaksdata2_anno_trns_exon=( Peaksdata2_anno[Peaksdata2_anno$ID%in%Junc_peaks,])
  CollapsedOut=as.data.frame(matrix(nrow = (1),ncol = ncol(Peaksdata2_anno_trns_exon)))
  colnames(CollapsedOut)=colnames(Peaksdata2_anno_trns_exon)
  CollapsedOut$IDmerge=NA
  
  for (x in 1:nrow(PGene_TBL2)) {
    Trnsc=PGene_TBL2[x,]
    TrnscUL=unlist(unique(str_split(Trnsc$JoinID,pattern = ",")))
    
    d2=Peaksdata2_anno_trns_exon[Peaksdata2_anno_trns_exon$ID%in%TrnscUL,]
    
    if (nrow(d2)==0) {
      next
    } else if (unique(d2$strand)=="+") {
      d5=d2[d2$start%in%min(d2$start),] 
      dmax=d2[d2$Counts_fracMM%in%max(d2$Counts_fracMM),]
      if (nrow(dmax)>1) {dmax=dmax[dmax$start%in%min(dmax$start),]}
    } else if (unique(d2$strand)=="-") {
      d5=d2[d2$start%in%max(d2$start),]
      dmax=d2[d2$Counts_fracMM%in%max(d2$Counts_fracMM),]
      if (nrow(dmax)>1) {dmax=dmax[dmax$start%in%max(dmax$start),]}
    }
    
    ### collapse all columns 
    d3= ( apply(d2 ,2, function(x){paste(unique(x[!is.na(x)]),collapse = ',')}))
    d3[d3==""]<-NA
    cols=c('Counts_Unique','Counts_fracMM','Length')
    d3=as.data.frame(t(d3))
    d3[1,cols]=t(as.data.frame(colSums(d2[,cols])))   
    d3[,cols]=as.numeric(d3[,cols])
    d3$IDmerge=d3$ID
    
    ## Select summary annotation based on Max(dmax) expression peak or 5'(d5) peak   
    danno=dmax # d5 dmax
    cols_AnnotSelect=c('Same_Comb_type_exon','Same_Comb_type_ncRNA','Oppo_Comb_type_exon',
                       'Oppo_Comb_type_ncRNA','Comb_type_exon_Oppo')
    d3[,cols_AnnotSelect]= (danno[,cols_AnnotSelect])
    
    # Select ID for 5' most read or max counts
    d3$start=d5$start
    d3$end=d5$end
    d3$ID=paste0(d5$chr,":",d5$start,"-",d5$end,"_",d5$strand)   
    
    ## if count < read_depth skip
    if (d3$Counts_fracMM<read_depth) {
      next
    } else{
      CollapsedOut=rbind(CollapsedOut[,colnames(CollapsedOut)],d3[1,colnames(CollapsedOut)])
    }
  }  
  CollapsedOut=CollapsedOut[is.na(CollapsedOut$ID)==F,]
  Peaksdata2_anno$IDmerge=NA
  Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$ID%in%Junc_peaks==F,]
  Peaksdata2_anno=(rbind(Peaksdata2_anno,CollapsedOut))
  
  
  ######################################################################
  #### Re filter reads after combing counts of spliced reads
  Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$Counts_fracMM>=read_depth,]
}else{
  
  #### collapse exon
  if (condense_exon==T) {
    print("Collapse exon")
    Peaksdata2_anno_trns=Peaksdata2_anno[(Peaksdata2_anno$Same_feature)%in%""==F,]
    Peaksdata2_anno_trns_exon=Peaksdata2_anno_trns[Peaksdata2_anno_trns$Comb_type_exon_Oppo%in%'protein_coding: exon',]
    trns=Peaksdata2_anno_trns_exon$Same_ensembl_gene_id
    trns=trns[grep(",",trns)]
    trns=c(Peaksdata2_anno_trns_exon$Same_ensembl_gene_id,unlist(strsplit(trns,",")))
    dupGeneName=unique(trns[duplicated(trns)])
    
    if(length(dupGeneName)>0){
      CollapsedOut=as.data.frame(matrix(nrow = 1,ncol=ncol(Peaksdata2_anno_trns_exon)+1));colnames(CollapsedOut)=c(colnames(Peaksdata2_anno_trns_exon),'IDmerge')
      for (x in 1:length(dupGeneName)) {
        d1=dupGeneName[x]
        d2=Peaksdata2_anno_trns_exon[grep(d1,Peaksdata2_anno_trns_exon$Same_ensembl_gene_id),]
        
        if (nrow(d2)==0) {
          next
        } else {
          d2$ID=paste0(d2$chr[1],":",min(d2$start),"-",max(d2$end),"_",unique(d2$strand))
          
          if (unique(d2$strand)=="+") {
            d3=d2[d2$start%in%min(d2$start),]
          } else if (unique(d2$strand)=="-") {
            d3=d2[d2$start%in%max(d2$start),]
          }
          
          d3[,'start']=min(d2$start)
          d3[,'end']=max(d2$end)
          d3$IDmerge=paste(d2$ID,collapse = ",")
          
          trns=trns[grep(",",d2$Same_ensembl_gene_id)]
          d3$Same_ensembl_gene_id=paste0(unique(c(d2$Same_ensembl_gene_id,unlist(strsplit(trns,",")))),collapse = ",")
          trns=trns[grep(",",d2$Same_gene_name_comb)]
          d3$Same_gene_name_comb=paste0(unique(c(d2$Same_gene_name_comb,unlist(strsplit(trns,",")))),collapse = ",")
          trns=trns[grep(",",d2$`Opposite Strand: Host_gene_ensembl_id`)]
          d3$Oppo_ensembl_gene_id=paste0(unique(c(d2$Oppo_ensembl_gene_id,unlist(strsplit(trns,",")))),collapse = ",")
          trns=trns[grep(",",d2$Oppo_gene_name_comb)]
          d3$Oppo_gene_name_comb=paste0(unique(c(d2$Oppo_gene_name_comb,unlist(strsplit(trns,",")))),collapse = ",")
          
          CollapsedOut=rbind(CollapsedOut,d3)
        }
      }
      
      rownames(CollapsedOut)=NULL
      CollapsedOut=CollapsedOut[-1,]
      CollapsedOutID= unique(unlist(strsplit(CollapsedOut$IDmerge,",")))
      Peaksdata2_anno$IDmerge=NA
      Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$ID%in%CollapsedOutID==F,]
      Peaksdata2_anno=(rbind(Peaksdata2_anno,CollapsedOut))
      
    } else { 
      print("Multiple peaks not found in single Gene")
    }
  }  
} 

## correct lnLc annotation - keep lnLc annotation above to hepl track source of lnc annotation
for (x in c('Same_gene_type','Same_gene_type_ALL','Same_type_simple_comb','Same_type_comb',
            'Oppo_gene_type','Oppo_gene_type_ALL','Same_type_simple_comb','Oppo_type_simple_comb','Oppo_type_comb')){
  
  Peaksdata2_anno[,x]=gsub("lnLcRNA","lncRNA",Peaksdata2_anno[,x])
  Peaksdata2_anno[,x]=gsub("linLcRNA","lincRNA",Peaksdata2_anno[,x])
}

#write out for junction annotation 
write.table(Peaksdata2_anno,paste0(out_dir,file_id,'peakannotation_complete.txt'),sep = "\t")

#write out for mapq
write.table(Peaksdata2_anno[,c('chr','start','end','strand','ID')],
            file=paste0(out_dir,file_id,"peakannotation_mapq_IN.txt"), 
            sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")

# print(session_info())
print(nrow(Peaksdata2_anno))

