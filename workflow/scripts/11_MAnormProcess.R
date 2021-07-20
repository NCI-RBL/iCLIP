#load libraries
library(tidyr)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(data.table)
library(argparse)

#set args
parser <- ArgumentParser()
parser$add_argument("samplename", action="store_true", default=FALSE, help="Sampleid for sample MAnorm comparion")
parser$add_argument("background", action="store_true", default=FALSE, help="Sampleid for background MAnorm comparion")
parser$add_argument("peak_anno_g1", action="store_true", default=FALSE, help="peakannotation_complete.txt for sample")
parser$add_argument("peak_anno_g2", action="store_true", default=FALSE, help="peakannotation_complete.txt for background")
parser$add_argument("pos_manorm", action="store_true", default=FALSE, help="Positive _MAvalues.xls for sample")
parser$add_argument("neg_manorm", action="store_true", default=FALSE, help="Negative _MAvalues.xls for sample")
parser$add_argument("output_file", action="store_true", default=FALSE, help="output text file name")

args <- parser$parse_args()
samplename = args$samplename
background = args$background
peak_anno_g1 = args$peak_anno_g1
peak_anno_g2 = args$peak_anno_g2
pos_manorm = args$pos_manorm
neg_manorm = args$neg_manorm
output_file = args$output_file

if(samplename==TRUE){
  samplename = "Ro_Clip"
  background = "Control_Clip"
  peak_anno_g1 = "~/../../Volumes/sevillas2/hg38_final/13_annotation/02_peaks/Ro_Clip_peakannotation_complete.txt"
  peak_anno_g2 =  "~/../../Volumes/sevillas2/hg38_final/13_annotation/02_peaks/Control_Clip_peakannotation_complete.txt"
  pos_manorm =  "~/../../Volumes/sevillas2/hg38_final/14_MAnorm/02_analysis/Ro_Clip_Control_Clip_P/Ro_Clip_MAvalues.xls"
  neg_manorm = "~/../../Volumes/sevillas2/hg38_final/14_MAnorm/02_analysis/Ro_Clip_vs_Control_Clip_N/Ro_Clip_MAvalues.xls"
  output_file = "~/../../Volumes/sevillas2/hg38_final/14_MAnorm/02_analysis/Ro_Clip_vs_Control_Clip/RO_Clip_vs_Control_Clip_post_processing.txt"
  
} 
############################################################################################################
# Sample1 Processing
SmplA_Peaks=fread(peak_anno_g1, header=T, sep="\t",stringsAsFactors = F,data.table=F)
SmplA_Peaks$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks$`Both Strand: RNA_type_Classification`, 
                                                          levels = c("ncRNA", "protein_coding: exon", "Repeat Element",
                                                                    "pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron",
                                                                    "lncRNA-intron","lncRNA","no Feature"))
SmplA_Peaks=separate(SmplA_Peaks,col = ID,into = c('chr','start'),sep = ':',remove = F)
SmplA_Peaks=separate(SmplA_Peaks,col = start,into = c('start','end'),sep = '-',remove = F)
SmplA_Peaks$start=as.numeric(SmplA_Peaks$start);SmplA_Peaks$end=as.numeric(SmplA_Peaks$end)
SmplA_Peaks$Length=SmplA_Peaks$end-SmplA_Peaks$start
SmplA_Peaks.GR <- GRanges(seqnames = as.character(SmplA_Peaks$chr), ranges=IRanges(start = as.numeric(SmplA_Peaks$start), end = as.numeric(SmplA_Peaks$end)),strand = SmplA_Peaks$strand,ID=SmplA_Peaks$ID )


############################################################################################################
# Sample2 Processing
SmplB_Peaks=fread(peak_anno_g2, header=T, sep="\t",stringsAsFactors = F,data.table=F)
SmplB_Peaks=separate(SmplB_Peaks,col = ID,into = c('chr','start'),sep = ':',remove = F)
SmplB_Peaks=separate(SmplB_Peaks,col = start,into = c('start','end'),sep = '-',remove = F)
SmplB_Peaks$start=as.numeric(SmplB_Peaks$start);SmplB_Peaks$end=as.numeric(SmplB_Peaks$end)
SmplB_Peaks$Length=SmplB_Peaks$end-SmplB_Peaks$start
SmplB_Peaks.GR <- GRanges(seqnames = as.character(SmplB_Peaks$chr), ranges=IRanges(start = as.numeric(SmplB_Peaks$start), end = as.numeric(SmplB_Peaks$end)),strand = SmplB_Peaks$strand,ID=SmplB_Peaks$ID )

############################################################################################################
#Input positive and negative
Neg_MAval=fread(neg_manorm, header=T, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(Neg_MAval)=gsub("_Neg","",colnames(Neg_MAval))
  Neg_MAval$strand="-"
Pos_MAval=fread(pos_manorm, header=T, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(Pos_MAval)=gsub("_Pos","",colnames(Pos_MAval))
  Pos_MAval$strand="+"

#bind values
MAval=rbind(Neg_MAval,Pos_MAval)
MAval$start=MAval$start-1
MAval$ID=paste0(MAval$chr,':',MAval$start,'-',MAval$end)

#format cols
MAval$Peak_Group=gsub("_Neg","",MAval$Peak_Group)
MAval$Peak_Group=gsub("_Pos","",MAval$Peak_Group)
colnames(MAval)=gsub('_50nt_peakDepth5',"",colnames(MAval))

############################################################################################################
## for spliced peaks select stats from peak with greatest read density
SmplA_PeaksJunc=SmplA_Peaks[(SmplA_Peaks$IDmerge)>1,]
for (x in 1:nrow(SmplA_PeaksJunc)) {
  val=MAval[MAval$ID%in%unlist(strsplit(SmplA_PeaksJunc[x,'IDmerge'],",")),]
  val=val[order(val[,paste0('normalized_read_density_in_',samplename)],decreasing = T),]
  SmplA_PeaksJunc[x,c('P_value',paste0('normalized_read_density_in_',samplename),
                    paste0('normalized_read_density_in_',background))]=val[1,c('P_value',paste0('normalized_read_density_in_',samplename),
                    paste0('normalized_read_density_in_',background))]
  }

SmplA_Peaks_MAval=SmplA_Peaks[SmplA_Peaks$ID%in%SmplA_PeaksJunc$ID==F,]
SmplA_Peaks_MAval=merge(SmplA_Peaks_MAval,MAval[,c("ID",'P_value',paste0('normalized_read_density_in_',samplename), paste0('normalized_read_density_in_',background))],by='ID',all.x=T)
SmplA_Peaks_MAval=rbind(SmplA_Peaks_MAval,SmplA_PeaksJunc)


############################################################################################################
SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$P_value),'P_value']=1e-10
SmplA_Peaks_MAval$FC=log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',samplename)]+1)-log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',background)]+1)
SmplA_Peaks_MAval$logPval=-log10(as.numeric(SmplA_Peaks_MAval$P_value))

SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`, levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))

############################################################################################################
#Prepare output
SmplA_Peaks_MAval_out=SmplA_Peaks_MAval[,c("ID","IDmerge","strand","chr","start","end" ,"Length",
                                           'Counts_Unique', 'Counts_Multimappers_Scaled',
                                           'P_value',paste0('normalized_read_density_in_',samplename),paste0('normalized_read_density_in_',background),'FC',
                                           "Same Strand: Host_gene_name","Same Strand: Host_gene_ensembl_id","Same Strand: RNA_type",'Same Strand: Repeat_Type',
                                           "Same Strand: Intron Exon",'Same Strand: RNA_type_Classification','Same Strand: ncRNA_Classification',"Opposite Strand: Host_gene_name",
                                           "Opposite Strand: Host_gene_ensembl_id","Opposite Strand: RNA_type",'Opposite Strand: Repeat_Type',"Opposite Strand: Intron Exon",
                                           'Opposite Strand: RNA_type_Classification','Opposite Strand: ncRNA_Classification','Both Strand: RNA_type_Classification'
)]

colnames(SmplA_Peaks_MAval_out)[colnames(SmplA_Peaks_MAval_out)%in%c(paste0('normalized_read_density_in_',samplename),
                                paste0('normalized_read_density_in_',background))]=c(paste0(samplename,'-normalized_Count'),paste0(background,'-normalized_Count'))
SmplA_Peaks_MAval_out=SmplA_Peaks_MAval_out[,colnames(SmplA_Peaks_MAval_out)[!colnames(SmplA_Peaks_MAval_out)%in%c('chr','start','end')]]
write.table(SmplA_Peaks_MAval_out,file=output_file, sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")




  







