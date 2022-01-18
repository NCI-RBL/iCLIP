#load libraries
suppressMessages(library(tidyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library("readxl"))
suppressMessages(library(argparse))

#set args
parser <- ArgumentParser()
parser$add_argument("-s","--samplename", dest="samplename", required=TRUE, help="Sampleid for sample MAnorm comparion")
parser$add_argument("-b","--background", dest="background", required=TRUE, help="Sampleid for background MAnorm comparion")
parser$add_argument("-p1","--peak_anno_g1", dest="peak_anno_g1", required=TRUE, help="peakannotation_complete.txt for sample")
parser$add_argument("-p2","--peak_anno_g2", dest="peak_anno_g2", required=TRUE, help="peakannotation_complete.txt for background")
parser$add_argument("-pos","--pos_manorm", dest="pos_manorm", required=TRUE, help="Positive _MAvalues.xls for sample")
parser$add_argument("-neg","--neg_manorm", dest="neg_manorm", required=TRUE, help="Negative _MAvalues.xls for sample")
parser$add_argument("-o","--output_file", dest="output_file", required=TRUE, help="output text file name")

args <- parser$parse_args()
samplename = args$samplename
background = args$background
peak_anno_g1 = args$peak_anno_g1
peak_anno_g2 = args$peak_anno_g2
pos_manorm = args$pos_manorm
neg_manorm = args$neg_manorm
output_file = args$output_file

# #testing
# setwd("/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/6-22-21-HaCaT_fCLIP")
# samplename = "Y5KO_fCLIP"
# background = "KO_fCLIP"
# peak_anno_g1 = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/6-22-21-HaCaT_fCLIP/13_annotation/Y5KO_fCLIP_peakannotation_final.txt"
# peak_anno_g2 =  "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/6-22-21-HaCaT_fCLIP/13_annotation/KO_fCLIP_peakannotation_final.txt"
# pos_manorm =  "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/6-22-21-HaCaT_fCLIP/14_MAnorm/02_analysis/Y5KO_fCLIPvsKO_fCLIP/Y5KO_fCLIP_Pos_MAvalues.xls"
# neg_manorm = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/6-22-21-HaCaT_fCLIP/14_MAnorm/02_analysis/Y5KO_fCLIPvsKO_fCLIP/Y5KO_fCLIP_Neg_MAvalues.xls"
# output_file = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/6-22-21-HaCaT_fCLIP/14_MAnorm/02_analysis/Y5KO_fCLIPvsKO_fCLIP/Y5KO_fCLIP_vs_KO_fCLIP_post_processing.txt"
# output_file = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/MAnorm_test/Y5KO_fCLIP_vs_KO_fCLIP_post_processing.txt"

############################################################################################################
# Sample1 Processing
SmplA_Peaks=fread(peak_anno_g1, header=T, sep="\t",stringsAsFactors = F,data.table=F)
SmplA_Peaks$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks$`Both Strand: RNA_type_Classification`, 
                                                          levels = c("ncRNA", "protein_coding: exon", "Repeat Element",
                                                                     "pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron",
                                                                     "lncRNA-intron","lncRNA","no Feature"))
SmplA_Peaks=separate(SmplA_Peaks,col = ID,into = c('start','strand'),sep = '_',remove = F)
SmplA_Peaks=separate(SmplA_Peaks,col = start,into = c('chr','start'),sep = ':',remove = F)
SmplA_Peaks=separate(SmplA_Peaks,col = start,into = c('start','end'),sep = '-',remove = F)
SmplA_Peaks$start=as.numeric(SmplA_Peaks$start);SmplA_Peaks$end=as.numeric(SmplA_Peaks$end)
SmplA_Peaks$Length=SmplA_Peaks$end-SmplA_Peaks$start
SmplA_Peaks.GR <- GRanges(seqnames = as.character(SmplA_Peaks$chr), ranges=IRanges(start = as.numeric(SmplA_Peaks$start), end = as.numeric(SmplA_Peaks$end)),strand = SmplA_Peaks$strand,ID=SmplA_Peaks$ID )


############################################################################################################
# Sample2 Processing
SmplB_Peaks=fread(peak_anno_g2, header=T, sep="\t",stringsAsFactors = F,data.table=F)

SmplB_Peaks=separate(SmplB_Peaks,col = ID,into = c('start','strand'),sep = '_',remove = F)
SmplB_Peaks=separate(SmplB_Peaks,col = start,into = c('chr','start'),sep = ':',remove = F)
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
MAval$ID=paste0(MAval$chr,':',MAval$start,'-',MAval$end,'_',MAval$strand)

#format cols
MAval$Peak_Group=gsub("_Neg","",MAval$Peak_Group)
MAval$Peak_Group=gsub("_Pos","",MAval$Peak_Group)
colnames(MAval)=gsub('_50nt_peakDepth5',"",colnames(MAval))

############################################################################################################
## for spliced peaks select stats from peak with greatest read density
SmplA_Peaks=merge(SmplA_Peaks,MAval[,c("ID",'P_value',paste0('normalized_read_density_in_',samplename), paste0('normalized_read_density_in_',background))],by='ID',all.x=T)
SmplA_Peaks=SmplA_Peaks[,!colnames(SmplA_Peaks)%in%'IDmerge']

if ('IDmerge'%in%colnames(SmplA_Peaks)) {
  if (nrow(SmplA_Peaks[(SmplA_Peaks$IDmerge)>1,])>0) {
    

SmplA_PeaksJunc=SmplA_Peaks[(SmplA_Peaks$IDmerge)>1,]
for (x in 1:nrow(SmplA_PeaksJunc)) {
  val=MAval[MAval$ID%in%unlist(strsplit(SmplA_PeaksJunc[x,'IDmerge'],",")),]
  val=val[order(val[,paste0('normalized_read_density_in_',samplename)],decreasing = T),]
  SmplA_PeaksJunc[x,c('P_value',paste0('normalized_read_density_in_',samplename),
                      paste0('normalized_read_density_in_',background))]=val[1,c('P_value',paste0('normalized_read_density_in_',samplename),
                                                                                 paste0('normalized_read_density_in_',background))]
}

SmplA_Peaks_MAval=SmplA_Peaks[SmplA_Peaks$ID%in%SmplA_PeaksJunc$ID==F,]
# SmplA_Peaks_MAval=merge(SmplA_Peaks_MAval,MAval[,c("ID",'P_value',paste0('normalized_read_density_in_',samplename), paste0('normalized_read_density_in_',background))],by='ID',all.x=T)
SmplA_Peaks_MAval=rbind(SmplA_Peaks_MAval,SmplA_PeaksJunc)

  }else{SmplA_Peaks_MAval=SmplA_Peaks;SmplA_Peaks_MAval$IDmerge=NA}
}else{SmplA_Peaks_MAval=SmplA_Peaks;SmplA_Peaks_MAval$IDmerge=NA}




############################################################################################################
SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$P_value),'P_value']=1e-10
SmplA_Peaks_MAval$FC=log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',samplename)]+1)-log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',background)]+1)
SmplA_Peaks_MAval$logPval=-log10(as.numeric(SmplA_Peaks_MAval$P_value))

SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`, levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))

############################################################################################################
#Prepare output
SmplA_Peaks_MAval_out=SmplA_Peaks_MAval[,c("ID","IDmerge","chr","start","end" ,"Length",
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
