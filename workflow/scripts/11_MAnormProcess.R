
library(tidyr)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
INpeaksAnno = args[1]
INMAnorm = args[2]
contrasts = args[3]
PeakIdnt=args[4]
samplename=args[5]
background=args[6]

if(length(args)==0){
  rm(list=setdiff(ls(), "params"))
  
INpeaksAnno= "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/14_annotation/peaks/"
INMAnorm= "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/15_MAnorm/"
contrasts= '/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/sam_test_master/contrasts.txt'
PeakIdnt= "ALL"
samplename="WT2"
background="KO2"
# compRow= 1
}

## Hardcoded path to markdown. Not sure if you want to make it an option
Rmd='/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/workflow/scripts/12_MAnormAnnotation.Rmd'


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


contrast=fread(contrasts, header=F, sep="\t",stringsAsFactors = F,data.table=F)


# for (compRow in 1:nrow(contrast)) {
# samplename=contrast[compRow,1]
# background=contrast[compRow,2]

INMAnormComp=paste0(INMAnorm,samplename,'vs',background,'/')
PosMAnorm= paste0(INMAnormComp,samplename,"_Pos_MAvalues.xls")
NegMAnorm= paste0(INMAnormComp,samplename,"_Neg_MAvalues.xls")



# ########################### ########################### ########################### ##########################


SmplA_Peaks=fread(paste0(INpeaksAnno,samplename,"_peakannotation_final.txt"), header=T, sep="\t",stringsAsFactors = F,data.table=F)
SmplA_Peaks$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks$`Both Strand: RNA_type_Classification`, levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))


SmplA_Peaks=separate(SmplA_Peaks,col = ID,into = c('chr','start'),sep = ':',remove = F)
SmplA_Peaks=separate(SmplA_Peaks,col = start,into = c('start','end'),sep = '-',remove = F)
SmplA_Peaks$start=as.numeric(SmplA_Peaks$start);SmplA_Peaks$end=as.numeric(SmplA_Peaks$end)
SmplA_Peaks$Length=SmplA_Peaks$end-SmplA_Peaks$start

SmplA_Peaks.GR <- GRanges(seqnames = as.character(SmplA_Peaks$chr), ranges=IRanges(start = as.numeric(SmplA_Peaks$start), end = as.numeric(SmplA_Peaks$end)),strand = SmplA_Peaks$strand,ID=SmplA_Peaks$ID )


# ##########################

SmplB_Peaks=fread(paste0(INpeaksAnno,background,"_peakannotation_final.txt"), header=T, sep="\t",stringsAsFactors = F,data.table=F)

SmplB_Peaks=separate(SmplB_Peaks,col = ID,into = c('chr','start'),sep = ':',remove = F)
SmplB_Peaks=separate(SmplB_Peaks,col = start,into = c('start','end'),sep = '-',remove = F)
SmplB_Peaks$start=as.numeric(SmplB_Peaks$start);SmplB_Peaks$end=as.numeric(SmplB_Peaks$end)
SmplB_Peaks$Length=SmplB_Peaks$end-SmplB_Peaks$start

SmplB_Peaks.GR <- GRanges(seqnames = as.character(SmplB_Peaks$chr), ranges=IRanges(start = as.numeric(SmplB_Peaks$start), end = as.numeric(SmplB_Peaks$end)),strand = SmplB_Peaks$strand,ID=SmplB_Peaks$ID )


##################################################################################################################################

Neg_MAval=fread(paste0( NegMAnorm), header=T, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(Neg_MAval)=gsub("_Neg","",colnames(Neg_MAval))
  Neg_MAval$strand="-"
Pos_MAval=fread(paste0(PosMAnorm), header=T, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(Pos_MAval)=gsub("_Pos","",colnames(Pos_MAval))
  Pos_MAval$strand="+"

MAval=rbind(Neg_MAval,Pos_MAval)
MAval$start=MAval$start-1
MAval$ID=paste0(MAval$chr,':',MAval$start,'-',MAval$end)

MAval$Peak_Group=gsub("_Neg","",MAval$Peak_Group)
MAval$Peak_Group=gsub("_Pos","",MAval$Peak_Group)
# colnames(MAval)=gsub(paste0('_',ntmerge,'|_peakDepth',readdepth,"|",NameAdd,"|Test"),"",colnames(MAval))
colnames(MAval)=gsub('_50nt_peakDepth5',"",colnames(MAval))


##################################################################################################################################
##################################################################################################################################

## for spliced peaks select stats from peak with greatest read density
SmplA_PeaksJunc=SmplA_Peaks[(SmplA_Peaks$IDmerge)>1,]

for (x in 1:nrow(SmplA_PeaksJunc)) {
  
  val=MAval[MAval$ID%in%unlist(strsplit(SmplA_PeaksJunc[x,'IDmerge'],",")),]
  val=val[order(val[,paste0('normalized_read_density_in_',samplename)],decreasing = T),]
  SmplA_PeaksJunc[x,c('P_value',paste0('normalized_read_density_in_',samplename), paste0('normalized_read_density_in_',background))]=val[1,c('P_value',paste0('normalized_read_density_in_',samplename), paste0('normalized_read_density_in_',background))]

  }

SmplA_Peaks_MAval=SmplA_Peaks[SmplA_Peaks$ID%in%SmplA_PeaksJunc$ID==F,]
SmplA_Peaks_MAval=merge(SmplA_Peaks_MAval,MAval[,c("ID",'P_value',paste0('normalized_read_density_in_',samplename), paste0('normalized_read_density_in_',background))],by='ID',all.x=T)

SmplA_Peaks_MAval=rbind(SmplA_Peaks_MAval,SmplA_PeaksJunc)

##################################################################################################################################
##################################################################################################################################


# SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$FC),'FC']=max(SmplA_Peaks_MAval$FC)
SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$P_value),'P_value']=1e-10

SmplA_Peaks_MAval$FC=log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',samplename)]+1)-log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',background)]+1)
SmplA_Peaks_MAval$logPval=-log10(as.numeric(SmplA_Peaks_MAval$P_value))

##########

SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`, levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))

# unique(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`)
# SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`),]
# SmplA_Peaks[is.na(SmplA_Peaks$`Both Strand: RNA_type_Classification`),]



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

SmplA_Peaks_MAval_out=SmplA_Peaks_MAval[,c("ID","IDmerge","strand","chr","start","end" ,"Length",
                                           # "Average MAPQ",
                                           # 'Peak Pvalue',
                                           'Counts_Unique', 'Counts_Multimappers_Scaled',#'Counts_Multimappers' ,'Counts_All', 'Counts_Multimappers_BestMapping', 'Multimaped Reads %',
                                           'P_value',paste0('normalized_read_density_in_',samplename),paste0('normalized_read_density_in_',background),'FC',
                                           "Same Strand: Host_gene_name","Same Strand: Host_gene_ensembl_id","Same Strand: RNA_type",'Same Strand: Repeat_Type',"Same Strand: Intron Exon",'Same Strand: RNA_type_Classification','Same Strand: ncRNA_Classification',"Opposite Strand: Host_gene_name","Opposite Strand: Host_gene_ensembl_id","Opposite Strand: RNA_type",'Opposite Strand: Repeat_Type',"Opposite Strand: Intron Exon",'Opposite Strand: RNA_type_Classification','Opposite Strand: ncRNA_Classification','Both Strand: RNA_type_Classification'#,
                                           # "MM_number","MM_Alt_Peaktype" ,"MM_Same_Peaktype" ,"MM_type","MM_Prec_Alt_Peaktype", "MM_Alt_Alignments"
)]

colnames(SmplA_Peaks_MAval_out)[colnames(SmplA_Peaks_MAval_out)%in%c(paste0('normalized_read_density_in_',samplename),paste0('normalized_read_density_in_',background))]=c(paste0(samplename,'-normalized_Count'),paste0(background,'-normalized_Count'))


SmplA_Peaks_MAval_out=SmplA_Peaks_MAval_out[,colnames(SmplA_Peaks_MAval_out)[!colnames(SmplA_Peaks_MAval_out)%in%c('chr','start','end')]]

write.table(SmplA_Peaks_MAval_out,file=paste0(INMAnormComp,samplename,"vs",background,"_MAnormPeaks.txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")


rmarkdown::render(Rmd,
                  output_file=paste0(INMAnormComp,samplename,"vs",background,"_MAnormPeaks.html"), 
                  params = list(
                    peak_in=INpeaksAnno,
                    peak_in=paste0(INMAnormComp,samplename,"vs",background,"_MAnormPeaks.txt"),
                    contrasts=contrasts,
                    samplename=samplename,
                    background=background,
                    PeakIdnt=PeakIdnt
                  )
)


# }






  







