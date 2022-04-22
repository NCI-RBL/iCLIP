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
parser$add_argument("-c2MM","--Smplpeak_bkgroundCount_MM", dest="background_countsMM", required=TRUE, help="Background feature Counts (MultiMapped+Unique) using Sample SAF")
parser$add_argument("-c2uniq","--Smplpeak_bkgroundCount_unique", dest="background_countsUnq", required=TRUE, help="Background feature Counts (Unique) using Sample SAF")
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
background_countsUnq=args$background_countsUnq
background_countsMM=args$background_countsMM

# # testing
# setwd("/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH")
# samplename = "YKO_fclip_wd"
# background = "FLAG_Ro_fclip_wd"
# peak_anno_g1 = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/04_annotation/YKO_fclip_annotation_final_table.txt"
# peak_anno_g2 =  "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/04_annotation/FLAG_Ro_fclip_annotation_final_table.txt"
# pos_manorm =  "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/YKO_fclip_vs_FLAG_Ro_fclip_P/YKO_fclip_wd_MAvalues.xls"
# neg_manorm = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/YKO_fclip_vs_FLAG_Ro_fclip_N/YKO_fclip_wd_MAvalues.xls"
# output_file = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/YKO_fclip_wd_vs_FLAG_Ro_fclip_wd_post_processingTest.txt"
# background_countsUnq="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/counts/FLAG_Ro_fclip_allUnique_YKO_fclip_Peaks.txt"
# background_countsMM="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/counts/FLAG_Ro_fclip_allFracMMCounts_YKO_fclip_Peaks.txt"
# background_countsMMall="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/counts/FLAG_Ro_fclip_allFracMMallCounts_YKO_fclip_Peaks.txt"

# #testing
# setwd("/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH")
# samplename = "FLAG_Ro_fclip_wd"
# background = "YKO_fclip_wd"
# peak_anno_g1 = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/04_annotation/FLAG_Ro_fclip_annotation_final_table.txt"
# peak_anno_g2 =  "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/04_annotation/YKO_fclip_annotation_final_table.txt"
# pos_manorm =  "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/YKO_fclip_vs_FLAG_Ro_fclip_P/FLAG_Ro_fclip_wd_MAvalues.xls"
# neg_manorm = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/YKO_fclip_vs_FLAG_Ro_fclip_N/FLAG_Ro_fclip_wd_MAvalues.xls"
# output_file = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/FLAG_Ro_fclip_wd_vs_YKO_fclip_wd_post_processingTest.txt"
# background_countsUnq="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/counts/YKO_fclip_allUnique_FLAG_Ro_fclip_Peaks.txt"
# background_countsMM="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/counts/YKO_fclip_allFracMMCounts_FLAG_Ro_fclip_Peaks.txt"
# background_countsMMall="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/05_demethod_PH/02_analysis/YKO_fclip_vs_FLAG_Ro_fclip/counts/FLAG_Ro_fclip_allFracMMallCounts_YKO_fclip_Peaks.txt"

############################################################################################################
# Sample1 Processing
SmplA_Peaks=fread(peak_anno_g1, header=T, sep="\t",stringsAsFactors = F,data.table=F)
SmplA_Peaks$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks$`Both Strand: RNA_type_Classification`, 
                                                          levels = c("ncRNA", "protein_coding: exon", "Repeat Element",
                                                                     "pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron",
                                                                     "lncRNA-intron","lncRNA","no Feature"))
####

SmplA_Peaks=tidyr::separate(SmplA_Peaks,col = 'ID',into = c('start','strand'),sep = '_',remove = F)
SmplA_Peaks=tidyr::separate(SmplA_Peaks,col = 'start',into = c('chr','start'),sep = ':',remove = F)
SmplA_Peaks=tidyr::separate(SmplA_Peaks,col = 'start',into = c('start','end'),sep = '-',remove = F)
SmplA_Peaks$start=as.numeric(SmplA_Peaks$start);SmplA_Peaks$end=as.numeric(SmplA_Peaks$end)


SmplA_Peaks$Length=SmplA_Peaks$end-SmplA_Peaks$start


SmplA_Peaks.GR <- GRanges(seqnames = as.character(SmplA_Peaks$chr), ranges=IRanges(start = as.numeric(SmplA_Peaks$start), end = as.numeric(SmplA_Peaks$end)),strand = SmplA_Peaks$strand,ID=SmplA_Peaks$ID )


############################################################################################################
# Sample2 Processing
SmplB_Peaks=fread(peak_anno_g2, header=T, sep="\t",stringsAsFactors = F,data.table=F)
SmplB_Peaks$`Both Strand: RNA_type_Classification`=factor(SmplB_Peaks$`Both Strand: RNA_type_Classification`, 
                                                          levels = c("ncRNA", "protein_coding: exon", "Repeat Element",
                                                                     "pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron",
                                                                     "lncRNA-intron","lncRNA","no Feature"))
####
SmplB_Peaks=tidyr::separate(SmplB_Peaks,col = 'ID',into = c('start','strand'),sep = '_',remove = F)
SmplB_Peaks=tidyr::separate(SmplB_Peaks,col = 'start',into = c('chr','start'),sep = ':',remove = F)
SmplB_Peaks=tidyr::separate(SmplB_Peaks,col = 'start',into = c('start','end'),sep = '-',remove = F)
SmplB_Peaks$start=as.numeric(SmplB_Peaks$start);SmplB_Peaks$end=as.numeric(SmplB_Peaks$end)

SmplB_Peaks$Peak_width=SmplB_Peaks$end-SmplB_Peaks$start
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
# SmplA_Peaks=SmplA_Peaks[,!colnames(SmplA_Peaks)%in%'IDmerge']

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
SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$P_value),c('P_value',paste0('normalized_read_density_in_',samplename),paste0('normalized_read_density_in_',background))]=NA
SmplA_Peaks_MAval$logFC=log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',samplename)]+1)-log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',background)]+1)
SmplA_Peaks_MAval$logPval=-log10(as.numeric(SmplA_Peaks_MAval$P_value))

SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`=factor(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`, levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))


############################################################################################################
## Add background counts 
BKcountMM=fread(background_countsMM, header=T, sep="\t",stringsAsFactors = F,data.table=F)
colnames(BKcountMM)[length(colnames(BKcountMM))]='Counts_Multimappers_Scaled'
BKcountUnq=fread(background_countsUnq, header=T, sep="\t",stringsAsFactors = F,data.table=F)
colnames(BKcountUnq)[length(colnames(BKcountUnq))]='Counts_Unique'

BKcount=merge(BKcountUnq[,!colnames(BKcountMM)%in%c('Chr','Start','End','Strand','Length')],BKcountMM[,!colnames(BKcountMM)%in%c('Chr','Start','End','Strand','Length')],by='Geneid')

SmplA_Peaks_MAval=merge(SmplA_Peaks_MAval,BKcount,by.x='ID',by.y='Geneid',all.x=T,suffixes=paste0('_',c(samplename,background)))


#### combine spliced backround reads

SmplA_Peaks_MAval[(SmplA_Peaks_MAval$IDmerge)%in%"",'Peak_width']=
  SmplA_Peaks_MAval[(SmplA_Peaks_MAval$IDmerge)%in%"",'end']-SmplA_Peaks_MAval[(SmplA_Peaks_MAval$IDmerge)%in%"",'start']


BKcount_dnsmple=BKcount
SmplA_Peaks_MAvalsplc=SmplA_Peaks_MAval[!(SmplA_Peaks_MAval$IDmerge)%in%"",]
splcID=str_split(SmplA_Peaks_MAvalsplc$IDmerge,",")%>%unlist()%>%c(.,SmplA_Peaks_MAvalsplc$ID)%>%unique()
BKcount_dnsmple=BKcount_dnsmple[BKcount_dnsmple$Geneid%in%splcID,]

for (x in 1:nrow(SmplA_Peaks_MAvalsplc)) {
  i=SmplA_Peaks_MAvalsplc[x,]
  
  iL=strsplit(i$IDmerge,",")%>%unlist()
  iL=c(iL,i$ID)%>%unique()
  
  iC=BKcount_dnsmple[BKcount_dnsmple$Geneid%in%iL,]
  BKcount_dnsmple=BKcount_dnsmple[!BKcount_dnsmple$Geneid%in%iL,]
  if (nrow(iC)<0) {Subsample_Counts_Error }
  
  # i[,paste0(c('Counts_Unique_','Counts_Multimappers_Scaled_'),background)]=colSums(iC[,c('Counts_Unique','Counts_Multimappers_Scaled')])
  SmplA_Peaks_MAvalsplc[SmplA_Peaks_MAvalsplc$ID%in%i$ID,paste0(c('Counts_Unique_'),background)]=sum(iC[,c('Counts_Unique')])
  SmplA_Peaks_MAvalsplc[SmplA_Peaks_MAvalsplc$ID%in%i$ID,paste0(c('Counts_Multimappers_Scaled_'),background)]=sum(iC[,c('Counts_Multimappers_Scaled')])
  
  iL=as.data.frame(iL)
  iL=separate(iL,"iL",into = c('chr','start'),sep=":",remove = F)
  iL=separate(iL,"start",into = c('start','strand'),sep="_",remove = T)
  iL=separate(iL,"start",into = c('start','end'),sep="-",remove = T)
  
  iL$length=as.numeric(iL$end)-as.numeric(iL$start)
  SmplA_Peaks_MAvalsplc[SmplA_Peaks_MAvalsplc$ID%in%i$ID,'Peak_width']=sum(iL$length)  
  
}

SmplA_Peaks_MAval=rbind(SmplA_Peaks_MAval[(SmplA_Peaks_MAval$IDmerge)%in%"",],SmplA_Peaks_MAvalsplc)

#################################################

smplnms=unique(c('start','end','Peak_width','Length','P_value','logFC','logPval',paste0(c('Counts_Unique_','Counts_Multimappers_Scaled_','normalized_read_density_in_'),samplename),paste0(c('Counts_Unique_','Counts_Multimappers_Scaled_','normalized_read_density_in_'),samplename)))
SmplA_Peaks_MAval[,smplnms]=sapply(SmplA_Peaks_MAval[,smplnms],as.numeric)




############################################################################################################
#Prepare output
SmplA_Peaks_MAval_out=SmplA_Peaks_MAval[,c("ID","IDmerge","chr","start","end" ,"Peak_width",#"Length",
                                           paste0(c('Counts_Unique_', 'Counts_Multimappers_Scaled_'),samplename),
                                           paste0(c('Counts_Unique_', 'Counts_Multimappers_Scaled_'),background),
                                           'P_value',paste0('normalized_read_density_in_',samplename),paste0('normalized_read_density_in_',background),'logFC',
                                           "Same Strand: Host_gene_name","Same Strand: Host_gene_ensembl_id","Same Strand: RNA_type",'Same Strand: Repeat_Type',
                                           "Same Strand: Intron Exon",'Same Strand: RNA_type_Classification','Same Strand: ncRNA_Classification',"Opposite Strand: Host_gene_name",
                                           "Opposite Strand: Host_gene_ensembl_id","Opposite Strand: RNA_type",'Opposite Strand: Repeat_Type',"Opposite Strand: Intron Exon",
                                           'Opposite Strand: RNA_type_Classification','Opposite Strand: ncRNA_Classification','Both Strand: RNA_type_Classification'
)]

colnames(SmplA_Peaks_MAval_out)[colnames(SmplA_Peaks_MAval_out)%in%c(paste0('normalized_read_density_in_',samplename),
                                                                     paste0('normalized_read_density_in_',background))]=
                                                                    c(paste0(samplename,'-normalized_Count'),
                                                                    paste0(background,'-normalized_Count'))

SmplA_Peaks_MAval_out=dplyr::rename(SmplA_Peaks_MAval_out,"{eval(paste0('ID_',samplename))}":=ID)

SmplA_Peaks_MAval_out=SmplA_Peaks_MAval_out[,colnames(SmplA_Peaks_MAval_out)[!colnames(SmplA_Peaks_MAval_out)%in%c('chr','start','end')]]





write.table(SmplA_Peaks_MAval_out,file=output_file, sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")



# # #############
# # ## testing
# idtest=str_split(SmplA_Peaks_MAval[SmplA_Peaks_MAval$ID%in%'chr1:152902508-152902509_-',"IDmerge"],",")%>%unlist
# 
# # View(BKcount[BKcount$Geneid%in%idtest,])
# 
# ftrCntMM='/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/mES_fclip_1_YL_012122/03_peaks/03_allreadpeaks/YKO_fclip_allFracMMCounts.txt'
# ftrCntMM=fread(ftrCntMM, header=T, sep="\t",stringsAsFactors = F,data.table=F)
# 
# # View(ftrCntMM[ftrCntMM$Geneid%in%idtest,])
# 
# # BKcountMMall=fread(background_countsMMall, header=T, sep="\t",stringsAsFactors = F,data.table=F)
# # colnames(BKcountMMall)[length(colnames(BKcountMMall))]='Counts_MMall'
# # SmplA_Peaks_MAval=merge(SmplA_Peaks_MAval,BKcountMMall[,!colnames(BKcountMMall)%in%c('Chr','Start','End','Strand','Length')],by.x='ID',by.y='Geneid',all.x=T,suffixes=paste0('_',c(samplename,background)))
# # 
# pdata=SmplA_Peaks_MAval[order(SmplA_Peaks_MAval$IDmerge,decreasing =T),]
# 
# 
# pdata$CountslogFC=log2(pdata[,paste0('Counts_Multimappers_Scaled_',samplename)]+1)-log2(pdata[,paste0('Counts_Multimappers_Scaled_',background)]+1)
# pdata$Spliced=pdata$IDmerge
# pdata[!pdata$Spliced%in%"",'Spliced']='Spliced'
# pdata[pdata$Spliced%in%"",'Spliced']='NotSpliced'
# pdata$Spliced=factor(pdata$Spliced,levels=c('Spliced','NotSpliced'))
# pdata=pdata[order(pdata$IDmerge),]
# 
# # ggplot(pdata,aes(x=Counts_MMall,y=`normalized_read_density_in_FLAG_Ro_fclip_wd`,color=Spliced,size=Spliced))+geom_point()+
# #   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
# #   scale_color_manual(values=c('blue','black'))+scale_size_manual(values = c(1,2))
# 
# # ggplot(pdata,aes(x=Counts_MMall,y=Counts_Multimappers_Scaled_FLAG_Ro_fclip_wd,color=Spliced,size=Spliced))+geom_point()+
# #   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
# #   scale_color_manual(values=c('blue','black'))+scale_size_manual(values = c(1,2))
# 
# ggplot(pdata,aes(x=Counts_Multimappers_Scaled_FLAG_Ro_fclip_wd,y=`normalized_read_density_in_FLAG_Ro_fclip_wd`,color=Spliced,size=Spliced))+geom_point()+
#   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
#   scale_color_manual(values=c('blue','black'))+scale_size_manual(values = c(1,2))
# 
# ggplot(pdata,aes(x=logFC,y=CountslogFC,color=Spliced,size=Spliced))+geom_point()+
#   scale_color_manual(values=c('blue','black'))+scale_size_manual(values = c(1,2))
#   # scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
# 
# ggplot(pdata,aes(x=logFC,y=CountslogFC,color=logPval))+geom_point()+
#   scale_colour_gradient(  low = ("blue"),  high = ("red"),na.value = "grey50",guide = "colourbar")
# 
# 
# pdatasplc=pdata[!pdata$IDmerge%in%"",]
# 
# # ggplot(pdatasplc,aes(x=Counts_MMall,y=`normalized_read_density_in_FLAG_Ro_fclip_wd`,color=logPval))+geom_point()+
# #   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
# #   scale_colour_gradient(  low = ("blue"),  high = ("red"),na.value = "grey50",guide = "colourbar")
# # 
# # ggplot(pdatasplc,aes(x=Counts_MMall,y=Counts_Multimappers_Scaled_FLAG_Ro_fclip_wd,color=logPval))+geom_point()+
# #   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
# #   scale_colour_gradient(  low = ("blue"),  high = ("red"),na.value = "grey50",guide = "colourbar")
# 
# ggplot(pdatasplc,aes(x=Counts_Multimappers_Scaled_FLAG_Ro_fclip_wd,y=`normalized_read_density_in_FLAG_Ro_fclip_wd`,color=logPval))+geom_point()+
#   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
#   scale_colour_gradient(  low = ("blue"),  high = ("red"),na.value = "grey50",guide = "colourbar")
# 
# ggplot(pdatasplc,aes(x=logFC,y=CountslogFC,color=logPval))+geom_point()+
#   scale_colour_gradient(  low = ("blue"),  high = ("red"),na.value = "grey50",guide = "colourbar")


