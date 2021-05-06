args
sample_peaks = argv[1]
background_peaks = argv[2]
manorm_pos = argv[3]
manorm_neg = argv[4]
output_file = argv[5]

if(is.na(sample_peaks)){
  sample_peaks = "WT_CLIP_Peaks_50nt_peakDepth_MD15.txt"
  background_peaks =  "KO_CLIP_Peaks_50nt_peakDepth_MD15.txt"
  manorm_pos = "WT_50nt_peakDepth5Test_Pos_MAvalues.xls"
  manorm_neg = "WT_50nt_peakDepth5Test_Neg_MAvalues.xls"
  output_file = "WT_KO_manorm_values.txt"
}

###phil - rewritten below
# SmplA_Peaks = fread(sample_peaks, header=T, sep="\t",stringsAsFactors = F,data.table=F)
# SmplA_Peaks$`Both Strand: RNA_type_Classification`= factor(SmplA_Peaks$`Both Strand: RNA_type_Classification`, 
#                                                            levels = c("ncRNA", "protein_coding: exon", 
#                                                                       "Repeat Element","pseudogene",
#                                                                       "lncRNA-exon","Antisense Feature",
#                                                                       "protein_coding: Intron","lncRNA-intron",
#                                                                       "lncRNA","no Feature"))
# SmplA_Peaks = separate(SmplA_Peaks, col = ID,into = c('chr','start'),sep = ':',remove = F) %>%
#   separate(col = start,into = c('start','end'),sep = '-',remove = F)
# 
# SmplA_Peaks$start = as.numeric(SmplA_Peaks$start)
# SmplA_Peaks$end = as.numeric(SmplA_Peaks$end)
# SmplA_Peaks$Length = SmplA_Peaks$end-SmplA_Peaks$start
# 
# SmplA_Peaks.GR <- GRanges(seqnames = as.character(SmplA_Peaks$chr), 
#                           ranges=IRanges(start = as.numeric(SmplA_Peaks$start), 
#                                          end = as.numeric(SmplA_Peaks$end)),
#                           strand = SmplA_Peaks$strand,ID=SmplA_Peaks$ID)
# 
# 
# SmplB_Peaks = fread(background_peaks, header=T, sep="\t",stringsAsFactors = F,data.table=F)
# 
# SmplB_Peaks=separate(SmplB_Peaks,col = ID,into = c('chr','start'),sep = ':',remove = F) %>%
#   separate(SmplB_Peaks,col = start,into = c('start','end'),sep = '-',remove = F)
# 
# SmplB_Peaks$start=as.numeric(SmplB_Peaks$start);SmplB_Peaks$end=as.numeric(SmplB_Peaks$end)
# SmplB_Peaks$Length=SmplB_Peaks$end-SmplB_Peaks$start
# 
# SmplB_Peaks.GR <- GRanges(seqnames = as.character(SmplB_Peaks$chr), 
#                           ranges=IRanges(start = as.numeric(SmplB_Peaks$start), 
#                                          end = as.numeric(SmplB_Peaks$end)),
#                           strand = SmplB_Peaks$strand,ID=SmplB_Peaks$ID )


GeneratePeakGR<- function(peak_in){
  sample_peaks = fread(peak_in, header=T, sep="\t",stringsAsFactors = F,data.table=F)
  
  ###phil - this isn't included in background code- do we need it?
  sample_peaks$`Both Strand: RNA_type_Classification`= factor(sample_peaks$`Both Strand: RNA_type_Classification`, 
                                                             levels = c("ncRNA", "protein_coding: exon", 
                                                                        "Repeat Element","pseudogene",
                                                                        "lncRNA-exon","Antisense Feature",
                                                                        "protein_coding: Intron","lncRNA-intron",
                                                                        "lncRNA","no Feature"))
  
  #separate ID and start
  sample_peaks=separate(sample_peaks,col = ID,into = c('chr','start'),sep = ':',remove = F) %>%
    separate(sample_peaks,col = start,into = c('start','end'),sep = '-',remove = F)
  
  ###phil redundant in GR
  sample_peaks$start=as.numeric(sample_peaks$start)
  sample_peaks$end=as.numeric(sample_peaks$end)
  sample_peaks$Length=sample_peaks$end-sample_peaks$start
  
  #generate GRanges object
  sample_peaks <- GRanges(seqnames = as.character(sample_peaks$chr), 
                            ranges=IRanges(start = as.numeric(sample_peaks$start), 
                                           end = as.numeric(sample_peaks$end)),
                            strand = sample_peaks$strand,ID=sample_peaks$ID )
  
  return(sample_peaks)
}

###fix rename sample_a to sample and sample_b to background OR change naming upstream to match
#read in peak files, create GR object
SmplA_Peaks.Gr = GeneratePeakGR(sample_peaks)
SmplB_Peaks.Gr = GeneratePeakGR(backgrond_peaks)

#read in positive and neg MANORM 
Neg_MAval=fread(manorm_neg, header=T, sep="\t",stringsAsFactors = F,data.table=F)
colnames(Neg_MAval)=gsub("_Neg","",colnames(Neg_MAval))
Neg_MAval$strand="-"

Pos_MAval=fread(paste0(PosMAnorm), header=T, sep="\t",stringsAsFactors = F,data.table=F)
colnames(Pos_MAval)=gsub("_Pos","",colnames(Pos_MAval))
Pos_MAval$strand="+"

#bind positive and neg
MAval=rbind(Neg_MAval,Pos_MAval)
MAval$start=MAval$start-1

#generate new id
MAval$ID=paste0(MAval$chr,':',MAval$start,'-',MAval$end)

###phil - i dont think this is necessary? We dont add this to the variables, only the columns
MAval$Peak_Group=gsub("_Neg","",MAval$Peak_Group)
MAval$Peak_Group=gsub("_Pos","",MAval$Peak_Group)

###phil this maybe not needed - test col?
colnames(MAval)=gsub(paste0('_',ntmerge,'_peakDepth',params$readdepth,params$NameAdd),"",colnames(MAval))
colnames(MAval)=gsub(paste0('_',ntmerge,'_peakDepth',params$readdepth,'Test'),"",colnames(MAval))

###phil - the MAval is already a merged SmplA_Peaks and SmpleB_Peaks - why are we merging again?
SmplA_Peaks_MAval=merge(SmplA_Peaks,MAval,by='ID',suffixes=c("","_MAnorm"),all.x=T)

SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$FC),'FC']=max(SmplA_Peaks_MAval$FC)
SmplA_Peaks_MAval[is.na(SmplA_Peaks_MAval$P_value),'P_value']=1e-10

SmplA_Peaks_MAval$FC = log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',samplename)]+1)-log2(SmplA_Peaks_MAval[,paste0('normalized_read_density_in_',background)]+1)
SmplA_Peaks_MAval$logPval=-log10(as.numeric(SmplA_Peaks_MAval$P_value))

###phil this col is added in the GR formation - we can remove above  
c$`Both Strand: RNA_type_Classification` = factor(SmplA_Peaks_MAval$`Both Strand: RNA_type_Classification`, 
                                                                  levels = c("ncRNA", "protein_coding: exon",
                                                                             "Repeat Element","pseudogene",
                                                                             "lncRNA-exon","Antisense Feature",
                                                                             "protein_coding: Intron","lncRNA-intron",
                                                                             "lncRNA","no Feature"))

#write out file
write.csv(SmplA_Peaks_MAval,output_file)