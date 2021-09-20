suppressMessages(library(dplyr))
suppressMessages(library(DiffBind))
suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-s","--samplename", dest="samplename", required=TRUE, help="Sampleid for sample group")
parser$add_argument("-b","--background", dest="background", required=TRUE, help="Sampleid for background group")
parser$add_argument("-d","--input_dir", dest="input_dir", required=TRUE, help="Parent input_dir for CLIP data")
parser$add_argument("-m","--samplemanifest", dest="samplemanifest", required=TRUE, help="path to sample_manifest.tsv")
parser$add_argument("-st","--strand", dest="strand", required=TRUE, help="strand for contrast use either P or N")
parser$add_argument("-ot","--output_table", dest="output_table", required=TRUE, help="output table file name")
parser$add_argument("-os","--output_summary", dest="output_summary", required=TRUE, help="output summary file name")
parser$add_argument("-of","--output_figures", dest="output_figures", required=TRUE, help="output dir for figures with base name")
parser$add_argument("-so","--sample_overlap", dest="sample_overlap", required=F, help="minimum number of sample_df a peak must be in to be included ")

args <- parser$parse_args()
samplename = as.character(args$samplename)
background = as.character(args$background)
input_dir = as.character(args$input_dir)
samplemanifest = as.character(args$samplemanifest)
strand = as.character(args$strand)
output_table = as.character(args$output_table)
output_summary = as.character(args$output_summary)
output_figures = as.character(args$output_figures)
sample_overlap = as.integer(args$sample_overlap)

#testing
testing="SSC"
if(testing=="Y"){
  rm(list=setdiff(ls(), "params"))
  
  input_dir= '/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/2-2-21_fCLIP_duplicates/'
  samplename='WT' # input from contrasts table
  background='KO'
  samplemanifest='sample_manifest.tsv'
  strand='N'
  output_table=paste0(input_dir,'/06_DEP/02_analysis/',samplename,'vs',background,'_DiffBind/',samplename,'vs',background,'_DiffBind','_',strand,'.txt')
  sample_overlap=1
} else if (testing=="SSC"){
  input_dir= "~/../../Volumes/data/diffbind/05_demethod/01_input/"
  samplename='WT'
  background='KO'
  samplemanifest="~/../../Volumes/data/diffbind/sample_manifest.tsv"
  strand='N'
  output_table="~/../../Volumes/data/diffbind/05_demethod/02_analysis/WT_vs_KO/WT_vs_KO_DIFFBINDTable_N.txt"
  output_summary="~/../../Volumes/data/diffbind/05_demethod/02_analysis/WT_vs_KO/WT_vs_KO_DIFFBINDSummary_N.txt"
  output_figures="~/../../Volumes/data/diffbind/05_demethod/02_analysis/WT_vs_KO/figures/WT_vs_KO_DIFFBIND"
  sample_overlap=1
}

#identify strand
if (strand=="P") {strand_sign='+'} else {strand_sign='-'}

#read in sample manifest
sample_df=read.table(samplemanifest,header = T)

################################################################################################
## create diffbind specific df
################################################################################################
#subset df for samples with group id
sample_df=rbind(sample_df[sample_df$group%in%samplename,],
                sample_df[sample_df$group%in%background,])

#create contrast df
contrast_df = as.data.frame(matrix(nrow=nrow(sample_df),ncol=5))%>%
  setNames(c("SampleID", "Condition", "bamReads", "Peaks", "PeakCaller"))

#add sampleid and condition cols
contrast_df[,'SampleID']=sample_df[,'sample']
contrast_df[,'Condition']=sample_df[,'group']

#add bam file location, bed file location
contrast_df[,'bamReads']=paste0(input_dir,sample_df[,'sample'],'_ReadsforDIFFBIND_',strand,'.bam')
contrast_df[,'Peaks']=paste0(input_dir,sample_df[,'sample'],'_PeaksforDIFFBIND_',strand,'.bed')

#add bed as peakcaller
contrast_df[,'PeakCaller']='bed'
################################################################################################
# run diff bind
################################################################################################
dba_out=dba(sampleSheet=contrast_df)

## will automatically detect fragment size
dba_out$config$fragmentSize=0
dba_out$config$minQCth=0
################################################################################################
## Counting Reads 
################################################################################################
## Identify peaks found in either sample sample_df
dba_out_consensus <- dba.peakset(dba_out,
                                 consensus=c(DBA_CONDITION),
                                 minOverlap=sample_overlap,
                                 bRetrieve=F)

#create sub dba for plotting downstream
dba_sub = dba(dba_out_consensus, mask=dba_out_consensus$masks[[samplename]])
dba_sub2 = dba(dba_out_consensus, mask=dba_out_consensus$masks$Consensus)

#run again to retrieve report
dba_out_consensus <- dba.peakset(dba_out_consensus,
                                 dba_out_consensus$masks[[samplename]],
                                 bRetrieve=T)

## Counts can identify all peaks from any sample. May consider modifying command
## below to use all peaks by commenting out peaks command and/or modifiying minOverlap
DBdataCounts <- dba.count(dba_out,
                          peaks = dba_out_consensus,
                          readFormat=DBA_READS_BAM,	                          
                          bRemoveDuplicates = FALSE,
                          bScaleControl = FALSE,
                          bSubControl=FALSE,
                          fragmentSize = 0,
                          minOverlap = sample_overlap,
                          filter = 0,
                          mapQCth=0,
                          minCount=0,
                          summits = FALSE,
                          bParallel=FALSE) %>% suppressWarnings()
################################################################################################
## Normalize data
################################################################################################
DBdatanorm <- dba.normalize(DBdataCounts,
                            method=DBA_DESEQ2,
                            normalize=DBA_NORM_NATIVE,
                            library=DBA_LIBSIZE_PEAKREADS,
                            background=FALSE,
                            control.subtract=FALSE)

#Pull df info from normalized data - not re-running normalization
DBdatanorminfo <- dba.normalize(DBdatanorm, bRetrieve=TRUE)

################################################################################################
## Perform contrast analysis
################################################################################################
## Establishing a model design and contrast
DBdatanormContrast <-dba.contrast(DBdatanorm,
                                  design="~Condition",
                                  contrast=c('Condition',samplename,background),
                                  categories=DBA_CONDITION,
                                  minMembers =2)

#set blacklist and greylist to false
#can be done with dba.analyze(bBlacklist=F, bGreylist=F)
DBdatanormContrast$config$doBlacklist=F
DBdatanormContrast$config$doGreylist=F

#current thresholds are set to:
#DBdatanormContrast$config$minQCth = 15
#DBdatanormContrast$config$mapQCth = 15

#want to include all peaks and reads regardless of quality
DBdatanormContrast$config$minQCth=0
DBdatanormContrast$config$mapQCth=0

#do not run done in parallel using multicore (one process for each contrast 
#for each method, plus an additional process per method)
DBdatanormContrast$config$multicoreInit=F

## Run Diff expression
DB=dba.analyze(DBdatanormContrast, 
               method = DBA_DESEQ2,
               bParallel=F,
               bRetrieveAnalysis=FALSE)

#generate report data
DBAReport <-
  dba.report(DB, 
             th=1, #set to include all sites, regardless of significance
             method=DBA_DESEQ2,
             bNormalized=T,
             bCalled = T,
             bCalledDetail=T,
             bCounts=T) %>% as.data.frame()

################################################################################################
## Create output tables
################################################################################################
#add sign and ID to Reportdf
DBAReport$strand=strand_sign
DBAReport$ID=paste0(DBAReport$seqnames,':',
                    DBAReport$start,'-',
                    DBAReport$end,'_',
                    DBAReport$strand)

#create final df from consensus
dba_final=as.data.frame(dba_out_consensus)
dba_final$strand=strand_sign
dba_final$ID=paste0(dba_final$seqnames,':',
                    dba_final$start,'-',
                    dba_final$end,'_',
                    dba_final$strand)

#remove cols
dba_final=select(dba_final, 
              -c(seqnames,start,end,width,strand))

#merge with DBAReport
dba_final=merge(DBAReport,dba_final,
                by='ID',suffixes=c("","_Occupancy"))

#write stats table
write.table(dba_final,file=output_table, sep = "\t", 
            row.names = FALSE, 
            col.names = T, append = F, quote= FALSE,na = "")

#write summary table
write.table(dba.show(DB),
            file=output_summary, 
            sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")

################################################################################################
## Create output figures
################################################################################################
### VENN DIAGRAMS
# max number of comparisons venn can perform is 4
if(nrow(sample_df[sample_df$group%in%samplename,])<=4){

  ## plot peak overlap
  jpeg(paste0(output_figures,'VennALL_',strand,'.jpeg'))
  dba.plotVenn(dba_sub,main='Sample Group Peak Overlap',
               mask=dba_sub$masks[[samplename]],
               bNotDB=T,bDB=F, bGain=F, bLoss=F, bAll=T)
  dev.off()
  
  ## plot Consensus peak overlap
  jpeg(paste0(output_figures,'VennConsensus_',strand,'.jpeg'))
  dba.plotVenn(dba_sub2,main='Group Concensus Peak Overlap',
               mask=dba_sub2$masks$Consensus,
               bNotDB=T,bDB=F, bGain=F, bLoss=F, bAll=T)
  dev.off()
  
} else{ print("max number of samples were reached, unable to perform Venn comparison")}

### PCA
jpeg(paste0(output_figures,'PCA_',strand,'.jpeg'))
dba.plotPCA(DB,DBA_CONDITION,label=DBA_CONDITION)
dev.off()

### MAplot
jpeg(paste0(output_figures,'MAplot_',strand,'.jpeg'))
dba.plotMA(DB)
dev.off()

### Box Plot
jpeg(paste0(output_figures,'Boxplot_',strand,'.jpeg'))
dba.plotBox(DB)
dev.off()

### Corr Heatmap
jpeg(paste0(output_figures,'Corrplot_',strand,'.jpeg'))
dba.plotHeatmap(DB,bLog=T)
dev.off() 

## Peak Clusters
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
jpeg(paste0(output_figures,'PeakHeatmap_',strand,'.jpeg'))
dba.plotHeatmap(DB, correlations=FALSE,
                scale="row", colScheme = hmap,maxSites=100000,th=1)
dev.off()

jpeg(paste0(output_figures,'PeakHeatmap_sig_',strand,'.jpeg'))
dba.plotHeatmap(DB, correlations=FALSE,contrast = 1,
                scale="row", colScheme = hmap,maxSites=100000,bUsePval = F,th=.01)
dev.off()
