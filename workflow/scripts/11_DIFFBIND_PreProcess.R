suppressMessages(library(dplyr))
suppressMessages(library(DiffBind))
suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-s","--samplename", dest="samplename", required=TRUE, help="Sampleid for sample group")
parser$add_argument("-b","--background", dest="background", required=TRUE, help="Sampleid for background group")
parser$add_argument("-d","--input_dir", dest="input_dir", required=TRUE, help="Parent input_dir for CLIP data")
parser$add_argument("-m","--samplemanifest", dest="samplemanifest", required=TRUE, help="path to sample_manifest.tsv")
parser$add_argument("-st","--strand", dest="strand", required=TRUE, help="strand for contrast use either P or N")
parser$add_argument("-o","--output_file", dest="output_file", required=TRUE, help="output text file name")
parser$add_argument("-so","--sample_overlap", dest="sample_overlap", required=F, help="minimum number of sample_df a peak must be in to be included ")

args <- parser$parse_args()
samplename = as.character(args$samplename)
background = as.character(args$background)
input_dir = as.character(args$input_dir)
samplemanifest = as.character(args$samplemanifest)
strand = as.character(args$strand)
output_file = as.character(args$output_file)
sample_overlap = as.integer(args$sample_overlap)

#testing
testing="N"
if(testing=="Y"){
  rm(list=setdiff(ls(), "params"))
  
  input_dir= '/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/2-2-21_fCLIP_duplicates/'
  samplename='WT' # input from contrasts table
  background='KO'
  samplemanifest='sample_manifest.tsv'
  strand='N'
  output_file=paste0(input_dir,'/06_DEP/02_analysis/',samplename,'vs',background,'_DiffBind/',samplename,'vs',background,'_DiffBind','_',strand,'.txt')
  sample_overlap=1
} else if (testing=="SSC"){
  input_dir= "~/../../Volumes/data/diffbind/05_demethod/01_input/"
  samplename='WT'
  background='KO'
  samplemanifest="~/../../Volumes/data/diffbind/sample_manifest.tsv"
  strand='N'
  output_file="~/../../Volumes/data/diffbind/05_demethod/02_analysis/WT_vs_KO/WT_vs_KO_DIFFBIND_N.txt"
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

################################################################################################
## Counting Reads 
################################################################################################
## Identify peaks found in either sample sample_df
dba_out_consensus <- dba.peakset(dba_out,
                                 consensus=c(DBA_CONDITION),
                                 minOverlap=sample_overlap,
                                 bRetrieve=F)

dba_out_consensus <- dba.peakset(dba_out_consensus,
                                 dba_out_consensus$masks[samplename][[1]],
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
## Create output
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

write.table(dba_final,file=output_file, sep = "\t", 
            row.names = FALSE, 
            col.names = T, append = F, quote= FALSE,na = "")