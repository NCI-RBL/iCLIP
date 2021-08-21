#Purpose:
## R script to generate read length statistics for unaligned and aligned reads
library(Rsamtools,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
library(dplyr)
library(ggplot2)
library(argparse)

#set args
parser <- ArgumentParser()
parser$add_argument("-s","--sampleid", dest="sampleid", required=TRUE, help="Sample id")
parser$add_argument("-ba","--bam_aligned", dest="bam_aligned", required=TRUE, help="Aligned bam file for sample id")
parser$add_argument("-bu","--bam_unaligned", dest="bam_unaligned", required=TRUE, help="Unaligned bam file for sample id")
parser$add_argument("-o","--output_dir", dest="output_dir", required=TRUE, help="output dir location")

args <- parser$parse_args()
sampleid = args$sampleid
bam_aligned = args$bam_aligned
bam_unaligned = args$bam_unaligned
output_dir = args$output_dir

Seq_Processing <- function(txt_in,cat){
  
  #read in seq file
  dat = read.table(txt_in)
  colnames(dat)=c("count","len")
  
  # set up cut-off values 
  breaks <- c(0,15,30,45,60,75,90,105,120,135,150,165)
  
  # specify interval/bin labels
  tags <- c("[0-15)","[15-30)", "[30-45)", "[45-60)", "[60-75)", "[75-90)","[90-105)", "[105-120)","[120-135)",
            "[135-150)","[150-165)")
  
  # bucketing values into bins
  dat$tags <- cut(dat$len, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=FALSE, 
                    labels=tags)
  
  #merge counts
  merged_df=data.frame()
  for (tagid in sort(unique(dat$tags))){
    merged_df[nrow(merged_df)+1,"tags"]=tagid
    merged_df[nrow(merged_df),"counts"]=sum(subset(dat,tags==tagid)$count)
    merged_df[nrow(merged_df),"percent"]=paste0(round(sum(subset(dat,
                                                                 tags==tagid)$count)/sum(dat$count)*100,1),"%")}
  
  #Create and save histogram of counts
  ggplot(data=merged_df, aes(x=tags, y=counts)) +
    geom_bar(stat="identity") +
    xlab('Sequence Length') + ylab('Number of Reads') +
    geom_text(aes(label = percent, y = counts, group = tags), vjust=-.5) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  file_save = paste(output_dir,sampleid,"_",cat,".png",sep="")
  ggsave(file_save,p)
  
  #Create and save text file summary 
  file_save = paste(output_dir,sampleid,"_",cat,".txt",sep="")
  fileConn<-file(file_save)
  line1 = paste("\n* SampleID: ",sampleid,"\n\n")
  line2 = paste("\t + Number of ",cat," reads: ",sum(dat$count),"\n",sep="")
  line3 = paste("\t + Average read length (std): ", format(round(mean(dat$len),2),nsmall=2), 
                " (", format(round(sd(dat$len),2),nsmall = 2),")\n",sep="")
  writeLines(paste(line1,line2,line3,sep=""), fileConn)
  close(fileConn)
}

# MAIN CODE
#process unaligned, aligned  
Seq_Processing(bam_unaligned,"unaligned")
Seq_Processing(bam_aligned,"aligned")
