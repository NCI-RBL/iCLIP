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
  dat = data.frame(len=scan(txt_in))
  
  # set up cut-off values 
  breaks <- c(0,15,30,45,60,75,90,105,120,135,150,165)
  
  # specify interval/bin labels
  tags <- c("[0-15)","[15-30)", "[30-45)", "[45-60)", "[60-75)", "[75-90)","[90-105)", "[105-120)","[120-135)",
            "[135-150)","[150-165)")
  
  # bucketing values into bins
  group_tags <- cut(dat$len, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=FALSE, 
                    labels=tags)
  
  #Create and save histogram of lengths
  p = ggplot(data = as_tibble(group_tags), mapping = aes(x=value)) + 
    geom_bar(fill="blue",color="white",alpha=0.7) + 
    stat_count(geom="text", aes(label=..count..), vjust=-0.5) +
    xlab('Sequence Length') + ylab('Number of Reads') +
    theme_minimal()
  p_final = p #+ ggtitle(paste(cat,"sequence lengths:",sampleid))
  file_save = paste(output_dir,sampleid,"_",cat,".png",sep="")
  ggsave(file_save,p_final)
  
  #Create and save text file summary 
  file_save = paste(output_dir,sampleid,"_",cat,".txt",sep="")
  fileConn<-file(file_save)
  line1 = paste("\n* SampleID: ",sampleid,"\n\n")
  line2 = paste("\t + Number of ",cat," reads: ",nrow(dat),"\n",sep="")
  line3 = paste("\t + Average read length (std): ", format(round(mean(dat$len),2),nsmall=2), 
                " (", format(round(sd(dat$len),2),nsmall = 2),")\n",sep="")
  writeLines(paste(line1,line2,line3,sep=""), fileConn)
  close(fileConn)
}

# MAIN CODE
#process unaligned, aligned  
Seq_Processing(bam_unaligned,"unaligned")
Seq_Processing(bam_aligned,"aligned")
