#Purpose:
## R script to generate read length statistics for unaligned and aligned reads
library(Rsamtools,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
sampleid = args[1]
bam_aligned = args[2]
bam_unaligned = args[3]
output_dir = args[4]

Seq_Processing <- function(txt_in,cat){
  
  #read in seq file
  dat = data.frame(len=scan(txt_in))
  
  # set up cut-off values 
  breaks <- c(0,10,20,30,40,50,60,70,80,90)
  
  # specify interval/bin labels
  tags <- c("[0-10)","[10-20)", "[20-30)", "[30-40)", "[40-50)", "[50-60)","[60-70)", "[70-80)","[80-90)")
  
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
  p_final = p + ggtitle(paste(cat,"sequence lengths:\n",sampleid))
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






##testing
sampleid = "WT_fCLIP"
bam_aligned = "/Volumes/data/iCLIP/marco/fCLIP_HS/00_qc_post/WT_fCLIP_aligned.bam"
bam_unaligned = "/Volumes/data/iCLIP/marco/fCLIP_HS/00_qc_post/WT_fCLIP_unaligned.bam"
output_dir = "/Volumes/sevillas2/git/iCLIP/testing/"
bam_file=bam_aligned
cat="aligned"

BamProcessing <- function(bam_file,cat){
  #Read in bam
  bam <- scanBam(bam_file,index=bam_file)
  
  #aligned workflow
  if(cat=="aligned"){
    
    #Convert seq to df, create length counts
    seq_df <- data.frame(x=bam[[1]]$qwidth,len = bam[[1]]$qwidth)
    remove(bam)
    
    # set up cut-off values 
    breaks <- c(0,10,20,30,40,50,60,70,80,90)
    
    # specify interval/bin labels
    tags <- c("[0-10)","[10-20)", "[20-30)", "[30-40)", "[40-50)", "[50-60)","[60-70)", "[70-80)","[80-90)")
    
    # bucketing values into bins
    group_tags <- cut(seq_df$len, 
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
    p_final = p + ggtitle(paste(cat,"sequence lengths:\n",sampleid))
    
    file_save = paste(output_dir,sampleid,"_",cat,".png",sep="")
    ggsave(file_save,p_final)
    
    #Create and save text file summary 
    file_save = paste(output_dir,sampleid,"_",cat,".txt",sep="")
    fileConn<-file(file_save)
    line1 = paste("\n* SampleID: ",sampleid,"\n\n")
    line2 = paste("\t + Number of ",cat," reads: ",length(bam[[1]]$qwidth),"\n",sep="")
    line3 = paste("\t + Average read length (std): ",format(round(mean(bam[[1]]$qwidth),2),nsmall = 2), 
                  " (", format(round(sd(bam[[1]]$qwidth),2),nsmall = 2),")\n",sep="")
    
    writeLines(paste(line1,line2,line3,sep=""), fileConn)
    close(fileConn)
    
  }
  #unaligned workflow
  else {
    #Convert seq to df, create length counts
    seq_df <- data.frame(x=bam[[1]]$seq,len = nchar(bam[[1]]$seq))
    remove(bam)
    
    # set up cut-off values 
    breaks <- c(0,10,20,30,40,50,60,70,80,90)
    
    # specify interval/bin labels
    tags <- c("[0-10)","[10-20)", "[20-30)", "[30-40)", "[40-50)", "[50-60)","[60-70)", "[70-80)","[80-90)")
    
    # bucketing values into bins
    group_tags <- cut(seq_df$len, 
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
    p_final = p + ggtitle(paste(cat,"sequence lengths:\n",sampleid))

    file_save = paste(output_dir,sampleid,"_",cat,".png",sep="")
    ggsave(file_save,p_final)
    
    #Create and save text file summary 
    file_save = paste(output_dir,sampleid,"_",cat,".txt",sep="")
    fileConn<-file(file_save)
    line1 = paste("\n* SampleID:",sampleid,"\n\n")
    line2 = paste("\t + Number of ",cat," reads: ",nrow(seq_df),"\n",sep="")
    line3 = paste("\t + Average read length (std): ",format(round(mean(seq_df$len),2),nsmall = 2), 
                  " (", format(round(sd(seq_df$len),2),nsmall = 2),")\n",sep="")
    writeLines(paste(line1,line2,line3,sep=""), fileConn)
    close(fileConn)
  }
}

# MAIN CODE
#process unaligned, aligned  
system.time(BamProcessing(bam_unaligned,"unaligned"))
system.time(BamProcessing(bam_aligned,"aligned"))
