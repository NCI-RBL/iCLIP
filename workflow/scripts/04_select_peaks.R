#Purpose:
##R script to

#Input (from command line):
## 1) Location of peaks text file ie. /data/sevillas2/Sample.txt
## 2) Dir for output ie. /data/sevillas/sample
## 3) minimum number of counts (int) ie. 2

library(GenomicRanges)
library('rtracklayer')

peaks=function(){
  #Intake args from command line
	args <- commandArgs(trailingOnly = TRUE)
	peaks_input <- args[1]
	out_dir <- args[2]
	minCount <- args [3]

  # Process input table
  peaks_all <- read.table(peaks_input, header=F, sep = "\t", check.names = FALSE,
                        col.names = c("chr", "start", "end", "count", "strand"))
  peaks_all$name <- paste(peaks_all$chr, ":", peaks_all$start, "-", peaks_all$end, sep = "")
  peaks_all$length <- peaks_all$end-peaks_all$start

  #Select loc with counts above minimum
  peaks_select <- peaks_all[peaks_all$count>=minCount,]
  peaks_select$ID <- paste(peaks_select$name,peaks_select$strand,sep = "_")
  peaks_select$Parent <- peaks_select$ID
  peaks_select$feature <- "gene"
  peaks_select.GR <- GRanges(seqnames = peaks_select$chr,
                             ranges=IRanges(start = peaks_select$start, end = peaks_select$end),
                             strand = peaks_select$strand,
                             mcols = peaks_select[,c("ID")],
                             type=peaks_select$feature )

  #out_dirput bed
  rtracklayer::export(peaks_select.GR, paste(out_dir, '.bed', sep = ""), format = 'bed')

  #gff processing
  peaks_gff <- peaks_select
  peaks_gff$feature <- 'exon'
  rownames(peaks_gff) <- NULL
  peaks_gff$name <- NA
  peaks_gff <- rbind(peaks_select,peaks_gff)
  peaks_gff <- peaks_gff[order(peaks_gff$ID,decreasing = T),]
  peaks_gff.GR <- GRanges(seqnames = peaks_gff$chr, ranges=IRanges(start = peaks_gff$start, end = peaks_gff$end),strand = peaks_gff$strand,mcols = peaks_gff[,c("Parent","ID","name")] ,type=peaks_gff$feature)

  #output gtf, gff3
  rtracklayer::export(peaks_gff.GR, paste(out_dir, '.gtf', sep = ""), format = 'gtf')
  rtracklayer::export(peaks_gff.GR, paste(out_dir, '.gff3', sep = ""), format = 'gff3')
}

peaks()
