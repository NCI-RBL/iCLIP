#input peak_junction_annotation
#input dedup.bam

#output mapq text file
#calcmapq true or false
args <- commandArgs(trailingOnly = TRUE)
bam_in = args[1]
output_file = args[9]

#testing
if(is.na(bam_in)){
  bam_in = paste0('/Volumes/data/iCLIP/confirm/TBD/KO_fClip_iCountcutadpt_all.unique.NH.mm.ddup.bam')
  output_file = "testing.csv"
}

glist=function(i){
  PL=i
  GRanges(seqnames = as.character(PL[1]), ranges=IRanges(start = as.numeric(PL[2]), end = as.numeric(PL[3])),strand = PL[4])
}

mapq=function(i){
  out1=as.data.frame(matrix(nrow=1,ncol=2));
  q=i
  s = Bam
  xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))
  sh=as.data.frame(s[xo$subjectHits],row.names = NULL)
  mean(sh$mapq,na.rm = T)
}

if (calcMAPQ==TRUE) {
  bam_param <- ScanBamParam(what=c("mapq", "flag"))
  Bam=readGAlignments(bam_in,
                      use.names=T,
                      index=paste0(bam_in,'.bai'), 
                                   param=bam_param)
  mcols(Bam)$qname=names(Bam)
  
  ###phil theres a prob here - you've commented out the variable Peaksdata2_RP
  #need to add it back
  gout=apply(Peaksdata2_RP[,c('chr','start','end','strand','ID')],1,glist)
  system.time({outl=mclapply(gout,mapq,mc.cores=3)}  )
  
  Peaksdata2_RP$Avg_mapq=unlist(outl )
  Peaksdata2_anno=Peaksdata2_RP
  
  write.table(Peaksdata2_RP[,c('ID','Avg_mapq')],
              file=output_file, 
              sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")
} else if (calcMAPQ==FALSE) {
  ###phil if we aren't running mapq, then we don't need any output, right?
  print ("TBD")
  #CLIP_Peaks_mapq=fread(paste0(inMAPQ,"CLIP_Peaks_mapq.txt"), 
               #         header=T, sep="\t",stringsAsFactors = F,data.table=F)
  #Peaksdata2_anno= merge(Peaksdata2_anno,CLIP_Peaks_mapq,by="ID",all.x = T)
}