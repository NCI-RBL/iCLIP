UntNts_CIGAR=function(){
  args <- commandArgs(trailingOnly = TRUE)
  peaks <- args[1] #/data/RBL_NCI/Wolin/Phil/mESC_small_RNA/Output_Custom/peaks/Comb.star_rg_added.sorted.dmark.tlen300.peaks.txt
  bamdir <- args[2] #/data/RBL_NCI/Wolin/Phil/mESC_small_RNA/Output_Custom/bam/
  # samtools<- args[3] #samtools
  # fafile<- args[4] #/data/RBL_NCI/Wolin/Phil/mESC_nascent/reference/fasta/mm10/fdb_igenomes_Mus_musculus_UCSC_mm10_Sequence_WholeGenomeFasta_genome
  bam=args[3] #w3.star_rg_added.sorted.dmark.bam
  # ident_end=args[6]
  # mode_End=args[7]
  # inflect_End=args[8]
  # rmvshrort=args[9]
  # WT_3p_ends=args[10]
  nameadj=args[4]
  Peakmin=args[5]

  samtools='samtools'
  # Rscript UntNts_CIGAR_V6.R ./data/Output_Custom/comb.star_rg_added.sorted.dmark.tlen300.peaks.bed \
  # ../small_RNA_hiseq/bam_novo/ /Users/homanpj/Documents/Tools/samtools/bin/samtools \
  # ./inst/extdata/fdb_igenomes_Mus_musculus_UCSC_mm10_Sequence_WholeGenomeFasta_genome.fa 0 0 w3.star_rg_added.sorted.dmark.bam
  print(peaks)
  print(bamdir)
  print(bam)
  print(nameadj)
  print(args)

  #
  # rm(list=ls())

  library(VariantAnnotation,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(GenomicRanges,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  #library(WhopGenome)
  #library(trackViewer)
  #library(vcfR)
  library(seqinr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(ggplot2,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library("viridis",quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(edgeR,quietly = T,verbose = F)
  library('GenomicFeatures',quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library('rtracklayer',quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  #library(GeneStructureTools)
  library(matrixStats,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(plyr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(tidyr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(fitdistrplus,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(stringr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(data.table)
  library(reshape)

  # library(pacman)
  # if(p_isinstalled(varhandle)==F){ p_install(varhandle) }
  # library(varhandle)
  #






  ################################################################################################################
  ################################################################################################################
  # CIGARadj=1
  # count_dup=1
  # count_multimap=1
  # End5p_count=0
  # unannot_count=0
  # min_nt_cnt_per_pk=-1
  # min_reads_per_nt=0
  # use_samtools=0
  # downsample=0
  # creatUniq=1
  # fig=0

  ################################################################################################################
  ######################################################################################################
  # peaks='./peaks/Ro_Clip_iCountcutadpt_all.unique.NH.mm.ddup.s.unique.bam.peaks.txt'
  # bamdir='./peaks'
  # bam='./peaks/Ro_Clip_iCountcutadpt_all.unique.NH.mm.ddup.s.unique.bam'
  # samtools='/Users/homanpj/Documents/Tools/samtools/bin/samtools'


  #' CIGARadj=args[8]
  # bam=paste(bamdir,bam,sep = "")

  ###############
  bamName=strsplit(bam,"/")[[1]]
  bamName=bamName[length(bamName)]
  bamOut=paste0(bamdir,bamName)

  ### check to see if anything was done to alter original bam file if not create a soft link to original file


  ###############

  peaks=read.table(peaks);
  # colnames(peaks)=c("chr","start","end",'V4','V5',"strand")
  colnames(peaks)=c("chr","start","end",'V4',"strand")
  # peaks=peaks[peaks$chr!="chr_rDNA",]
  peaks=peaks[peaks$V4>=Peakmin,]



  ################################################################################################################
  ######################################################################################################

  col=c('Location','Type','UnannotNT_Start','total_reads')


  outAnnotRNA_3p=as.data.frame(matrix(nrow=nrow(peaks),ncol = 4))
  colnames(outAnnotRNA_3p)=col



  ################################################################################################################
  #### Subset counts for annotated peaks
  ################################################################################################################



  peaks.GR <- GRanges(seqnames = as.character(peaks$chr), ranges=IRanges(start = as.numeric(peaks$start), end = as.numeric(peaks$end)),strand = peaks$strand )



  rmvfiles=function(){
    if(file.exists(paste(bam,".bam",sep = ""))){file.remove(paste(bam,".bam",sep = ""))}
    if(file.exists(paste(bam,".bam.bai",sep = ""))){file.remove(paste(bam,".bam.bai",sep = ""))}
    if(file.exists(paste(bam,".uniq.bam.subset.bam",sep = ""))){file.remove(paste(bam,".uniq.bam.subset.bam",sep = ""))}

    if(file.exists(paste(bam,'.peaktemp.bed',sep=""))){file.remove(paste(bam,'.peaktemp.bed',sep=""))}
    if(file.exists(paste(bam,'.fastatemp.fa',sep=""))){file.remove(paste(bam,'.fastatemp.fa',sep=""))}
    if(file.exists(paste(bam,".SCadj.tmp",sep = ""))){file.remove(paste(bam,".SCadj.tmp",sep = ""))}
    if(file.exists(paste(bam,".SCadj.tmp.bam",sep = ""))){file.remove(paste(bam,".SCadj.tmp.bam",sep = ""))}
    if(file.exists(paste(bam,".SCadj.tmp.bam.bai",sep = ""))){file.remove(paste(bam,".SCadj.tmp.bam.bai",sep = ""))}

    if(file.exists(paste(bam,".subset.sam",sep = ""))){file.remove(paste(bam,".subset.sam",sep = ""))}
    if(file.exists(paste(bam,".subset.bam",sep = ""))){file.remove(paste(bam,".subset.bam",sep = ""))}
    if(file.exists(paste(bam,".subset.bam.bai",sep = ""))){file.remove(paste(bam,".subset.bam.bai",sep = ""))}
    if(file.exists(paste(bam,".subset.adjflag.txt",sep = ""))){file.remove(paste(bam,".subset.adjflag.txt",sep = ""))}
    if(file.exists(paste(bam,".subset.adjflag2.txt",sep = ""))){file.remove(paste(bam,".subset.adjflag2.txt",sep = ""))}
    if(file.exists(paste(bam,".subset.adjflagC.txt",sep = ""))){file.remove(paste(bam,".subset.adjflagC.txt",sep = ""))}
    if(file.exists(paste(bam,".subset.bam.bai",sep = ""))){file.remove(paste(bam,".subset.bam.bai",sep = ""))}
    if(file.exists(paste(bam,".subset.bam.subblast.bam",sep = ""))){file.remove(paste(bam,".subset.bam.subblast.bam",sep = ""))}
    if(file.exists(paste(bam,".subset.bam.subblast.bam.blast.txt",sep = ""))){file.remove(paste(bam,".subset.bam.subblast.bam.blast.txt",sep = ""))}

    if(file.exists(paste(bam,".ntcnt.tmp",sep = ""))){file.remove(paste(bam,".ntcnt.tmp",sep = ""))}
    if(file.exists(paste(bam,".uniq.bam.subset.sam",sep = ""))){file.remove(paste(bam,".uniq.bam.subset.sam",sep = ""))}
    if(file.exists(paste(bam,".headr.tmp",sep = ""))){file.remove(paste(bam,".headr.tmp",sep = ""))}

    if(file.exists(paste(bam,".blast.txt",sep = ""))){file.remove(paste(bam,".blast.txt",sep = ""))}


  }


  for (nx in (1:nrow(peaks))) {
    # for (nx in (nx):nrow(peaks)) {
    l=ls()
    rm(list=l[l%in%c('nx','bam',"bamU",'peaks','nameadj','f','blastn','blastdb','bamOut',
                     'outAnnotRNA_3p','outAnnotRNA_3p_short','outAnnotRNA_3p_long','outAnnotRNA_3p_0','outAnnotRNA_3p_longU','outAnnotRNA_3p_0U',
                     'bamdir','samtools','bedtools','fafile','count_dup','count_multimap','rmvshrort','fig','CIGARadj','ident_end','inflect_End',
                     'mode_End','End5p_count','unannot_count','min_nt_cnt_per_pk','min_reads_per_nt','small_RNA.GR','use_samtools','rmvfiles','alignr','bamQ')==F])
    # print(nx)

    if ((nx/1000)%in%c(1:50000)) {
      print(bam)
      print(paste0(round((nx/nrow(peaks)*100),1)," %"))
      print(nx)}
    #### Select Peak locations
    p=peaks[nx,]

    rmvfiles()

    write.table(p,file=paste(bam,'.peaktemp.bed',sep=""),sep="\t",row.names = F,col.names=F,quote=F)

    pl=paste(p$chr,p$start,sep = ":")
    pl=paste(pl,p$end,sep="-")

    bamout=paste(bam,".out.sam",sep = "")


    ################################################################################################################
    p.GR <- GRanges(seqnames = as.character(p$chr), ranges=IRanges(start = as.numeric(p$start), end = as.numeric(p$end)),strand = as.character(p$strand) )




    ##########################################################################################################
    #### determine if mode was identified and select mode from file
    ################################################################################################################

    # if (ident_end==0){
    #
    #   WT_3p_ends_sub=WT_3p_ends[grep(paste0(p$chr[1],":",p$start[1],"-",p$end[1]),WT_3p_ends$Location ),]
    #
    #   if (nrow(WT_3p_ends_sub)>1) {print("Error too many ends selected"); print(xxxx)}
    #
    #   m=WT_3p_ends_sub$mode
    #   m2=WT_3p_ends_sub[,'mode.1']
    #
    # }

    ##########################################################################################################
    #### Subset BAM
    ################################################################################################################

    samH=paste(samtools,"view -H", bam,' >',paste(bam,".headr.tmp",sep = ""),sep = " ")
    system(samH,wait=T)

    ####################################
    ## select peak region
    ####################################
    ## select only reads that match strand of peak
    if (p$strand=="+") {

      sn=system(paste0(samtools," view -h -c ",bam," ",pl), intern = T)
    }
    if (p$strand=="-") {

      sn=system(paste0(samtools," view -h -c ",bam," ",pl), intern = T)
    }




    outAnnotRNA_3p[nx,1]=paste0(p$chr[1],":",p$start[1],"-",p$end[1],"_",p$strand[1])
    outAnnotRNA_3p[nx,2]=NA
    outAnnotRNA_3p[nx,3]=NA
    outAnnotRNA_3p[nx,4]=sn

    ##################################################################################################


    rmvfiles()

    # if (nx==2) {break}

  }
  # outAnnotRNA_3p=outAnnotRNA_3p[,colnames(outAnnotRNA_3p)%in%c('dist_length','dist_UTnt')==F]
  # outAnnotRNA_3p=outAnnotRNA_3p[is.na(outAnnotRNA_3p[,'UnannotNT_Start'])==F,]
  write.table(outAnnotRNA_3p,file=paste0(bamOut,nameadj,".txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE)


  rmvfiles()
}
UntNts_CIGAR()
