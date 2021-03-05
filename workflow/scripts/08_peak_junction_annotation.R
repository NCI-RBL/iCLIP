if (params$JoinJunc==TRUE) {
  
  ###recreate original code (line 226)
  #Junc_peaks=unique(c(Junc_PL_merged$ID_Peaks1,Junc_PL_merged$ID_Peaks2))
  Junc_peaks = subset(Peaksdata2,TBD = "TBD")
  
  #select peaks with junction
  Peaksdata2_anno_trns_exon=( Peaksdata2_anno[Peaksdata2_anno$ID%in%Junc_peaks,])
  
  #create new df
  CollapsedOut=as.data.frame(matrix(nrow = (1),ncol = ncol(Peaksdata2_anno_trns_exon)))
  colnames(CollapsedOut)=colnames(Peaksdata2_anno_trns_exon)
  CollapsedOut$IDmerge=NA
  
  ###recreate original code (174)
  #PGene_TBL2=(matrix(nrow=nrow(PGene_TBL_unique),ncol=ncol(PGene_TBL_unique)+4))
  PGene_TBL2 = subset(Peaksdata2,TBD = "TBD")
  
  #for each unique peak
  for (x in 1:nrow(PGene_TBL2)) {
    
    #select row
    Trnsc=PGene_TBL2[x,]
    
    #split ID
    TrnscUL=unlist(unique(str_split(Trnsc$ID_comb,pattern = ",")))
    
    #create df of id annotation information
    d2=Peaksdata2_anno_trns_exon[Peaksdata2_anno_trns_exon$ID%in%TrnscUL,]
    
    #if the line isn't in the annotation file
    ###phil - we should prob do some error handling here; we shouldn't have instances
    #where this is true, unless there are no peaks? and then this whole thing should stop, right? 
    if (nrow(d2)==0) {
      next
    } else{
      print("ERROR")
    }
    
    #Positive strand
    ###fix variables
    if (unique(d2$strand)=="+") {
      d5=d2[d2$start%in%min(d2$start),] 
      dmax=d2[d2$Counts_fracMM%in%max(d2$Counts_fracMM),]
      
      ###phil need error  handling here
      if (nrow(dmax)>1) {
        dmax=dmax[dmax$start%in%min(dmax$start),]
      } else{
        print ("TBD")
      }
      #negative strand
    } else if (unique(d2$strand)=="-") {
      d5=d2[d2$start%in%max(d2$start),]
      dmax=d2[d2$Counts_fracMM%in%max(d2$Counts_fracMM),]
      
      ###phil need error handling here
      if (nrow(dmax)>1) {
        dmax=dmax[dmax$start%in%max(dmax$start),]
      } else{
        print("TBD")
      }
    }
    
    ### collapse all columns 
    d3 = (apply(d2 ,2, function(x){paste(unique(x[!is.na(x)]),collapse = ',')}))
    d3[d3==""]<-NA
    
    cols=c('Counts_Unique','Counts_fracMM','Length')
    d3=as.data.frame(t(d3))
    d3[1,cols] = t(as.data.frame(colSums(d2[,cols])))   
    d3[,cols] = as.numeric(d3[,cols])
    d3$IDmerge=d3$ID
    
    ## Select summary annotation based on Max(dmax) expression peak or 5'(d5) peak   
    danno=dmax # d5 dmax
    cols_AnnotSelect=c('Same_Comb_type_exon', 'Same_Comb_type_ncRNA',
                       'Oppo_Comb_type_exon','Oppo_Comb_type_ncRNA','Comb_type_exon_Oppo')
    d3[,cols_AnnotSelect]= (danno[,cols_AnnotSelect])
    
    # Select ID for 5' most read or max counts
    d3$start=d5$start
    d3$end=d5$end
    d3$ID=paste0(dmax$chr,":",dmax$start,"-",dmax$end)   
    
    ## if count < readdepth skip
    if (d3$Counts_fracMM<params$readdepth) {
      next
    }
    
    CollapsedOut=rbind(CollapsedOut[,colnames(CollapsedOut)],d3[1,colnames(CollapsedOut)])
  }
  
  #remove rows without id
  CollapsedOut = complete.cases(CollapsedOut[,"ID"])
  
  ###fix variable we dont want to overwrite the original
  #subset for all non-junction peaks
  Peaksdata2_anno$IDmerge=NA
  Peaksdata2_anno=subset(Peaksdata2_anno, !(ID%in%Junc_peaks))
  
  #merge non-junction peaks (Peaksdata2_anno) with junction peaks (CollapsedOut)
  Peaksdata2_anno=(rbind(Peaksdata2_anno,CollapsedOut))
  
  #### Re filter reads after combing counts of spliced reads
  Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$Counts_fracMM>=params$readdepth,]
  
} else if (params$JoinJunc==FALSE) {
  
  ###phil need to fill else statement
  # collapse exon
  if (params$Condense==TRUE) {
    Peaksdata2_anno_trns_exon = subset(Peaksdata2_anno, !(Same_feature == "") &
                                         Comb_type_exon_Oppo == 'protein_coding: exon')
    #create list
    trns=Peaksdata2_anno_trns_exon$Same_ensembl_gene_id
    
    ###fix variable naming
    trns=trns[grep(",",trns)]
    trns=c(Peaksdata2_anno_trns_exon$Same_ensembl_gene_id,unlist(strsplit(trns,",")))
    
    #remove dups
    dupGeneName=unique(trns[duplicated(trns)])
    
    #create df
    CollapsedOut = as.data.frame(matrix( nrow = 1,
                                         ncol = ncol(Peaksdata2_anno_trns_exon)+1))
    colnames(CollapsedOut)=c(colnames(Peaksdata2_anno_trns_exon),'IDmerge')
    
    #for each unique val
    for (x in 1:length(dupGeneName)) {
      d1=dupGeneName[x]
      d2=Peaksdata2_anno_trns_exon[grep(d1,Peaksdata2_anno_trns_exon$Same_ensembl_gene_id),]
      d2ID=paste0(d2$chr[1],":",min(d2$start),"-",max(d2$end))
      
      #positive strand
      if (unique(d2$strand)=="+") {
        d3=d2[d2$start%in%min(d2$start),]
        #negative strand
      } else if (unique(d2$strand)=="-") {
        d3=d2[d2$start%in%max(d2$start),]
      }
      
      d3[,'start']=min(d2$start)
      d3[,'end']=max(d2$end)
      d3$IDmerge=paste(d2$ID,collapse = ",")
      
      trns=trns[grep(",",d2$Same_ensembl_gene_id)]
      d3$Same_ensembl_gene_id = paste0(unique(c(d2$Same_ensembl_gene_id,unlist(strsplit(trns,",")))),
                                       collapse = ",")
      trns=trns[grep(",",d2$Same_gene_name_comb)]
      d3$Same_gene_name_comb=paste0(unique(c(d2$Same_gene_name_comb,unlist(strsplit(trns,",")))),
                                    collapse = ",")
      trns=trns[grep(",",d2$`Opposite Strand: Host_gene_ensembl_id`)]
      d3$Oppo_ensembl_gene_id=paste0(unique(c(d2$Oppo_ensembl_gene_id,unlist(strsplit(trns,",")))),
                                     collapse = ",")
      trns=trns[grep(",",d2$Oppo_gene_name_comb)]
      d3$Oppo_gene_name_comb=paste0(unique(c(d2$Oppo_gene_name_comb,unlist(strsplit(trns,",")))),
                                    collapse = ",")
      CollapsedOut=rbind(CollapsedOut,d3)
    }
    
    #remove rownames, col
    rownames(CollapsedOut)=NULL
    CollapsedOut=CollapsedOut[-1,]
    CollapsedOutID= unique(unlist(strsplit(CollapsedOut$IDmerge,",")))
    
    ###fix we dont want to overwrite original df
    Peaksdata2_anno$IDmerge=NA
    Peaksdata2_anno=Peaksdata2_anno[Peaksdata2_anno$ID%in%CollapsedOutID==F,]
    Peaksdata2_anno=(rbind(Peaksdata2_anno,CollapsedOut))
  }
}

#write out for junction annotation 
write.csv(Peaksdata2_anno,"peakjunction_annotation.csv")

#write out for mapq
write.table(Peaksdata2_anno[,c('chr','start','end','strand','ID')],
            file=paste0(params$OutDIR,"/MAPQ/CLIP_",samplename,"_",
                        params$PeakIdnt,"Peaks",ntmerge,'_peakDepth',
                        params$readdepth,params$NameAdd,"_mapq_IN.txt"), 
            sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE,na = "")