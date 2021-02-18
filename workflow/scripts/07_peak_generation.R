library(dplyr)
library(data.table)

testing = "Y"

if (testing=="Y"){
  
  out_dir = "/Volumes/data/iCLIP/marco/14_peaks/"
  anno_dir = "/Volumes/data/iCLIP/marco/15_annotation/"
  
  peaks_in = (paste0(out_dir,"KO_fClip_50nt_peakjunction.txt"))
  
  ref_species="hg38"
  gencode_path = (paste0(anno_dir,"ref_gencode.csv")) 
  refseq_path = (paste0(anno_dir,"ref_refseq.csv")) 
  lncra_path = (paste0(anno_dir,"ref_lncRNA.csv")) 
  alias_path = (paste0(anno_dir,"ref_alias.csv")) 
  
  sample_id = "KO_fClip"
  nt_merge = "50nt"
  
}
##########################################################################################
############### Peak info
##########################################################################################
#read in peaks and alias file generated from 06_peak_junction.R
peaks = read.csv(peaks_in)
alias_anno = read.csv(alias_path)

#merge peaks with ref
peaks_alias = merge(peaks,alias_anno[,c('chr','aliasNCBI2')],by.x='chr',by.y='aliasNCBI2',all.x=T)

#clean
peaks_alias = select(peaks_alias, -c("chr"))

###phil should this be to chr.y? there is no col chr anymore (deleted above)
peaks_alias=peaks_alias[is.na(peaks_alias$chr)==F,]

#rename col
colnames(peaks_alias)=gsub('chr.y','chr',colnames(peaks_alias))

###phil this col is the same as ID - is it necessary? 
#add IDmod
peaks_alias$IDmod=paste0(peaks_alias$chr,":",peaks_alias$start,"-",peaks_alias$end)

###phil why are we switching positive and negative?
peaks_oppo=peaks 
peaks_oppo$strand=gsub("\\+","pos",peaks_oppo$strand) 
peaks_oppo$strand=gsub("\\-","+",peaks_oppo$strand)
peaks_oppo$strand=gsub("pos","-",peaks_oppo$strand)


##################################################################################################
##################################################################################################
############### reference, annotation table
##################################################################################################
#read in references
ref_lncra = read.csv(lncra_path)
ref_refseq = read.csv(ref_refseq)
ref_gencode = read.csv(ref_gencode)

### Add Column that combines all additional annotations
YRNA_rmsk = fread(paste0(anno_dir, "yRNA.bed"))
srpRNA_rmsk = fread(paste0(anno_dir, "srpRNA.bed"))
SKRNA_rmsk = fread(paste0(anno_dir, "SKRNA.bed"))
scRNA_rmsk = fread(paste0(anno_dir, "scRNA.bed"))

#we should think about how to organize this - we should pull in from anno_config? 
#can also call the annotation script into this if we dont need to create BEDS
###phil this assumes that SY flag - is there no option here? if so need to update previous code
annocol=c('chr','start','end','strand','type','name')
if(ref_species == "mm10"){
  #pull in references
  tRNA_sy = fread(paste0(anno_dir, "tRNA.bed"))
  sncRNA_sy = fread(paste0(anno_dir, "sncRNA.bed"))
  rRNA_BK00964 = fread(paste0(anno_dir, "rRNA.bed"))
  
  Anno_RNA_comb=rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],
                      tRNA_sy[,annocol],sncRNA_sy[,annocol],rRNA_BK00964[,annocol],
                      SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],rRNA_rmsk[,annocol])
} else if (ref_species == "hg38"){
  #pull in references
  rRNA_rmsk = fread(paste0(anno_dir, "rRNA.bed"))
  tRNA_rmsk = fread(paste0(anno_dir, "tRNA.bed"))
  
  Anno_RNA_comb = rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],
                        SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],
                        rRNA_rmsk[,annocol],tRNA_rmsk[,annocol])
}

##################################################################################################
############### call peaks
##################################################################################################
#merge references depending on species, run annotation
bam_anno2=function(peaksTable,Annotable,ColumnName,pass_n){
  #testing
  #peaksTable = peaks[,c('chr','start','end','ID','ID2','strand')]
  #anno_dir = "/data/sevillas2/iCLIP/marco/15_annotation/"
  #ColumnName = colSelect
  #Annotable = rbind(ref_gencode[,colMerge],
  #      ref_refseq[,colMerge],
  #      ref_lncra[,colMerge])
  
  ###phil why are we adding multiple chr cols?
  #anno_output = Annotable[,c('chr','start','end','chr','chr','strand',ColumnName)] %>%
  anno_output = Annotable[,c('chr','start','end','strand',ColumnName)] %>%
    rename(
      chr_anno = chr,
      start_anno = start,
      end_anno = end,
      strand_anno = strand
    )
  
  ###phil why are we clearing these?
  #clear chr cols 4,5
  #anno_output[,4:5]=""
  
  #set file prefix
  file_prefix = paste0(sample_id,"_",nt_merge,"_",pass_n)
  
  #output bed files
  write.table(anno_output,
              file = paste0(anno_dir,file_prefix,"annotable.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  write.table(peaksTable,
              file = paste0(anno_dir,file_prefix,"peakstable.bed"), 
              sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  
  #remove lines with illegal characters, remove NAs
  cmd = dQuote("NA")
  system(paste0("cat ",anno_dir,file_prefix,"peakstable.bed | awk '$2 !~ /e/' | awk '$2 != ", cmd, "' > ",
                anno_dir,file_prefix,"peakstable_clean.bed"))
  
  #merge bedfiles into output text
  system(paste0('bedtools intersect -a ', 
                anno_dir, file_prefix,'peakstable_clean.bed -b ', 
                anno_dir, file_prefix,'annotable.bed -wao -s > ',
                out_dir, file_prefix,'peaks_OL.txt'))
  
  ###phil need to figure out where the differences are - must be an error in either the peakstable or annotable
  #read in merged file
  ab_OL=fread(paste0(out_dir,file_prefix,"peaks_OL.txt"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
  colnames(ab_OL)=c(paste0(colnames(peaksTable)),paste0(colnames(anno_output)),'ntOL')
  
  #add width
  ab_OL$width = ab_OL$end-ab_OL$start
  ab_OL$width_anno = ab_OL$end_anno-ab_OL$start_anno
  
  ###phil db is all 0... everything is removed? 
  #subset 
  ab_OL=ab_OL[ab_OL$ntOL>0,]
  
  ab_OL$OLper_anno = ab_OL$ntOL/ab_OL$width_anno
  ab_OL$OLper = ab_OL$ntOL/ab_OL$width
  
  #subset
  ab_OL=ab_OL[(ab_OL$OLper_anno>.75 | ab_OL$OLper>.51), ]
  
  #ID duplicates and uniques
  dup=unique(ab_OL[duplicated(ab_OL$ID),'ID'])
  ab_OL_single=ab_OL[!(ab_OL$ID%in%dup),]
  ab_OL_double=ab_OL[(ab_OL$ID%in%dup),]
  
  u=unique(ab_OL_double$ID)
  ab_OL_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(ab_OL_double)));
  colnames(ab_OL_colapsed)=colnames(ab_OL_double)
  
  #for every unique id
  for(x in 1:length(u)){
    peaksTable = u[x]
    pam = ab_OL_double[ab_OL_double$ID%in%peaksTable,]
    ab_OL_colapsed[x,colnames(ab_OL_double)%in%c(ColumnName)==F]=(pam[1,colnames(ab_OL_double)%in%c(ColumnName)==F,drop=T])
    
    for (cx in 1:length(ColumnName)) {
      ab_OL_colapsed[x,ColumnName[cx]]=paste(sort(unique((pam[,ColumnName[cx]]))), collapse =",")
    }
  }
  
  ab_OL_colapsed=rbind(ab_OL_colapsed,ab_OL_single)
  ab_OL_colapsed=ab_OL_colapsed[!is.na(ab_OL_colapsed$ID),]
  ab_OL_colapsed[(ab_OL_colapsed$ntOL<=0),ColumnName]=NA
  ab_OL=ab_OL_colapsed
  
  return(ab_OL) 
}

peak_calling<-function(peak_in,xopp,nmeprfix){
  colMerge = c('chr','start','end','strand','ensembl_gene_id','transcript_id','external_gene_name','gene_type','gene_type_ALL')
  colSelect = c('ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')
  if (ref_species=='hg38'){
    p=bam_anno2(peaks,
                rbind(ref_gencode[,colMerge],
                      ref_refseq[,colMerge],
                      ref_lncra[,colMerge]), 
                colSelect,
                "pass1")
  } else if (ref_species=='mm10'){
    ###phil removed rbind here - not binding anything - should it be? (1147)
    p=bam_anno2(peaks,
                ref_gencode[,colMerge],
                colSelect,
                "pass1")
  }
  
  #merge peaks_jucntion with p file just created
  ###phil isn't this is redundant? 
  PeaksdataOut = merge(peaks,
                       p[,c('ID','ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')],
                       by='ID',all.x=T)
  
  # PeaksdataOut=bam_anno('RNA_anno',Anno_RNA_comb,PeaksdataOut)
  ### fix variable naming is confusing - back and forth between p and PeaksdataOut
  ###phil same as above... seems redundant especially when merging back?
  p=bam_anno2(PeaksdataOut,
              Anno_RNA_comb,
              c('type','name'),
              "pass2")
  PeaksdataOut = merge(PeaksdataOut,
                       p[,c('ID','type','name')],by='ID',all.x=T)
  
  ###fix p variable swap
  p=PeaksdataOut
  
  ##########################################################################################
  ############### Clean up
  ##########################################################################################
  # change lincRNA,rRNA to rRNA only
  # do not change the name if there is actually a lincRNA name separate from the rRNA grepl(',',p$name)==F)
  p[(grepl(',',p$name)==F)&(p$type%in%'lincRNA-exon,rRNA'),'type']='rRNA'
  p[(grepl(',',p$name)==F)&(p$type%in%'lincRNA-intron,rRNA'),'type']='rRNA'
  p[(grepl(',',p$name)==F)&(p$type%in%'lncRNA,rRNA'),'type']='rRNA'
  
  # change all psueudogene classes to pseudogene
  p[grep('pseudogene',p$gene_type_ALL),'gene_type_ALL']='pseudogene'
  
  # change all sc RNA to yRNA names
  p$type=gsub('scRNA,yRNA','yRNA',p$type)
  
  ###fix another p variable swap
  PeaksdataOut=p
  
  # Combine peak RNAtype info
  ###fix another p variable swap
  p=PeaksdataOut
  
  ###phil can we just remove these cols? 
  p$type_simple_comb=NA
  p$type_comb=NA
  p$gene_name_comb=NA
  
  ##########################################################################################
  ############### duplicates / unique 
  ##########################################################################################
  #handle duplicates
  dup=p[(is.na(p$gene_type_ALL)==F)|(is.na(p$type)==F),'ID']
  peaksTable_single=p[!(p$ID%in%dup),]
  peaksTable_double=p[(p$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  #for each unique id create...?
  ###phil help explain what is being created here
  for(x in 1:length(u)){
    p=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    
    #Comb Annotation
    pam_1=pam$gene_type_ALL ### type from Gencode 
    pam_1=as.data.frame(strsplit(pam_1,','));colnames(pam_1)='a';pam_1$a=as.character(pam_1$a)
    
    pam_2=pam$type ### Bitype from Annotation
    pam_2=as.data.frame(strsplit(pam_2,','));colnames(pam_2)='a';pam_2$a=as.character(pam_2$a)
    
    #remove misc_RNA catagory : these are covered by additional annotations (mostly yRNA)
    if (grepl('misc_RNA',pam_1)) {pam_1[pam_1$a%in%'misc_RNA','a']=NA}
    
    #lincRNA are from an older Gencode version (VM18 or older) so don't use    
    if (grepl('lncRNA',pam_1)&grepl('lincRNA',pam_2)) {pam_1[pam_2$a%in%'lincRNA',]='lncRNA'}
    
    #change ribozyme (gencode) to RNA type from additional anno
    if (grepl('ribozyme',pam_1)) {pam_1[pam_1$a%in%'ribozyme',]=unique(pam_2$a)}
    
    #combine all rna types subtypes into ncRNA
    pam_c=rbind(pam_1,pam_2)
    pam_c=pam_c[!is.na(pam_c$a),,drop=F]
    
    #### Annotation from Gencode : unique(sort(mm10$gene_type_ALL))
    ReplaceType <-function(df_in,col_in,list_in,replace_id){
      for (id in list_in){
        df_in[,col_in]=gsub(id,replace_id,df_in[,col_in])
      } 
      return(df_in)
    }
    
    type_list = c("miRNA","miscRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA",
                  "ribozyme","scRNA", "sRNA","scaRNA", "yRNA", "srpRNA", "tRNA", "7SKRNA",
                  "sncRNA")
    pam_c2 = ReplaceType(pam_c,"a", type_list,"ncRNA")
    
    #phil what if the length isn't greater than 1? keep the full file? 
    if (length(pam_c$a)>1) {pam_c=pam_c[order(pam_c$a),,drop=F]}
    if (length(pam_c2$a)>1) {pam_c2=pam_c2[order(pam_c2$a),,drop=F]}
    
    pam_c=paste(unique(pam_c$a),collapse = ',')
    pam_c2=paste(unique(pam_c2$a),collapse = ',')
    
    #Comb GeneName 
    pam_G1=pam$external_gene_name ### gene name from gencode
    pam_G1=as.data.frame(strsplit(pam_G1,','));colnames(pam_G1)='a';pam_G1$a=as.character(pam_G1$a)
    
    pam_G2=pam$name ### gene name from annotation
    pam_G2=as.data.frame(strsplit(pam_G2,','));colnames(pam_G2)='a';pam_G2$a=as.character(pam_G2$a)
    
    if (grepl('sk',pam_G1)&grepl('7SK',pam_G2)) {
      pam_G2[pam_G2$a%in%'7SK',]=NA
    }
    
    pam_gmerge=rbind(pam_G1,pam_G2)
    pam_gmerge=pam_gmerge[!is.na(pam_gmerge$a),,drop=F]
    pam_gmerge=paste(unique(pam_gmerge$a),collapse = ',')
    
    ### Create Table 
    peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
    peaksTable_colapsed[x,'type_simple_comb']=pam_c2
    peaksTable_colapsed[x,'type_comb']=pam_c
    peaksTable_colapsed[x,'gene_name_comb']=pam_gmerge
    
    remove('pam_1','pam_2','pam_c','pam_c2','pam_G1','pam_G2','pam_gmerge','pam')
  }
  
  ###fix p variable swap
  p=peaksTable_colapsed
  rnatype=p[,c('ID','ensembl_gene_id','external_gene_name','gene_type',"gene_type_ALL",
               'type','name','type_simple_comb','type_comb','gene_name_comb')]
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  colnames(peaksTable_colapsed)[colnames(peaksTable_colapsed) %in% c("ensembl_gene_id",
                                                                     "external_gene_name","gene_type",
                                                                     "gene_type_ALL","type","name",
                                                                     "type_simple_comb","type_comb","gene_name_comb")] =
    paste0(nmeprfix,c("ensembl_gene_id","external_gene_name","gene_type","gene_type_ALL",
                      "type","name","type_simple_comb","type_comb","gene_name_comb"))
  
  ###fix variable switch
  ###phil this is never used again
  PeaksdataOut=peaksTable_colapsed
  
  ###phil why are we creating NA cols?
  peaksTable[,paste0(nmeprfix,'feature')]=NA
  peaksTable[,paste0(nmeprfix,'exon_number')]=NA
  peaksTable[,paste0(nmeprfix,'exon_LargeOL')]=NA
  peaksTable[,paste0(nmeprfix,'intron_number')]=NA
  peaksTable[,paste0(nmeprfix,'intron_LargeOL')]=NA
  peaksTable[,paste0(nmeprfix,'intron_5pStart')]=NA
  peaksTable[,paste0(nmeprfix,'intron_length')]=NA
  peaksTable[,paste0(nmeprfix,'exon_length')]=NA
  return(peaksTable)
}

### fix designations - 1 = "same" 2 = "oppo"? same strand and complementary strand?
peak_same = peak_calling(peaks,1,"same_")
peak_oppo = peak_calling(peaks,2,"oppo_")

##########################################################################################
############### INTRON EXON ANNOTATION

### Identify if CLIP peak overlaps with Intron or Exonic region   

#   Using GTF file from GENCODE v23 -mm10
#   Using GTF file from GENCODE v32 -hg38

# Peaks were annotated by whether they overlap with Host gene intron/exon region
# Intron coordinates were calculated from GTF file.
# 
# A second column was added to idenify if the peak also overlapped with the 5'UTR 3'UTR or CDS (Column: Featrue 2)

##########################################################################################
IE_calling <- function(peak_in,xopp,nmeprfix){
  #create annotation table
  ColumnName = c("feature","exon_number")
  colSelect = c('ensembl_gene_id','external_gene_name','gene_type','gene_type_ALL')
  anno_tmp = intron_exon[grep('protein_coding',intron_exon$gene_type),] %>%
    select(c('chr','start','end','transcript_id','transcript_id','strand',ColumnName))
  colnames(anno_tmp)[colnames(anno_tmp)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')
  
  #run bam_anno 
  exoninof = bam_anno2(PeaksdataOut[,c('chr','start','end','ID','ID2','strand')],
                       anno_tmp,
                       colSelect,
                       "pass3")
  
  ###phil was there a reason for not using bam_anno2?  
  # Annotable=intron_exon[grep('protein_coding',intron_exon$gene_type),]
  # 
  # ###fix variable swap
  # peaksTable=PeaksdataOut
  # 
  ###phil col name is being re-written
  # ColumnName=paste0(c('ensembl_gene_id','external_gene_name'))
  # ColumnName=c("feature","exon_number")
  # 
  # ###fix variable swap
  # p=peaksTable[,c('chr','start','end','ID','ID2','strand')]
  # a=Annotable[,c('chr','start','end','transcript_id','transcript_id','strand',ColumnName)]
  # 
  # ###phil why are we clearing this column 
  # a[,5]=""
  # colnames(a)[colnames(a)%in%c('chr','start','end','strand')]=paste0(c('chr','start','end','strand'),'_anno')
  # 
  # ###phil - why are we writing the same files over again? 
  # write.table(a,file=paste0(out_dir,"/",misc,"/annotable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  # write.table(p,file=paste0(out_dir,"/",misc,"/peakstable.bed"), sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE)
  # 
  # system(paste0('bedtools intersect -a ',gsub(" ","\\\\ ",out_dir),'/',misc,'/peakstable.bed -b ',gsub(" ","\\\\ ",out_dir),'/',misc,'/annotable.bed -wao -s  >',gsub(" ","\\\\ ",out_dir),'/',misc,'/peaks_OL.txt'))
  # 
  # ### fix - writing same files again
  ###phil are there differences between this and the bam_anno2 function? If not we can run together
  # exoninof=fread(paste0(out_dir,"/",misc,"/peaks_OL.txt"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
  # colnames(exoninof)=c(paste0(colnames(p)),paste0(colnames(a)),'ntOL')
  # exoninof=exoninof[exoninof$ntOL>0,]
  # exoninof$width_anno=exoninof$end_anno-exoninof$start_anno
  # exoninof$OLper_anno=exoninof$ntOL/exoninof$width_anno
  # exoninof$width=exoninof$end-exoninof$start
  # exoninof$OLper=exoninof$ntOL/exoninof$width
  # 
  # exoninof=exoninof[(exoninof$OLper_anno>.75 | exoninof$OLper>.51), ]
  
  ###phil what if this is opp? we ignore?
  if (xopp==1) {
    #Calc distance of 5' peak to 5' intron/exondist
    exoninof_pos=exoninof[exoninof$strand%in%'+',]
    exoninof_neg=exoninof[exoninof$strand%in%'-',]
    
    #run positive
    exoninof_pos[,paste0(nmeprfix,'feature_Distance')]=((exoninof_pos$start-exoninof_pos$start_anno)/abs(exoninof_pos$end_anno-exoninof_pos$start_anno))*100
    exoninof_pos[,paste0(nmeprfix,'feature_5pStart')]=exoninof_pos$start_anno
    exoninof_pos[,paste0(nmeprfix,'feature_length')]=abs(exoninof_pos$end_anno-exoninof_pos$start_anno)
    
    #run negative
    exoninof_neg[,paste0(nmeprfix,'feature_Distance')]=((exoninof_neg$end_anno-exoninof_neg$end)/abs(exoninof_neg$end_anno-exoninof_neg$start_anno))*100
    exoninof_neg[,paste0(nmeprfix,'feature_5pStart')]=exoninof_neg$end_anno
    exoninof_neg[,paste0(nmeprfix,'feature_length')]=abs(exoninof_neg$end_anno-exoninof_neg$start_anno)
    
    #bind positive and negative
    exoninof=rbind(exoninof_pos,exoninof_neg);
    
    ###phil why are we pushing NA cols?
    peaksTable[,paste0(nmeprfix,'Exn_start_dist')]=NA
    peaksTable[,paste0(nmeprfix,'Intron_start_dist')]=NA
    peaksTable[,paste0(nmeprfix,'Intron_5pStart')]=NA
    peaksTable[,paste0(nmeprfix,'Exn_5pStart')]=NA
    
  }
}


rpmsk_anno=function(ColumnName,Annotable,peaksTable){
  s=peaksTable
  s=separate(s,ID,into=c('chr','start'),sep=":",remove=F)
  s=separate(s,start,into=c('start','end'),sep="-",remove=F)
  
  s.GR <- GRanges(seqnames = as.character(s$chr), ranges=IRanges(start = as.numeric(s$start), end = as.numeric(s$end)),strand = s$strand,ID=s$ID )
  
  anno.GR <- GRanges(seqnames = as.character(Annotable$genoName), ranges=IRanges(start = as.numeric(Annotable$genoStart), end = as.numeric(Annotable$genoEnd)),strand = Annotable$strand,repClass=Annotable$repClass,repName=Annotable$repName )
  
  q =s.GR
  s=anno.GR
  xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))
  
  qh=as.data.frame(q[xo$queryHits],row.names = NULL)
  sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_repeat")
  rmskinfo=cbind(qh,sh)
  rmskinfo=rmskinfo[,c('ID','repClass_repeat')];colnames(rmskinfo)[colnames(rmskinfo)%in%'repClass_repeat']= ColumnName      
  
  peaksTable=merge(peaksTable[,!colnames(peaksTable)%in%ColumnName],rmskinfo,by='ID',all.x=T)
  
  dup=unique(peaksTable[duplicated(peaksTable$ID),'ID'])
  peaksTable_single=peaksTable[!(peaksTable$ID%in%dup),]
  peaksTable_double=peaksTable[(peaksTable$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  for(x in 1:length(u)){
    peaksTable=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%peaksTable,]
    peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F]=(pam[1,colnames(peaksTable_double)%in%c(ColumnName)==F,drop=T])
    
    peaksTable_colapsed[x,ColumnName]=paste(sort(unique((pam[,ColumnName]))), collapse =" | ")
  }
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  peaksTable_colapsed[is.na(peaksTable_colapsed[,ColumnName]),ColumnName]=NA
  peaksTable_colapsed[peaksTable_colapsed[,ColumnName]%in%"",ColumnName]=NA
  peaksTable=peaksTable_colapsed
  
  return(peaksTable) 
}

### fix designations - 1 = "same" 2 = "oppo"? same strand and complementary strand?
xopp = 1 #same
xopp = 2 #oppo




for (x in 1:nrow(peaksTable)) {
  l=peaksTable[x,c('ID','ID2')]
  g=exoninof[exoninof$ID%in%l$ID,]
  
  if (nrow(g)>0) {
    g_e=g[g$feature%in%'exon',]
    g_i=g[g$feature%in%'intron',]
    
    
    gname=paste(unique(g$feature),collapse = ",")
    gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
    peaksTable[x,paste0(nmeprfix,'feature')]=gsub('NA,',"",gname)
    
    if (nrow(g_e)>0) {
      gname=paste(unique(g_e$exon_number),collapse = ",")
      gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
      peaksTable[x,paste0(nmeprfix,'exon_number')]=gsub('NA,',"",gname)
      
      peaksTable[x,paste0(nmeprfix,'exon_LargeOL')]=max(g_e$ntOL/(g_e$end-g_e$start))
      
      if (xopp==1) {
        g_e[g_e[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
        peaksTable[x,paste0(nmeprfix,'Exn_start_dist')]=mean(as.numeric(g_e[,paste0(nmeprfix,'feature_Distance')]))
      }
    }
    
    if (nrow(g_i)>0) {
      gname=paste(unique(g_i$exon_number),collapse = ",")
      gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
      peaksTable[x,paste0(nmeprfix,'intron_number')]=gsub('NA,',"",gname)
      
      peaksTable[x,paste0(nmeprfix,'intron_LargeOL')]=max(g_i$ntOL/(g_i$end-g_i$start))
      
      if (xopp==1) {
        g_i[g_i[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
        peaksTable[x,paste0(nmeprfix,'Intron_start_dist')]=mean(as.numeric(g_i[,paste0(nmeprfix,'feature_Distance')]))
      }
    }
  }
  if (nrow(g)==0) {
    if (xopp==1) {peaksTable[x,paste0(nmeprfix,c("feature",'exon_number','intron_number','Exn_start_dist','Intron_start_dist'))]=NA}
    if (xopp==2) {peaksTable[x,paste0(nmeprfix,c("feature",'exon_number','intron_number'))]=NA}
  } 
}

############################## to fix this Get introns from refseqNCBI
peaksTable[is.na(peaksTable[,paste0(nmeprfix,'feature')])&peaksTable[,paste0(nmeprfix,'gene_type')]%in%'protein_coding',paste0(nmeprfix,'feature')]='exon'
##############################

peaksTable_inex=peaksTable[grep(',',peaksTable[,paste0(nmeprfix,'feature')]),]

PeaksdataOut=peaksTable

##########################################################################################
############### TO DO 
##########################################################################################

todo<-function(){

    
    #############################################################################################################
    #############################################################################################################
    ### IDENTIFY PEAKS IN REPEAT REGIONS   
    
    # Annotate all repeat regions/Classes identified in Repeatmasker Annotation file (UCSC Table browser)  
    # Data was not filtered based on any of the identified Repeats.  
    #  1) LINE/SINE   
    #  2) LTR   
    #  3) DNA   
    #  4) Satalites   
    #  5) Simple Repeats   
    #  6) Low Complexity   
    #  7) Other   
    #  8) Unknown  
    #############################################################################################################
    #############################################################################################################
    
    PeaksdataOut=PeaksdataOut[!is.na(PeaksdataOut$ID),]
    
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LINE_SINE'),rmsk_GRCm38_LISI,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LTR'),rmsk_GRCm38_LTR,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_DNA'),rmsk_GRCm38_DNA,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Satalites'),rmsk_GRCm38_sat,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Simple_Repeats'),rmsk_GRCm38_SR,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Low_Complexity'),rmsk_GRCm38_LowComplx,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Other'),rmsk_GRCm38_Other,PeaksdataOut)
    PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Unknown'),rmsk_GRCm38_unknown,PeaksdataOut)
    
    
    p=PeaksdataOut
    p[,paste0(nmeprfix,'Repeat_comb')]=NA
    repcol=paste0(nmeprfix,c('Repeat_LINE_SINE','Repeat_LTR','Repeat_DNA','Repeat_Satalites','Repeat_Simple_Repeats','Repeat_Low_Complexity','Repeat_Other','Repeat_Unknown'))
    
    dup=p[rowSums((p[,repcol]>0),na.rm = T)>0,'ID']
    
    peaksTable_single=p[!(p$ID%in%dup),]
    peaksTable_double=p[(p$ID%in%dup),]
    
    u=unique(peaksTable_double$ID)
    peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
    colnames(peaksTable_colapsed)=colnames(peaksTable_double)
    
    
    for(x in 1:length(u)){
      p=u[x]
      pam=peaksTable_double[peaksTable_double$ID%in%p,]
      pam_c=((pam[1,repcol,drop=F]))
      
      pmat=as.data.frame(matrix(nrow=length(pam_c),ncol=1));colnames(pmat)='a'
      pmat$a=t(pam_c)
      pmat=pmat[is.na(pmat)==F,]
      pmat=paste(unique(pmat),collapse = ',')
      
      ######################     
      
      peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
      peaksTable_colapsed[x,paste0(nmeprfix,'Repeat_comb')]=pmat
      
      remove('p','pam','pam_c')
    }
    
    
    peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
    peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
    
    p=peaksTable_colapsed
    PeaksdataOut=peaksTable_colapsed
    
    
    #############################################################################################################
    #############################################################################################################
    ### Asigning Clip peak attributes   
    
    # Not all Peaks overlap with a single feature so peak assignments were assigned by priority:  
    # 
    # ncRNA > Protein coding : Exonic > repeats > Pseudogene > Antisense Feature > Protein Coding : Intronic > lncRNA > no Feature  
    # 
    # All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table.   
    
    #############################################################################################################
    #############################################################################################################
    ### fix with loop
    PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]=NA
    PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')] )
    PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')] )
    PeaksdataOut[,paste0(nmeprfix,'type_comb')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_comb')] )
    PeaksdataOut[,paste0(nmeprfix,'type_comb')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'type_comb')] )
    PeaksdataOut[,paste0(nmeprfix,'gene_type')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type')] )
    PeaksdataOut[,paste0(nmeprfix,'gene_type')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type')] )
    PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')]= gsub('lincRNA',"linLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')] )
    PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')]= gsub('lncRNA',"lnLcRNA",PeaksdataOut[,paste0(nmeprfix,'gene_type_ALL')] )   
    
    # 1. ncRNA        
    comp=( grepl('ncRNA',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) ) 
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='ncRNA'
    
    # 2. protein coding - Exonic
    comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')])&(grepl('exon',PeaksdataOut[,paste0(nmeprfix,'feature')])) )
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='protein_coding: exon'
    
    # 3. repeats
    comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& (is.na(PeaksdataOut[,paste0(nmeprfix,'Repeat_comb')])==F) )
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]="Repeat Element"
    
    # 4. Pseudogene
    comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& grepl('pseudogene',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='pseudogene'
    
    # 5. lncRNA-exon
    comp1=is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])  
    comp2=grepl('linLcRNA-exon',PeaksdataOut[,paste0(nmeprfix,'type_comb')]) | grepl('lnLcRNA',PeaksdataOut[,paste0(nmeprfix,'type_comb')])
    comp=comp1&comp2       
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='lnLcRNA-exon'
    
    # 6. intron
    comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')])& grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='protein_coding: Intron'
    
    # 7. lncRNA-intron
    comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]) & (grepl('linLcRNA-intron',PeaksdataOut[,paste0(nmeprfix,'type_comb')]) ))
    PeaksdataOut[comp,paste0(nmeprfix,'Comb_type_exon')]='lnLcRNA-intron'
    PeaksdataOut[is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]),paste0(nmeprfix,'Comb_type_exon')]='no Feature'
    PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')]=factor(PeaksdataOut[,paste0(nmeprfix,'Comb_type_exon')], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lnLcRNA-exon","Antisense Feature","protein_coding: Intron","lnLcRNA-intron","no Feature"))
    
    ##############################
    #### RNA subtypes
    ##############################
    
    p=PeaksdataOut
    p[,paste0(nmeprfix,'Comb_type_ncRNA')]=NA
    p1=p[p[,paste0(nmeprfix,'Comb_type_exon')]%in%'ncRNA',]
    p2=p[!p[,paste0(nmeprfix,'Comb_type_exon')]%in%'ncRNA',]
    
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=p1[,paste0(nmeprfix,'type_comb')]
    
    ### Protein coding + ncRNA annotations -> ncRNA subtype
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('protein_coding,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(',protein_coding',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    
    ### any ncRNA takes priority over pseudogene
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('pseudogene,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub(',pseudogene',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    
    ### any double annotations with lncRNA become second annotation only
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('lnLcRNA,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[p1[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-intron,linLcRNA-exon',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
    p1[p1[,paste0(nmeprfix,'Comb_type_ncRNA')]%in%'linLcRNA-exon,linLcRNA-intron',paste0(nmeprfix,'Comb_type_ncRNA')]="linLcRNA"
    
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA-exon,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('linLcRNA-intron,',"",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    
    ### fix with loop
    ### any double annotations with yRNA become yRNA only
    p1[grep('yRNA',p1[,paste0(nmeprfix,'Comb_type_ncRNA')]),paste0(nmeprfix,'Comb_type_ncRNA')]='yRNA'
    
    ### tRNA takes priority over miRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,tRNA',"miRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    ### miRNA takes priority over rRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,rRNA',"miRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    ### snRNA takes priority over 7SKRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,snRNA',"snRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('snRNA,7SKRNA',"snRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    ### snRNA takes priority over 7SKRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('7SKRNA,snRNA',"snRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    ### tRNA takes priority over rRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('rRNA,tRNA',"tRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('tRNA,rRNA,',"tRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    ### snoRNA takes priority over miRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,snoRNA',"snoRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('snoRNA,miRNA,',"snoRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    ### scaRNA takes priority over miRNA
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('miRNA,scaRNA',"scaRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    p1[,paste0(nmeprfix,'Comb_type_ncRNA')]=gsub('scaRNA,miRNA,',"scaRNA",p1[,paste0(nmeprfix,'Comb_type_ncRNA')])
    
    PeaksdataOut=rbind(p1,p2)
    
    if (xopp==1) {PeaksdataOut_same=PeaksdataOut}
    if (xopp==2) {PeaksdataOut=merge(PeaksdataOut_same,PeaksdataOut[,colnames(PeaksdataOut)[!colnames(PeaksdataOut)%in%c("chr","start","end","strand","ID2" )]],by='ID')}
    
  }### for pos or neg anno
  
  
  PeaksdataOut$Comb_type_exon_Oppo=NA
  # PeaksdataOut$type_simple_comb
  
  # 1. ncRNA        
  comp=(grepl('ncRNA',PeaksdataOut[,paste0('Same_','type_simple_comb')],fixed = T) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='ncRNA'
  
  # 2. protein coding - Exonic
  comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])&grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')])&(grepl('exon',PeaksdataOut[,paste0('Same_','feature')])) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']=paste0('protein_coding: exon')
  
  # 3. repeats
  comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])&(is.na(PeaksdataOut[,paste0('Same_','Repeat_comb')])==F) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']="Repeat Element"
  
  # 4. Pseudogene
  comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& grepl('pseudogene',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='pseudogene'
  
  # 5. Antisense - non LncRNA
  comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'
  
  # 6. lncRNA-exon
  comp1=is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])
  comp2=((PeaksdataOut[,paste0('Same_','type_comb')]%in%c('linLcRNA-exon')) | (PeaksdataOut[,paste0('Same_','type_comb')]%in%c('lnLcRNA')) )
  comp=comp1&comp2       
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='lnLcRNA-exon'
  
  # 7. Antisense - LncRNA - exon
  comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA-intron',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'
  
  # 8. intron
  comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
  PeaksdataOut[comp,'Comb_type_exon_Oppo']=paste0('protein_coding: Intron')
  
  # 9. Antisense - LncRNA - intron
  comp=( is.na(PeaksdataOut[,'Comb_type_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_type_exon')]%in%'no Feature')==F) & (grepl('lnLcRNA-exon',PeaksdataOut[,paste0('Oppo_','Comb_type_exon')])==F) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='Antisense Feature'
  
  # 10. lncRNA -intron
  comp=(is.na(PeaksdataOut[,'Comb_type_exon_Oppo']) & ((PeaksdataOut[,paste0('Same_','type_comb')]%in%'linLcRNA-intron')==T) )
  PeaksdataOut[comp,'Comb_type_exon_Oppo']='lnLcRNA-intron'
  PeaksdataOut[is.na(PeaksdataOut[,'Comb_type_exon_Oppo']),'Comb_type_exon_Oppo']='no Feature'
  
  
  ### only distances for final annotation_output
  ## distances to any protein coading annotation
  PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Exn_start_dist")]=NA
  PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Intron_start_dist")]=NA
  
  ################################################  
  ### Revert lnLcRNA to lncRNA  
  
  PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
  PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
  PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Same_'),'Comb_type_exon')] )
  PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')]= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,paste0(c('Oppo_'),'Comb_type_exon')] )
  PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-exon','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
  PeaksdataOut[,'Comb_type_exon_Oppo']= gsub('lnLcRNA-intron','lncRNA',PeaksdataOut[,'Comb_type_exon_Oppo'] )
  PeaksdataOut[,'Comb_type_exon_Oppo']=factor(PeaksdataOut[,'Comb_type_exon_Oppo'], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","lncRNA-exon","Antisense Feature","protein_coding: Intron","lncRNA-intron","lncRNA","no Feature"))
  
  unlink(misc, recursive = TRUE)  
  return(PeaksdataOut)
  }
