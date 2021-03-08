#library
library(tidyr)
library(GenomicRanges)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
peak_type = args[1]
peak_unique = args[2]
peak_all = args[3]
join_junction = args[4]
read_depth = args[5]
DEmethod = args [6]
sample_id = args[7]
nt_merge = args[8]
out_dir = args[9]

testing="Y"
if(testing=="Y"){
  peak_type= "ALL"
  peak_unique = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/NextSeq_fCLIP_0218/13_counts/uniquereadpeaks/WT1_fCLIP_50nt_uniqueCounts.txt"
  peak_all = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/NextSeq_fCLIP_0218/13_counts/uniquereadpeaks/WT1_fCLIP_50nt_allFracMMCounts.txt"
  join_junction = "TRUE"
  read_depth = "5"
  DEmethod = "MANORM"
  sample_id = "WT1"
  nt_merge = "50nt"
  out_dir = "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/NextSeq_fCLIP_0218/junction/"
}

##########################################################################################
############### unique reads
##########################################################################################
FtrCount_uniq=read.delim(peak_unique, header=T, sep="\t",
                         stringsAsFactors = F,comment.char = '#')
colnames(FtrCount_uniq)[7]='Counts'
FtrCount_uniq$Start=as.numeric(FtrCount_uniq$Start)
FtrCount_uniq$End=as.numeric(FtrCount_uniq$End)
FtrCount_uniq=FtrCount_uniq[!is.na(FtrCount_uniq$Start),]
FtrCount_uniq$ID=paste0(FtrCount_uniq$Chr,":",FtrCount_uniq$Start,"-",FtrCount_uniq$End)
FtrCount_uniq$ID2=paste0(FtrCount_uniq$Chr,":",FtrCount_uniq$Start,"-",FtrCount_uniq$End,"_",FtrCount_uniq$Strand)

##########################################################################################
############### all reads
##########################################################################################
FtrCount_frac=read.delim(peak_all, header=T, sep="\t",
                         stringsAsFactors = F,comment.char = '#')
colnames(FtrCount_frac)[7]='Counts'
FtrCount_frac$Start=as.numeric(FtrCount_frac$Start)
FtrCount_frac$End=as.numeric(FtrCount_frac$End)
FtrCount_frac=FtrCount_frac[!is.na(FtrCount_frac$Start),]
FtrCount_frac$ID=paste0(FtrCount_frac$Chr,":",FtrCount_frac$Start,"-",FtrCount_frac$End)

##########################################################################################
############### merge
##########################################################################################
FtrCount_merged=merge(FtrCount_uniq[,(colnames(FtrCount_uniq) %in% "Geneid")==F],
               FtrCount_frac[,c('ID','Counts')],
               by='ID',suffixes=c("_Unique","_fracMM"),all.x=T)
colnames(FtrCount_merged)[which(colnames(FtrCount_merged)%in%c('Chr','Start','End','Strand'))]=c('chr','start','end','strand')
FtrCount_merged=FtrCount_merged[duplicated(FtrCount_merged)==F,]

##########################################################################################
############### splice junctions
##########################################################################################
if (join_junction==TRUE) {
  
  print("Running Join Junction")
  
  #read in jcounts file
  if (peak_type=="UNIQUE") {
    FtrCount_fracJCount = read.delim(paste0(peak_unique,".jcounts"),
                                     header=TRUE,sep="\t",stringsAsFactors = FALSE,comment.char = '#') 
  }  else {
    FtrCount_fracJCount = read.delim(paste0(peak_all,".jcounts"), 
                                   header=TRUE,sep="\t",stringsAsFactors = FALSE,comment.char = '#')
  } 
  
  #rename last col
  colnames(FtrCount_fracJCount)[ncol(FtrCount_fracJCount)]='counts'
  
  #clear na's
  FtrCount_fracJCount = FtrCount_fracJCount[is.na(FtrCount_fracJCount$PrimaryGene)==FALSE,]
  
  ###phil need to determine what to do if there are no splice junctions identified
  #if there are junctions identifitied merge
  if(nrow(FtrCount_fracJCount)>0){
    print("Junctions were identitified")
    FtrCount_fracJCount_merged = merge(FtrCount_fracJCount, 
                                       FtrCount_merged[,c("ID","strand")],by.x="PrimaryGene",by.y="ID",all.x=T)
    #split primary gene col
    FtrCount_fracJCount_merged = separate(FtrCount_fracJCount_merged,
                                          PrimaryGene,into=c('chr','start','end'),sep = ":|-",remove = F)
    
    ### remove junctions where splicing is within peak
    site_1 = FtrCount_fracJCount_merged$Site1_location >=FtrCount_fracJCount_merged$start&FtrCount_fracJCount_merged$Site1_location<=FtrCount_fracJCount_merged$end
    site_2 = FtrCount_fracJCount_merged$Site2_location>=FtrCount_fracJCount_merged$start&FtrCount_fracJCount_merged$Site2_location<=FtrCount_fracJCount_merged$end
    FtrCount_fracJCount_merged=FtrCount_fracJCount_merged[site_1!=site_2,]
    FtrCount_fracJCount_merged$JunctionID = paste0(FtrCount_fracJCount_merged$Site1_chr,':',
                                                   FtrCount_fracJCount_merged$Site1_location,'-', 
                                                   FtrCount_fracJCount_merged$Site2_location)
    
    #create Unique Row/junction ID
    FtrCount_fracJCount$rID=seq(1,nrow(FtrCount_fracJCount))
    
    
    jcount1.GR = GRanges(seqnames = as.character(FtrCount_fracJCount_merged$Site1_chr), 
                         ranges=IRanges(start = as.numeric(FtrCount_fracJCount_merged$Site1_location), 
                                        end = as.numeric(FtrCount_fracJCount_merged$Site1_location)),
                         strand = FtrCount_fracJCount_merged$strand,
                         PrimaryGene=FtrCount_fracJCount_merged$PrimaryGene,
                         count=FtrCount_fracJCount_merged$counts,
                         JunctionID=FtrCount_fracJCount_merged$JunctionID,
                         rID=FtrCount_fracJCount$rID)
    jcount2.GR = GRanges(seqnames = as.character(FtrCount_fracJCount_merged$Site2_chr),
                         ranges=IRanges(start = as.numeric(FtrCount_fracJCount_merged$Site2_location),
                                        end = as.numeric(FtrCount_fracJCount_merged$Site2_location)),
                         strand = FtrCount_fracJCount_merged$strand,
                         PrimaryGene=FtrCount_fracJCount_merged$PrimaryGene,
                         count=FtrCount_fracJCount_merged$counts,
                         JunctionID=FtrCount_fracJCount_merged$JunctionID,
                         rID=FtrCount_fracJCount$rID)
    
    FtrCount.GR=GRanges(seqnames = as.character(FtrCount_merged$chr), 
                        ranges=IRanges(start = as.numeric(FtrCount_merged$start), 
                                       end = as.numeric(FtrCount_merged$end)),
                        strand = FtrCount_merged$strand,
                        ID=FtrCount_merged$ID)
    
    #Add Junction ID
    junctionID <- function(Ftr_in, j_in, j_name, p_name){
      xo=as.data.frame(GenomicRanges::findOverlaps(j_in, Ftr_in, 
                                                   type = "any",ignore.strand=F))
      qh=as.data.frame(j_in[xo$queryHits],row.names = NULL)
      colnames(qh)=paste0(colnames(qh),j_name)

      sh=as.data.frame(Ftr_in[xo$subjectHits],row.names = NULL)
      colnames(sh)=paste0(colnames(sh),p_name)
      
      return(cbind(qh,sh))
    }
    
    ###phil - this doesn't work - the number of rows are different
    Junc_PL1 = junctionID(FtrCount.GR, jcount1.GR,"_Junction1","_Peaks1")
    Junc_PL2 = junctionID(FtrCount.GR, jcount2.GR,"_Junction2","_Peaks2")
    #Junc_PL_merged = cbind(Junc_PL1, Junc_PL2)
    Junc_PL=merge(Junc_PL1,Junc_PL2,by.x='rID_Junction1',by.y='rID_Junction2')
    Junc_PL=Junc_PL[(is.na(Junc_PL$ID_Peaks1)|is.na(Junc_PL$ID_Peaks2))==F,]
    
    # identify PEAK ID for Site1 and Site2 location of splice junction
    Junc_PLOut = Junc_PL_merged[,c('PrimaryGene_Junction1','JunctionID_Junction1',
                                   'strand_Peaks1','ID_Peaks1','ID_Peaks2')]
    
    # Identify peaks that are linked between splice junctions
    # Remove Peaks where PrimaryGene - Peak1, Peak2 are duplicated
    # Junction ID can be off by a few NT because different splice start/end in peak
    Junc_PLOut=Junc_PLOut[duplicated(Junc_PLOut[,c(1,4,5)])==F,]
    
    # Combine by primaryGene ID - Junction ID acts as ID to distinguish all possible junctions
    # get all unique PrimaryGene_Junction
    prim_junct=unique(Junc_PLOut$PrimaryGene_Junction1)
    PGene_TBL = setNames(data.frame(matrix(nrow=length(prim_junct), ncol=3)),
                           c('PrimaryGene_Junction1','JunctionID_Junction1','ID_comb'))
    
    for (junctions in 1:length(prim_junct)) {
      PrimGene=prim_junct[junctions]
      g=Junc_PL_merged[Junc_PL_merged$PrimaryGene_Junction1%in%PrimGene,]
      PGene_TBL[junctions,'PrimaryGene_Junction1']=unique(g$PrimaryGene_Junction1)
      PGene_TBL[junctions,'JunctionID_Junction1']=paste(unique(sort(g$JunctionID_Junction1)),collapse = ",")
      PGene_TBL[junctions,'ID_comb']=paste(unique(c(g$ID_Peaks1,g$ID_Peaks2)),collapse = ",")
    }
    #PGene_TBL=as.data.frame(PGene_TBL)
    #remove duplicated junctions
    PGene_TBL_unique=PGene_TBL[duplicated(PGene_TBL)==FALSE,]
    
    ##### Combine all connected peaks
    # get all unique PrimaryGene_Junction
    prim_junct = unique(PGene_TBL_unique$PrimaryGene_Junction1)
    PGene_TBL2=(matrix(nrow=nrow(PGene_TBL_unique),ncol=ncol(PGene_TBL_unique)+4))
    colnames(PGene_TBL2)=c(colnames(PGene_TBL_unique),'chr','start','end','linkedID')
    
    for (junction in 1:length(prim_junct)) {
      a = PGene_TBL[junction,]
      id = unlist(str_split (a$ID_comb,pattern = ","))
      idcomb = PGene_TBL[grep(paste(id,collapse = "|"),PGene_TBL$ID_comb),]
      
      ## look again to see if more locations with found connected peaks
      id2 = sort(unique(unlist(str_split (idcomb$ID_comb,pattern = ","))))
      
      ## If new peaks come up check until no new peaks  
      while (FALSE%in%(id2%in%id)) {
        id=id2
        
        # Look Again
        idcomb=PGene_TBL[grep(paste(id2,collapse = "|"),PGene_TBL$ID_comb),]
        
        # xxx3 Get Peak IDS to end while loop
        id2=sort(unique(unlist(str_split (idcomb$ID_comb,pattern = ","))))
      }
      
      id2_coord = separate(as.data.frame(id2),1,sep = ":|-",
                           into = c('chr','start','end'))
      id2_coord$start=as.numeric(id2_coord$star)
      id2_coord$end=as.numeric(id2_coord$end)
      PGene_TBL2[junction,'PrimaryGene_Junction1']=paste(sort(unique(idcomb$PrimaryGene_Junction1)),collapse = ',')
      PGene_TBL2[junction,'JunctionID_Junction1']=paste(sort(unique(idcomb$JunctionID_Junction1)),collapse = ',')
      PGene_TBL2[junction,'ID_comb']=paste(sort(unique(idcomb$ID_comb)),collapse = ',')
      PGene_TBL2[junction,'chr']=paste(unique(id2_coord$chr),collapse=",") ###phil added paste here
      PGene_TBL2[junction,'start']=min((id2_coord[,c('start','end')]))
      PGene_TBL2[junction,'end']=max((id2_coord[,c('start','end')]))
    }
    
    PGene_TBL2=as.data.frame(PGene_TBL2)
    PGene_TBL2=PGene_TBL2[duplicated(PGene_TBL2)==F,]
    PGene_TBL2$linkedID=paste0(PGene_TBL2$chr,":",PGene_TBL2$start,"-",PGene_TBL2$end)
    
    ##### check no reeated peaks 
    ##### get all unique PrimaryGene_Junction
    d3=unique(PGene_TBL2$PrimaryGene_Junction1)
    PGene_TBL3 = (matrix(nrow=nrow(PGene_TBL2),
                         ncol=ncol(PGene_TBL2)+1));colnames(PGene_TBL3)=c(colnames(PGene_TBL2),"nrows")
    
    for (x in 1:length(d3)) {
      a2=PGene_TBL2[x,]
      id2=unlist(str_split (a2$ID_comb,pattern = ","))
      idcomb2=PGene_TBL2[grep(paste(id2,collapse = "|"),PGene_TBL2$ID_comb),]
      
      PGene_TBL3[x,1:3]=as.matrix(a2[,1:3])
      PGene_TBL3[x,'nrows']=nrow(idcomb2)
    }
    PGene_TBL3=as.data.frame(PGene_TBL3)
    
    # ##### Get all peaks with junctions
    Junc_peaks=unique(c(Junc_PL_merged$ID_Peaks1,Junc_PL_merged$ID_Peaks2))
    FtrCount_Junc_peaks=FtrCount_merged[FtrCount_merged$ID%in%Junc_peaks,]
    FtrCount_merged_junction=FtrCount_merged[FtrCount_merged$ID%in%Junc_peaks==F,]
    
    ### Trim peaks without Splicing
    FtrCount_nosplice=FtrCount_merged_junction[FtrCount_merged_junction$Counts_fracMM>=read_depth,]
    FtrCount_splice_junc=rbind(FtrCount_nosplice,FtrCount_Junc_peaks)
    FtrCount_splice_junc=FtrCount_splice_junc[duplicated(FtrCount_splice_junc)==F,]
    FtrCount_splice_junc=FtrCount_splice_junc[!FtrCount_splice_junc$chr%in%c('chrM'),]
    
    if (DEmethod=='MANORM') {
      print("Running MANORM")
      FtrCount_splice_junc$Counts_fracMM >
      write.table(FtrCount_splice_junc[, c('chr','start','end','ID','ID','strand')],
                    file=paste0(out_dir, sample_id,"_", nt_merge,"_peakDepth.bed"), 
                    sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
      write.table(FtrCount_splice_junc[FtrCount_splice_junc$strand%in%"+", c('chr','start','end','ID','ID','strand')],
                  file=paste0(out_dir,sample_id,"_",nt_merge,"_peakDepth_P.bed"), 
                  sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
      write.table(FtrCount_splice_junc[FtrCount_splice_junc$strand%in%"-", c('chr','start','end','ID','ID','strand')],
                  file=paste0(out_dir,sample_id,"_",nt_merge,"_peakDepth_N.bed"), 
                  sep = "\t", row.names = FALSE, col.names = F, append = F, quote= FALSE,na = "")
    } else{
      print("MANORM skipped")
    } 
    
    #format final table
    colnames(FtrCount_splice_junc)[which(colnames(FtrCount_splice_junc) %in%
                                           c("V1","V2","V3","V6"))]=c('chr','start','end','strand')
  } else{
    print("No junctions were identified")
    ###phil if there are no junctoins, then just use the merge?
    FtrCount_splice_junc=FtrCount_merged[c("chr","start","end","strand")]
  }
}  else {
  print("Junctions joining not selected")
  FtrCount_splice_junc=FtrCount_merged[c("chr","start","end","strand")]
}

#write final junction output
write.csv(FtrCount_splice_junc,paste0(out_dir,sample_id,"_",nt_merge,"_peakjunction.txt"))
