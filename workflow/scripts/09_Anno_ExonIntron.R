#library
suppressMessages(library(tidyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
suppressMessages(library(parallel))



#set args
parser <- ArgumentParser()
parser$add_argument("-r","--rscript", dest="rscript", required=TRUE, help="rscript for functions")
parser$add_argument("-pt","--peak_type", dest="peak_type", required=TRUE, help="peak type (all, unique)")
# parser$add_argument("-pkt","--peak_table", dest="peak_table", required=TRUE, help="Processed Peak Table")
# parser$add_argument("-pu","--peak_unique", dest="peak_unique", required=TRUE, help="path for unique count peaks")
# parser$add_argument("-pa","--peak_all", dest="peak_all", required=TRUE, help="path MM count peaks")
# parser$add_argument("-ptot","--peak_total", dest="peak_total", required=TRUE, help="path total count peaks")
# parser$add_argument("-jj","--join_junction", dest="join_junction", required=TRUE, help="include junctions (T or F)")
parser$add_argument("-a","--anno_anchor", dest="anno_anchor", required=TRUE, help="5prime or max. Annotation for Spliced peaks will be based on 5` most region or region with most reads")
# parser$add_argument("-ce","--condense_exon", dest="condense_exon", required=TRUE, help="condense exons (T or F)")
parser$add_argument("-rdepth","--read_depth", dest="read_depth", required=TRUE, help="read depth filtering parameter")
# parser$add_argument("-de","--demethod", dest="demethod", required=TRUE, help="DE method (MAnorm, DiffBind or none)")
parser$add_argument("-s","--sample_id", dest="sample_id", required=TRUE, help="sample id")
parser$add_argument("-sp","--ref_species", dest="ref_species", required=TRUE, help="reference species")
parser$add_argument("-anno","--anno_dir", dest="anno_dir", required=TRUE, help="path for annotation dir")
# parser$add_argument("-st","--splice_table", dest="splice_table", required=TRUE, help="path for Splice Table")
parser$add_argument("-reft","--reftable_path", dest="reftable_path", required=TRUE, help="path for reftable")
parser$add_argument("-g","--gencode_path", dest="gencode_path", required=TRUE, help="path for gencode")
parser$add_argument("-i","--intron_path", dest="intron_path", required=TRUE, help="path for intron")
parser$add_argument("-rmsk","--rmsk_path", dest="rmsk_path", required=TRUE, help="path for rmsk")
parser$add_argument("-tmp","--tmp_dir", dest="tmp_dir", required=TRUE, help="path for tmp dir")
parser$add_argument("-od","--out_dir", dest="out_dir", required=TRUE, help="path for output dir")
parser$add_argument("-of","--out_file", dest="out_file", required=TRUE, help="path for output dir")
# parser$add_argument("-odm","--out_dir_DEP", dest="out_dir_DEP", required=FALSE, help="path for manorm output dir")
# parser$add_argument("-o","--output_file_error", dest="output_file_error", required=FALSE, help="path for output error file")
parser$add_argument("-as","--anno_strand",dest="anno_strand",required=T,help="Annotion for regions on same strand or opposite strand")


args <- parser$parse_args()
rscript = args$rscript
peak_type = args$peak_type
# peak_unique = args$peak_unique
# peak_all = args$peak_all
# peak_total = args$peak_total
# join_junction = as.logical(args$join_junction)
# condense_exon = as.logical(args$condense_exon)
anno_anchor=args$anno_anchor
read_depth = as.numeric(args$read_depth)
# DEmethod = args$demethod
sample_id = args$sample_id
ref_species = args$ref_species
anno_dir = args$anno_dir
splice_table=args$splice_table
reftable_path = args$reftable_path
gencode_path = args$gencode_path
intron_path = args$intron_path
rmsk_path = args$rmsk_path
tmp_dir = args$tmp_dir
out_dir = args$out_dir
out_file = args$out_file
# out_dir_DEP = args$out_dir_DEP
# output_file_error = args$output_file_error
anno_strand = args$anno_strand
# subset=as.numeric(args$subset)


test=F
if(test){
  rm(list=setdiff(ls(), "params"))
  
  '/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput2/04_annotation/02_peaks/FLAG_Ro_noPFA_ALLreadPeaks_AllRegions.txt'

  wd="/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/testing/testoutput2/"
  setwd(wd)
  wd="."
  
  
  sample_id = "FLAG_Ro_noPFA"
  
  peak_type= "ALL"
  anno_dir = paste0(wd,"/04_annotation/01_project/")
    ref_path = "/Users/homanpj/Documents/Resources/ref/"
    tmp_dir= paste0(wd,"/04_annotation/tmp/")
    
    reftable_path = "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/CLIPpipeline/iCLIP/config/annotation_config.txt"
  
    #output
    out_dir = paste0(wd,"/04_annotation/02_peaks/")
    out_dir_DEP =paste0(wd,"/05_MAnorm/")
    
    #project annotation files
    anno_dir = paste0(wd,"/04_annotation/01_project/")
    ref_path = "/Users/homanpj/Documents/Resources/ref/"
    tmp_dir= paste0(wd,"/04_annotation/02_peaks/")
    
    
    #feature information
    join_junction = "TRUE"
    anno_anchor="max_total"
    read_depth = 3
    ref_species="mm10"
    anno_strand="SameStrand"
    
    rscript="/Volumes/RBL_NCI/Wolin/Phil/test/testoutput2/09_peak_annotation_functions.R"
    
    if(ref_species == "mm10"){
      gencode_path = paste0(ref_path, "mm10/Gencode_VM23/fromGencode/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.txt")
      intron_path = paste0(ref_path, "mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed")
      rmsk_path = paste0(ref_path,"mm10/repeatmasker/rmsk_GRCm38.txt")
      
    } else if (ref_species == "hg38"){
      gencode_path = paste0(ref_path,"hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt")
      intron_path = paste0(ref_path,"hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed")
      rmsk_path = paste0(ref_path,"hg38/repeatmasker/rmsk_GRCh38.txt")
    }
  }

##########################################################################################
#source R script with functions
source(rscript)

#annotation paths
gencode_transc_path = paste0(anno_dir,"ref_gencode.txt")
lncra_path = paste0(anno_dir,"lncRNA_gencode.txt")
YRNA_path = paste0(anno_dir, "yRNA_repeatmasker.bed")
srpRNA_path = paste0(anno_dir, "srpRNA_repeatmasker.bed")
SKRNA_path = paste0(anno_dir, "7SKRNA_repeatmasker.bed")
scRNA_path = paste0(anno_dir, "scRNA_repeatmasker.bed")
tRNA_path = paste0(anno_dir, "tRNA_repeatmasker.bed")
sncRNA_path = paste0(anno_dir, "sncRNA_gencode.bed")
rRNA_BK00964_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
rRNA_rmsk_path = paste0(anno_dir, "rRNA_repeatmasker.bed")
tRNA_rmsk_path = paste0(anno_dir, "tRNA_repeatmasker.bed")

#set id for files
file_id = paste0(sample_id,"_")
# file_id = sample_id

##########################################################################################
## read file
##########################################################################################


FtrCount_trimmed=read.delim(paste0(out_dir,file_id,peak_type,'readPeaks_AllRegions.txt'), 
                               header=T,sep="\t",stringsAsFactors = F,comment.char = '#') 

#prep peaks
peaks_Oppo=FtrCount_trimmed[,c('chr','start','end','ID','ID','strand')] 
peaks_Oppo$strand=gsub("\\+","pos",peaks_Oppo$strand) 
peaks_Oppo$strand=gsub("\\-","+",peaks_Oppo$strand)
peaks_Oppo$strand=gsub("pos","-",peaks_Oppo$strand)

##########################################################################################
############### INTRON EXON ANNOTATION
##########################################################################################
### Identify if CLIP peak overlaps with Intron or Exonic region   

#   Using GTF file from GENCODE v23 -mm10
#   Using GTF file from GENCODE v32 -hg38

# Peaks were annotated by whether they overlap with Host gene intron/exon region
# Intron coordinates were calculated from GTF file.
# 
# A second column was added to idenify if the peak also overlapped with the 5'UTR 3'UTR or CDS (Column: Featrue 2)
#########################################################################################3
print("Intron/exon annotation")

intron_exon=IntronExon_prep(gencode_path,intron_path,gencode_transc_path)%>%suppressMessages()

##########################################################################################
############### Intron Exon Calling
##########################################################################################

print("IE Calling")


#run functions
if (anno_strand=="SameStrand") {
  
system.time({ print('Same strand')
  peaks_SameAnno = IE_calling(FtrCount_trimmed[,c('chr','start','end','ID','ID','strand')],'Same',"Same_")
  write.table(peaks_SameAnno,paste0(out_file),sep = "\t",row.names = F)
  
  
})

  }else if (anno_strand=="OppoStrand"){

system.time({ print('Oppo strand')
  peaks_OppoAnno = IE_calling(peaks_Oppo,'Oppo',"Oppo_")
  write.table(peaks_OppoAnno,paste0(out_file),sep = "\t",row.names = F)
  
})

}







