### Local
# Rscript CLIP_Out_anno.R \
# WT \
# hg38 \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.iPeaks50nt.FCountUnique.txt" \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.iPeaks50nt.FCountall_frac.txt" \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/" \
# 5 \
# all \
# MAnorm \
# Test #_SpliceAware_FullGenome_Transcme.SplicTrans_CollapseExons
# 
# Rscript CLIP_Out_anno.R \
# KO \
# hg38 \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/KO/all/KO_CLIP.dedup.iPeaks50nt.FCountUnique.txt" \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/KO/all/KO_CLIP.dedup.iPeaks50nt.FCountall_frac.txt" \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/" \
# 5 \
# all \
# MAnorm \
# Test #_SpliceAware_FullGenome_Transcme.SplicTrans_CollapseExons



### Biowulf
# module load R
# module load bedtools
# module load samtools

# Rscript CLIP_Out_anno.R \
# WT \
# human \ # change to hg38 or mm10
# /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.iPeaks50nt.FCountUnique.txt \
# /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.iPeaks50nt.FCountall_frac.txt \
# /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/ \
# 5 \
# all \
# MAnorm \
# _SpliceAware_FullGenome_Transcme.SplicTrans_CollapseExons

###########
args <- commandArgs(trailingOnly = TRUE)
samplename= as.character(args[1]) # currently WT or KO
species=as.character(args[2]) # human or mouse or mm10 or hg38
PeaksUniq=as.character(args[3]) # path to peaks counts of unique reads
PeaksFracMM=as.character(args[4]) # path to peaks counts of Unique + fracMM reads
OutDIR=as.character(args[5]) # path to place Annotation table 
readdepth=as.numeric(args[6]) ## minimum number of reads to evaluate a peak
PeakIdnt=as.character(args[7]) ## Unique or all
DEmethod=as.character(args[8]) ## DE peaks method for now only MAnorm works or none
#NameAdd=as.character(args[9]) ## additonal information to add to file name
ntmerge=as.character(args[9]) #nt merge length
outName=as.character(args[10]) #output file name
spliceaware=as.character(args[11])
OutDIR2=as.character(args[12]) # path to place Annotation table 

PeaksUniq=gsub("\\\\ "," ",PeaksUniq)
PeaksFracMM=gsub("\\\\ "," ",PeaksFracMM)
OutDIR=gsub("\\\\ "," ",OutDIR)

PeaksUniq

Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

#File output : params$OutDIR,samplename,"_CLIP_",params$PeakIdnt,"Peaks",ntmerge,'_peakDepth',params$readdepth,params$NameAdd,"_MD15.txt"
rmarkdown::render('workflow/scripts/CLIP_MD15_BW.Rmd',
                  encoding=encoding,
                  output_file=outName,
                  params = list(
                    samplename= samplename,
                    # INbam = paste0(bamDIR,'ddup/'),
                    # INpeaks= paste0(PeaksDIR),
                    PeaksUniq=PeaksUniq,
                    PeaksFracMM=PeaksFracMM,
                    # inMAPQ: '/MAPQ/',
                    OutDIR=OutDIR,
                    OUT= paste0(OutDIR,"/annotation/"),
                    readdepth= readdepth,
                    PeakIdnt=PeakIdnt,
                    #NameAdd=NameAdd,
                    DEmethod=DEmethod,
                    species= species,
                    ntmerge = ntmerge,
                    JoinJunc = spliceaware,
                    OutDIR2 = OutDIR2
                  )
)
