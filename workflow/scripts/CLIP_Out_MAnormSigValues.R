### Local
# Rscript CLIP_Out_MAnormSigValues.R \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/" \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/peaks_MAnorm/" \
# "/Users/homanpj/OneDrive\ -\ National\ Institutes\ of\ Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/contrasts.txt" \
# "/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/peaks_MAnorm/" \
# 5 \
# all \
# Test #_SpliceAware_FullGenome_Transcme.SplicTrans_CollapseExons

# contrasts=as.character('/Users/homanpj/OneDrive - National Institutes of Health/Loaner/Wolin/CLIP/HaCat_fCLIP2/contrasts.txt')

### Biowulf
# module load R
# module load samtools
# module load bedtools



###########
args <- commandArgs(trailingOnly = TRUE)
INpeaksAnno= as.character(args[1]) # path to folder with tables of annotated peaks
INMAnorm=as.character(args[2]) # path to MAnorm output folder with single comparison
contrasts=(args[3]) # contrasts table
outDIR=as.character(args[4]) # path to place MAnorm report and table
readdepth=as.character(args[5]) # minimum number of reads to evaluate a peak 
PeakIdnt=as.character(args[6]) ## Unique or all
NameAdd=as.character(args[7]) ## append name to output files
ntmerge='50nt'

library(data.table)
INpeaksAnno=gsub("\\\\ "," ",INpeaksAnno)
INMAnorm=gsub("\\\\ "," ",INMAnorm)
outDIR=gsub("\\\\ "," ",outDIR)
contrasts=gsub("\\\\ "," ",contrasts)
contrastF2=read.table(file=contrasts, header=F, sep="\t",stringsAsFactors = F,quote = "")
samplename=contrastF2[1,1]
background=contrastF2[1,2]

Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

rmarkdown::render('CLIP_MD15_SigPeaks_MAnorm_BW.Rmd',
                  encoding=encoding,
                  output_file=paste0(outDIR,samplename,"vs",background,"_CLIP_MAnorm_",PeakIdnt,"Peaks",ntmerge,'_peakDepth',readdepth,NameAdd,"_MD15.html"), 
                  params = list(
                    INpeaksAnno= INpeaksAnno,
                    INMAnorm=INMAnorm,
                    contrasts=contrasts,
                    outDIR=outDIR,
                    readdepth= readdepth,
                    PeakIdnt=PeakIdnt,
                    NameAdd=NameAdd,
                    ntmerge= ntmerge
                  )
)
