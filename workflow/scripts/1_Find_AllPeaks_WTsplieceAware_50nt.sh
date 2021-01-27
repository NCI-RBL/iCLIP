#swarm -f 1_Find_AllPeaks_WTsplieceAware_50nt.sh --verbose 1 -g 25 -t 32 -b 25 --time 6:00:00 --gres=lscratch:200 --job-name allpeaks_50nt_WT_recomb --module samtools,R,bedtools,subread
  

###############################
### Peak counting step
###############################


bedtools bamtobed -split -i /data/RBL_NCI/Wolin/CLIP_Pipeline/iCLIP/fCLIP/fCLIP_HS/10_dedup_bam/WT_fCLIP.dedup.i.bam | bedtools sort -i - > /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/WT_CLIP.dedup.i.bam.bed

bedtools merge -c 6 -o count,distinct -bed -s -d 50 -i /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/WT_CLIP.dedup.i.bam.bed | \
awk '{OFS="\t"; print $1":"$2"-"$3,$1,$2,$3,$5}'| awk 'BEGIN{print "ID","Chr","Start","End","Strand"}1' > /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.i.bam.allpeaks50nt.bedtools.SAF



#### Feature Counts to Count reads under peak

#### Unique reads (fractional counts correctly count splice reads for each peak. When peaks counts are combined for peaks connected by splicing in Rscript)
featureCounts -F SAF \
-a /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.i.bam.allpeaks50nt.bedtools.SAF \
-O \
-J \
--fraction \
--minOverlap 1 \
-s 1 \
-T 8 \
-o /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.iPeaks50nt.FCountUnique.txt /data/RBL_NCI/Wolin/CLIP_Pipeline/iCLIP/fCLIP/fCLIP_HS/10_dedup_bam/WT_fCLIP.dedup.i.bam



#### Include Multimap reads - MM reads given fractional count based on # of mapping locations. All spliced reads also get fractional count. So Unique reads can get fractional count when spliced peaks combined in R script the summed counts give whole count for the unique alignement in combined peak.

featureCounts -F SAF \
-a /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.i.bam.allpeaks50nt.bedtools.SAF \
-M \
-O \
-J \
--fraction \
--minOverlap 1 \
-s 1 \
-T 8 \
-o /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/all/WT_CLIP.dedup.iPeaks50nt.FCountall_frac.txt /data/RBL_NCI/Wolin/CLIP_Pipeline/iCLIP/fCLIP/fCLIP_HS/10_dedup_bam/WT_fCLIP.dedup.i.bam

#### Count Unique + MM reads MM reads only get counted once based on primary alignment of read - don't think this is necisary
#  featureCounts -F SAF \
-a /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/All/WT_CLIP.dedup.i.bam.Allpeaks50nt.bedtools.SAF \
-s 1 \
-M \
-O \
-J \
--primary \
--fraction \
--minOverlap 1 \
-T 8 \
-R BAM \
--Rpath /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT \
-o /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/All/WT_CLIP.dedup.iPeaks50nt.FCountAll_FracPrime.txt /data/RBL_NCI/Wolin/CLIP_Pipeline/iCLIP/CLIP/CLIP_HS/10_dedup_bam/WT_CLIP.dedup.i.bam
 

#### Count Unique + MM reads - must include fractional count to correctly count splicing, don't think this is necisary
# featureCounts -F SAF \
-a /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/All/WT_CLIP.dedup.i.bam.Allpeaks50nt.bedtools.SAF \
-s 1 \
-M \
-O \
-J \
--primary \
--minOverlap 1 \
-T 8 \
-R BAM \
--Rpath /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT \
-o /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/All/WT_CLIP.dedup.iPeaks50nt.FCountAll_primary.txt /data/RBL_NCI/Wolin/CLIP_Pipeline/iCLIP/CLIP/CLIP_HS/10_dedup_bam/WT_CLIP.dedup.i.bam





