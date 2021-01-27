#swarm -f MAnorm.sh --verbose 1 -g 50 -t 32 -b 3 --time 6:00:00 --gres=lscratch:200 --job-name PeaksMAnorm --module manorm/1.1.4
#mkdir /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/WTvsKO


awk '$6=="+"' /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/WT_CLIP.dedup.i.bam.bed >                   /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/WT_CLIP.dedup.i.bam.P.bed

awk '$6=="+"' /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/KO/KO_CLIP.dedup.i.bam.bed > /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/KO_CLIP.dedup.i.bam.P.bed

manorm \
--p1 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/WT_allPeaksforMAnrom_50nt_peakDepth5Test.P.bed \
--p2 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/KO_allPeaksforMAnrom_50nt_peakDepth5Test.P.bed \
--r1 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/WT_CLIP.dedup.i.bam.P.bed \
--r2 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/KO_CLIP.dedup.i.bam.P.bed \
--s1 0 \
--s2 0 \
-p 1 \
-d 25 \
-n 10000 \
-s \
-o /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/WTvsKO/ \
--name1 WT_50nt_peakDepth5Test_Pos \
--name2 KO_50nt_peakDepth5Test_Pos

###########################


awk '$6=="-"' /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/WT/WT_CLIP.dedup.i.bam.bed >                   /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/WT_CLIP.dedup.i.bam.N.bed


awk '$6=="-"' /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_FullGenome_Transcme.SplicTrans2/KO/KO_CLIP.dedup.i.bam.bed > /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/KO_CLIP.dedup.i.bam.N.bed

manorm \
--p1 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/WT_allPeaksforMAnrom_50nt_peakDepth5Test.N.bed \
--p2 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/KO_allPeaksforMAnrom_50nt_peakDepth5Test.N.bed \
--r1 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/WT_CLIP.dedup.i.bam.N.bed \
--r2 /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/input/KO_CLIP.dedup.i.bam.N.bed \
--s1 0 \
--s2 0 \
-p 1 \
-d 25 \
-n 10000 \
-s \
-o /data/RBL_NCI/Wolin/Phil/HaCat_fCLIP2/peaks_MAnorm/WTvsKO \
--name1 WT_50nt_peakDepth5Test_Neg \
--name2 KO_50nt_peakDepth5Test_Neg


