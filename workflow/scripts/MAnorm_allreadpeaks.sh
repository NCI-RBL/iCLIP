#swarm -f MAnorm_allreadpeaks.sh --verbose 1 -g 50 -t 32 -b 3 --time 6:00:00 --gres=lscratch:200 --job-name PeaksMAnorm --module manorm/1.1.4

## one input should indicate comparison of interest
## lets make the input a table where column 1 is sample of interest and column2 is bacgroud eg contrasts.txt

mkdir outdir/15_MAnorm/Samplevsbackground

## subset Peaks file (sample + background) into positive peaks
## for unique read peaks use outdir/10_bed/Sample_Clip_unique.bed
awk '$6=="+"' outdir/10_bed/Sample_Clip_all.bed > outdir/10_bed/Sample_Clip_all.P.bed ## can be tempfile
awk '$6=="+"' outdir/10_bed/background_Clip_all.bed  > outdir/10_bed/background_Clip_all.P.bed ## can be tempfile

## 09_peak_annotation.R creates --p input files for each sample outdir/15_MAnorm/input/
manorm \
--p1 outdir/15_MAnorm/input/Sample_50nt_ALL_PeaksforMAnrom_P.bed \
--p2 outdir/15_MAnorm/input/background_50nt_ALL_PeaksforMAnrom_P.bed \
--r1 outdir/10_bed/Sample_Clip_all.P.bed \
--r2 outdir/10_bed/background_Clip_all.P.bed \
--s1 0 \
--s2 0 \
-p 1 \
-d 25 \
-n 10000 \
-s \
-o outdir/15_MAnorm/Samplevsbackground/ \
--name1 Sample_Pos \
--name2 background_Pos

###########################


awk '$6=="-"' outdir/10_bed/Sample_Clip_all.bed > outdir/10_bed/Sample_Clip_all.N.bed ## can be tempfile
awk '$6=="-"' outdir/10_bed/background_Clip_all.bed  > outdir/10_bed/background_Clip_all.N.bed ## can be tempfile

## 09_peak_annotation.R creates --p input files for each sample outdir/15_MAnorm/input/
manorm \
--p1 outdir/15_MAnorm/input/Sample_50nt_ALL_PeaksforMAnrom_N.bed \
--p2 outdir/15_MAnorm/input/background_50nt_ALL_PeaksforMAnrom_N.bed \
--r1 outdir/10_bed/Sample_Clip_all.N.bed \
--r2 outdir/10_bed/background_Clip_all.N.bed \
--s1 0 \
--s2 0 \
-p 1 \
-d 25 \
-n 10000 \
-s \
-o outdir/15_MAnorm/Samplevsbackground/ \
--name1 Sample_Neg \
--name2 background_Neg

