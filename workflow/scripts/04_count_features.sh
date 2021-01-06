num_of_threads=$1
anno_file=$2
bam_file=$3
out_dir=$4
sample_name=$5

#unique features
cmd="featureCounts \
$bam_file \
-T $num_of_threads \
-a $anno_file -F SAF \
-O -s 1 \
-R BAM --Rpath $out_dir \
-o "${out_dir}${sample_name}_FCount_unique.txt""
$cmd

#all features
cmd="featureCounts \
$bam_file \
-T $num_of_threads \
-a $anno_file -F SAF \
-s 1 -M -O --minOverlap 1 \
-R BAM --Rpath $out_dir \
-o "${out_dir}${sample_name}_FCount_all.txt""
$cmd

#all fraction
cmd="featureCounts \
$bam_file \
-T $num_of_threads \
-a $anno_file -F SAF \
-s 1 -M -O --minOverlap 1 --fraction \
-R BAM --Rpath $out_dir \
-o "${out_dir}${sample_name}_FCount_AllFrac.txt""
$cmd

#all primary
cmd="featureCounts \
$bam_file \
-T $num_of_threads \
-a $anno_file -F SAF \
-s 1 -M -O --minOverlap 1 --primary \
-R BAM --Rpath $out_dir \
-o "${out_dir}${sample_name}_FCount_AllPrim.txt""
$cmd

#all fraction primary
cmd="featureCounts \
$bam_file \
-T $num_of_threads \
-a $anno_file -F SAF \
-s 1 -M -O --minOverlap 1 --fraction --primary \
-R BAM --Rpath $out_dir \
-o "${out_dir}${sample_name}_FCount_AllFracPrim.txt""
$cmd
