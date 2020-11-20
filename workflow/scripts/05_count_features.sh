num_of_threads=$1
anno_file=$2
bam_file=$3
sample_name=$4
output_dir=$5

#options
#-O = overlaping metafeatures
#-s 1 = stranded
#-T = number of threads

featureCounts \
-a $anno_file \
-o "${sample_name}_FCountUnique.txt"
$bam_file \
-F SAF \
-O \
-s 1 \
-T $num_of_threads \
-R BAM --Rpath $out_dir \

featureCounts \
-a $anno_file \
-o "${sample_name}_FCountAll.txt"
$bam_file \
-F SAF \
-s 1 \
-M \
-O \
--minOverlap 1 \
-T $num_of_threads \
-R BAM --Rpath $out_dir \

featureCounts \
-a $anno_file \
-o "${sample_name}_FCountAll_frac.txt"
$bam_file \
-F SAF \
-s 1 \
-M \
-O \
--fraction \
--minOverlap 1 \
-T $num_of_threads \
-R BAM --Rpath $out_dir \

featureCounts \
-a $anno_file \
-o "${sample_name}_FCountAll_primary.txt"
$bam_file \
-F SAF \
-s 1 \
-M \
-O \
--primary \
--minOverlap 1 \
-T $num_of_threads \
-R BAM --Rpath $out_dir \

featureCounts \
-a $anno_file \
-o "${sample_name}_FCountAll_FracPrime.txt"
$bam_file \
-F SAF \
-s 1 \
-M \
-O \
--primary \
--fraction \
--minOverlap 1 \
-T $num_of_threads \
-R BAM --Rpath $out_dir \
