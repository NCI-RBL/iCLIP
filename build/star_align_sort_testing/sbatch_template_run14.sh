#!/bin/sh
module load samtools

# set tmp dir
tmp_dir="/lscratch/${SLURM_JOB_ID}"
export tmp_dir

# set other variables
sample_id="FILL_SAMPLE"
run_id="FILL_RUN"
output_dir="/data/sevillas2/marco_star/$run_id"
output_bam="$output_dir/$run_id/$sample_id.bam"
output_log="$output_dir/logs/${run_id}_${sample_id}.log"
output_unmapped="$output_dir/$run_id/${sample_id}_unmapped.bam"

# Sort and Index
#samtools sort -m 80G -T $tmp_dir $output_dir/tmp/${sample_id}_Aligned.out.bam -o $output_dir/${sample_id}_Aligned.sortedByCoord.out.bam
samtools index -@ 32 $output_dir/${sample_id}_Aligned.sortedByCoord.out.bam
