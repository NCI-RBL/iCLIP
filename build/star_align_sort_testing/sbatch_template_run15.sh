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
samtools sort --threads 8 -m 100G -T $tmp_dir $output_dir/tmp/${sample_id}_Aligned.out.bam -o $output_dir/${sample_id}_Aligned.sortedByCoord.out.bam
samtools index -@ 32 $tmp_dir/${sample_id}_Aligned.sortedByCoord.out.bam
#samtools sort --threads 32 -m 40G -T $tmp_dir $output_dir/tmp/${sample_id}_Aligned.out.bam -o $output_dir/${sample_id}_Aligned.sortedByCoord.out.bam
#samtools index -@ 32 $tmp_dir/${sample_id}_Aligned.sortedByCoord.out.bam

# move STAR files and final log file to output
mv $tmp_dir/${sample_id}_Aligned.sortedByCoord.out.bam $output_bam
mv $tmp_dir/${sample_id}_Log.final.out $output_log
        
# move mates to unmapped file
touch $output_unmapped
for f in $tmp_dir/${sample_id}_Unmapped.out.mate*; do cat $f >> $output_unmapped; done