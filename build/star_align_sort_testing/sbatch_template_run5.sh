#!/bin/sh
module load STAR

# set tmp dir
tmp_dir="/lscratch/${SLURM_JOB_ID}"
export tmp_dir

# set other variables
sample_id="FILL_SAMPLE"
run_id="FILL_RUN"
output_dir="/data/sevillas2/marco_star"

# check dirs exits
if [[ ! -d $output_dir/$run_id ]]; then mkdir -p $output_dir/$run_id; fi

# output files
output_bam="$output_dir/$run_id/$sample_id.bam"
output_log="$output_dir/logs/${run_id}_${sample_id}.log"
output_unmapped="$output_dir/$run_id/${sample_id}_unmapped.bam"

# run star
STAR --runMode alignReads --genomeDir /data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/hg38/index \
--sjdbGTFfile /data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/hg38/ref/gencode.v32.annotation.gtf \
--readFilesCommand zcat \
--readFilesIn /data/RBL_NCI/Wolin/8-09-21-HaCaT_fCLIP_v2.0_2/01_preprocess/01_fastq/$sample_id.fastq.gz \
--outFileNamePrefix $tmp_dir/${sample_id}_ \
--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --alignEndsType Local --alignIntronMax 50000 \
--alignSJDBoverhangMin 3 --alignSJoverhangMin 5 --alignTranscriptsPerReadNmax 10000 --alignWindowsPerReadNmax 10000 \
--outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.9 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
--outFilterMultimapNmax 10000 --outFilterMultimapScoreRange 0 --outFilterScoreMin 0 --outFilterType Normal --outSAMattributes All \
--outSAMunmapped None --outSJfilterCountTotalMin 3 1 1 1 --outSJfilterOverhangMin 30 12 12 12 --outSJfilterReads All --seedMultimapNmax 10000 \
--seedNoneLociPerWindow 20 --seedPerReadNmax 10000 --seedPerWindowNmax 500 --sjdbScore 2 --winAnchorMultimapNmax 500 --outBAMsortingBinsN 600

# move STAR files and final log file to output
mv $tmp_dir/${sample_id}_Aligned.sortedByCoord.out.bam $output_bam
mv $tmp_dir/${sample_id}_Log.final.out $output_log
        
# move mates to unmapped file
touch $output_unmapped
for f in $tmp_dir/${sample_id}_Unmapped.out.mate*; do cat $f >> $output_unmapped; done