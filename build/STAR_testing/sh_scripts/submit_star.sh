#!/bin/bash
module load STAR
sample_id="fill_sample"
project_id="fill_project"
seedPerWindowNmax="fill_seed"
winAnchorMultimapNmax="fill_anchor"
output_dir="/data/sevillas2/star/${sample_id}/fill_trial"

input_fq="/data/RBL_NCI/Wolin/$project_id/01_preprocess/01_fastq/${sample_id}.fastq.gz"; 
sample_id_prefix="${sample_id}_"

tmp_dir="/lscratch/${SLURM_JOB_ID}"
export tmp_dir

STAR --runMode alignReads --genomeDir /data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/mm10/index \
--sjdbGTFfile /data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/mm10/ref/gencode.vM23.annotation.gtf --readFilesCommand zcat \
--readFilesIn $input_fq \
--outFileNamePrefix $tmp_dir/${sample_id_prefix} \
--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --alignEndsType Local --alignIntronMax 50000 --alignSJDBoverhangMin 3 \
--alignSJoverhangMin 5 --alignTranscriptsPerReadNmax 10000 --alignWindowsPerReadNmax 10000 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.9 \
--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 10000 --outFilterMultimapScoreRange 0 --outFilterScoreMin 0 \
--outFilterType Normal --outSAMattributes All --outSAMunmapped None --outSJfilterCountTotalMin 3 1 1 1 --outSJfilterOverhangMin 30 12 12 12 \
--outSJfilterReads All --seedMultimapNmax 10000 --seedNoneLociPerWindow 20 --seedPerReadNmax 10000 --seedPerWindowNmax $seedPerWindowNmax --sjdbScore 2 --winAnchorMultimapNmax $winAnchorMultimapNmax

# move STAR files and final log file to output
mv $tmp_dir/${sample_id_prefix}Aligned.sortedByCoord.out.bam ${output_dir}/${sample_id}.bam;
mv $tmp_dir/${sample_id_prefix}Log.final.out ${output_dir}/${sample_id}.out; 
    
# move mates to unmapped file
touch ${output_dir}/${sample_id}.unmapped.out; 
for f in $tmp_dir/${sample_id_prefix}Unmapped.out.mate*; do cat $f >> ${output_dir}/${sample_id}.unmapped.out; done