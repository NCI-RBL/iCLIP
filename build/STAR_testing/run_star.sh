#!/bin/bash

################################
#purpose
################################
# Benchmark STAR alignment parameters using iCLIP methods paper (Busch) as reference.
# Vary these initial parameters and compare to default, and ENCODE parameters.

#References
#https://www.sciencedirect.com/science/article/pii/S1046202318304948
#https://www.encodeproject.org/documents/739ca190-8d43-4a68-90ce-1a0ddfffc6fd/@@download/attachment/eCLIP_analysisSOP_v2.2.pdf
#https://github.com/alexdobin/STAR/releases

# Variables
#1 unique reads only
## a) directly from Busch paper
## b) increase mismatch threshold
## c) decrease mismatch threshold

#2 multimapped reads included
## a) directly from Busch paper
## b) increase mismatch threshold
## c) decrease mismatch threshold

#set dirs, files
parent_dir="/data/RBL_NCI/Wolin/mES_fclip_1_YL_012122/alignment_analysis/star"
input_dir=$parent_dir/input
output_dir=$parent_dir/output
ref_dir=$parent_dir/ref
index_dir=$ref_dir/index
ref_gen=$ref_dir/GRCm38.p6.genome.fa
ref_bk=$ref_dir/BK000964.3_TPA_exp.fa
ref_gtf=$ref_dir/gencode.vM23.annotation.gtf
log_dir=$parent_dir/log
variable_file=$parent_dir/docs/STAR_variations.txt

novo_list=("FLAG_Ro_fclip.dedup.si.Sptan1.unique.parameter_original.useq.bam" "FLAG_Ro_fclip.dedup.si.Sptan1.unique.parameter_set1.useq.bam" "FLAG_Ro_fclip.dedup.si.Sptan1.unique.parameter_set2.useq.bam" "FLAG_Ro_fclip.dedup.si.Sptan1.unique.parameter_set3.useq.bam")

#controls
flag_download="N"
flag_index="N"
flag_variables="N"

flag_align_partial="N"
flag_gap_partial="N"

flag_align_complete="N"
flag_gap_complete="N"
flag_align_stats="N"

flag_subset_rnu="N"
flag_align_rnu="N"
flag_align_stats_rnu="N"

flag_subset_bam_multiple_genes="N"
flag_subset_fq_multiple_genes="N"
flag_align_overhang="N"
flag_gap_overhang="Y"
flag_align_stats_overhang="Y"

#load
module load STAR/2.7.8a python samtools R

#functions
function star_index()
{
    local star_command="STAR \
        --runThreadN 32 \
        --runMode genomeGenerate \
        --genomeDir $index_dir \
        --genomeFastaFiles $ref_gen $ref_bk"
    echo "$star_command"
}

function star_run()
{
    if [[ ! -d $output_dir/$round_of_test/$version_id ]]; then mkdir -p $output_dir/$round_of_test/$version_id; fi

    local star_command="STAR \
        --runMode alignReads \
        --genomeDir $index_dir \
        --sjdbGTFfile $ref_gtf \
        --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax \
        --outFilterMismatchNmax $outFilterMismatchNmax \
        --outFilterMultimapNmax $outFilterMultimapNmax \
        --outSJfilterReads $outSJfilterReads \
        --outSJfilterCountTotalMin $outSJfilterCountTotalMin \
        --alignEndsType $alignEndsType \
        --outSAMunmapped $outSAMunmapped \
        --outSAMattributes $outSAMattributes \
        --outFilterType $outFilterType \
        --outFilterScoreMin $outFilterScoreMin \
        --outFilterMultimapScoreRange $outFilterMultimapScoreRange \
        --sjdbScore $sjdbScore \
        --outSJfilterOverhangMin $outSJfilterOverhangMin \
        --alignIntronMax $alignIntronMax \
        --outFilterMatchNmin $outFilterMatchNmin \
        --outFilterMatchNminOverLread $outFilterMatchNminOverLread \
        --alignSJDBoverhangMin $alignSJDBoverhangMin \
        --alignSJoverhangMin $alignSJoverhangMin \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn $input_dir/$input_fq \
        --outFileNamePrefix $output_dir/$round_of_test/$version_id/${version_id}_"
    echo "$star_command"
}

function star_run_expanded()
{
    if [[ ! -d $output_dir/$round_of_test/$version_id ]]; then mkdir -p $output_dir/$round_of_test/$version_id; fi

    local star_command="STAR \
        --runMode alignReads \
        --genomeDir $index_dir \
        --sjdbGTFfile $ref_gtf \
        --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax \
        --outFilterMismatchNmax $outFilterMismatchNmax \
        --outFilterMultimapNmax $outFilterMultimapNmax \
        --outSJfilterReads $outSJfilterReads \
        --outSJfilterCountTotalMin $outSJfilterCountTotalMin \
        --alignEndsType $alignEndsType \
        --outSAMunmapped $outSAMunmapped \
        --outSAMattributes $outSAMattributes \
        --outFilterType $outFilterType \
        --outFilterScoreMin $outFilterScoreMin \
        --outFilterMultimapScoreRange $outFilterMultimapScoreRange \
        --sjdbScore $sjdbScore \
        --outSJfilterOverhangMin $outSJfilterOverhangMin \
        --alignIntronMax $alignIntronMax \
        --outFilterMatchNmin $outFilterMatchNmin \
        --outFilterMatchNminOverLread $outFilterMatchNminOverLread \
        --alignSJDBoverhangMin $alignSJDBoverhangMin \
        --alignSJoverhangMin $alignSJoverhangMin \
        --winAnchorMultimapNmax=$winAnchorMultimapNmax \
        --seedPerReadNmax=$seedPerReadNmax \
        --seedPerWindowNmax=$seedPerWindowNmax \
        --alignWindowsPerReadNmax=$alignWindowsPerReadNmax \
        --alignTranscriptsPerReadNmax=$alignTranscriptsPerReadNmax \
        --seedMultimapNmax=$seedMultimapNmax \
        --seedNoneLociPerWindow=$seedNoneLociPerWindow \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn $input_dir/$input_fq \
        --outFileNamePrefix $output_dir/$round_of_test/$version_id/${version_id}_"
    echo "$star_command"
}

#analysis
if [[ "$flag_download" == "Y" ]]; then
    echo "*** Running download ***"

    #download / copy ref files
    ## mm10 v23 fa
    if [[ ! -f $ref_gen ]]; then
        wget -P $ref_dir http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz
        gunzip ${ref_gen}.gz
    fi

    ## additional fa not included in mouse
    if [[ ! -f $ref_bk ]]; then
        cp /data/CCBR_Pipeliner/iCLIP/index/active/2021_0607/mm10/01_source/BK000964.3_TPA_exp.fa $ref_bk
    fi

    ## gtf file
    if [[ ! -f $ref_gtf ]]; then
        wget -P $ref_dir http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
        gunzip ${ref_gtf}.gz
    fi
fi 

if [[ "$flag_index" == "Y" ]]; then
    echo "*** Indexing ***"

    #create index
    version_id="index"
    cmd=`star_index`
    echo $cmd > $log_dir/index.sh
    swarm -f $log_dir/index.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

fi

if [[ "$flag_variables" == "Y" ]]; then
    echo "*** Setting variables *** "

    if [[ -f variable_set.txt ]]; then rm variable_set.txt; fi
    touch variable_set.txt

    for i in {3..36}; do
        echo "version_id=\"$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)\"
        outFilterMismatchNoverReadLmax=\"$(awk -F'\t' -v i=$i 'NR == 2 {print $i}' $variable_file)\"
        outFilterMismatchNmax=\"$(awk -F'\t' -v i=$i 'NR == 3 {print $i}' $variable_file)\"
        outFilterMultimapNmax=\"$(awk -F'\t' -v i=$i 'NR == 4 {print $i}' $variable_file)\"
        outSJfilterReads=\"$(awk -F'\t' -v i=$i 'NR == 5 {print $i}' $variable_file)\"
        outSJfilterCountTotalMin=\"$(awk -F'\t' -v i=$i 'NR == 6 {print $i}' $variable_file)\"
        alignEndsType=\"$(awk -F'\t' -v i=$i 'NR == 7 {print $i}' $variable_file)\"
        outSAMunmapped=\"$(awk -F'\t' -v i=$i 'NR == 8 {print $i}' $variable_file)\"
        outSAMattributes=\"$(awk -F'\t' -v i=$i 'NR == 9 {print $i}' $variable_file)\"
        outFilterType=\"$(awk -F'\t' -v i=$i 'NR == 10 {print $i}' $variable_file)\"
        outFilterScoreMin=\"$(awk -F'\t' -v i=$i 'NR == 11 {print $i}' $variable_file)\"
        outFilterMultimapScoreRange=\"$(awk -F'\t' -v i=$i 'NR == 12 {print $i}' $variable_file)\"
        sjdbScore=\"$(awk -F'\t' -v i=$i 'NR == 13 {print $i}' $variable_file)\"
        outSJfilterOverhangMin=\"$(awk -F'\t' -v i=$i 'NR == 14 {print $i}' $variable_file)\"
        alignIntronMax=\"$(awk -F'\t' -v i=$i 'NR == 15 {print $i}' $variable_file)\"
        outFilterMatchNmin=\"$(awk -F'\t' -v i=$i 'NR == 16 {print $i}' $variable_file)\"
        outFilterMatchNminOverLread=\"$(awk -F'\t' -v i=$i 'NR == 17 {print $i}' $variable_file)\"
        alignSJDBoverhangMin=\"$(awk -F'\t' -v i=$i 'NR == 18 {print $i}' $variable_file)\"
        alignSJoverhangMin=\"$(awk -F'\t' -v i=$i 'NR == 19 {print $i}' $variable_file)\"
        winAnchorMultimapNmax=\"$(awk -F'\t' -v i=$i 'NR == 31 {print $i}' $variable_file)\"
        seedPerReadNmax=\"$(awk -F'\t' -v i=$i 'NR == 32 {print $i}' $variable_file)\"
        seedPerWindowNmax=\"$(awk -F'\t' -v i=$i 'NR == 33 {print $i}' $variable_file)\"
        alignWindowsPerReadNmax=\"$(awk -F'\t' -v i=$i 'NR == 34 {print $i}' $variable_file)\"
        alignTranscriptsPerReadNmax=\"$(awk -F'\t' -v i=$i 'NR == 35 {print $i}' $variable_file)\"
        seedMultimapNmax=\"$(awk -F'\t' -v i=$i 'NR == 36 {print $i}' $variable_file)\"
        seedNoneLociPerWindow=\"$(awk -F'\t' -v i=$i 'NR == 37 {print $i}' $variable_file)\"
        cmd=\`star_run\`
        echo -e \"--\$version_id:\n$cmd\"
        swarm -f \$log_dir/complete_\$version_id.sh --job-name \$version_id --merge-output --logdir \$log_dir --partition=norm -g 40 -t 4 --time 12:00:00
        " >> variable_set.txt

        cp variable_set.txt $output_dir/variable_set.txt
    done
fi

if [[ "$flag_align_partial" == "Y" ]]; then
    echo "*** Running Variations Partial***"
    
    #set group variables
    round_of_test="partial_sample"
    input_fq="FLAG_Ro_fclip.dedup.si.Sptan1.unique.fastq"

    # version_id="Default"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="20"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="3 1 1 1"
    #     alignEndsType="Local"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="2"
    #     outSJfilterOverhangMin="30 12 12 12"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="3"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="1A"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="1"
    #     outSJfilterReads="Unique"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="1b"
    #     outFilterMismatchNoverReadLmax="0.05"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="1"
    #     outSJfilterReads="Unique"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="1c"
    #     outFilterMismatchNoverReadLmax="0.02"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="1"
    #     outSJfilterReads="Unique"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2a"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2b"
    #     outFilterMismatchNoverReadLmax="0.05"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2c"
    #     outFilterMismatchNoverReadLmax="0.02"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2d"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2e"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="0"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2f"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2g"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2h"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="20 10 10 10 10"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="2i"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="EndToEnd"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
        
    # version_id="Encode1"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="30"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="EndToEnd"
    #     outSAMunmapped="Within"
    #     outSAMattributes="All"
    #     outFilterType="BySJout"
    #     outFilterScoreMin="10"
    #     outFilterMultimapScoreRange="1"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="1"
    #     alignSJoverhangMin="8"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd
    
    # # version_id="Encode2"
    # #         outFilterMismatchNoverReadLmax="0.04"
    # #         outFilterMismatchNmax="999"
    # #         outFilterMultimapNmax="1"
    # #         outSJfilterReads="All"
    # #         outSJfilterCountTotalMin="10 5 5 5"
    # #         alignEndsType="EndToEnd"
    # #         outSAMunmapped="Within"
    # #         outSAMattributes="All"
    # #         outFilterType="BySJout"
    # #         outFilterScoreMin="10"
    # #         outFilterMultimapScoreRange="1"
    # #         sjdbScore="0"
    # #         outSJfilterOverhangMin="15 6 6 6 6"
    # #         alignIntronMax="1000000"
    # #         outFilterMatchNmin="10"
    # #         outFilterMatchNminOverLread="0.6"
    # #         alignSJDBoverhangMin="1"
    # #         alignSJoverhangMin="8"
    # #         cmd=`star_run`
    # #         echo -e "--$version_id:\n"
    # #         $cmd
        
    # version_id="ClipSeqTools_v1"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="20"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="3 1 1 1"
    #     alignEndsType="Local"
    #     outSAMunmapped="None"
    #     outSAMattributes="All"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="0"
    #     sjdbScore="2"
    #     outSJfilterOverhangMin="30 12 12 12"
    #     alignIntronMax="50000"
    #     outFilterMatchNmin="15"
    #     outFilterMatchNminOverLread="0.9"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd

    # version_id="ClipSeqTools_v2"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="20"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="3 1 1 1"
    #     alignEndsType="Local"
    #     outSAMunmapped="None"
    #     outSAMattributes="All"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="0"
    #     sjdbScore="2"
    #     outSJfilterOverhangMin="30 12 12 12"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="15"
    #     outFilterMatchNminOverLread="0.9"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     $cmd

    version_id="clipv1_double_raiseq"
            outFilterMismatchNoverReadLmax="0.1"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="30000"
            alignTranscriptsPerReadNmax="30000"
            seedMultimapNmax="30000"
            seedNoneLociPerWindow="30"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            $cmd
            samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam   

fi

if [[ "$flag_pysam_partial_sample_old" == "Y" ]]; then
    echo "*** Running pysam ***"

    cigar_dir=$output_dir/partial_sample/cigar

    #remove prev counts
    rm $cigar_dir/*.csv    

    #skip encode2
    for i in {3..16} 18 19; do
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"
        bam_file=$output_dir/partial_sample/$version_id/${version_id}_Aligned.sortedByCoord.out.bam

        python cigar_pysam.py star_$version_id $bam_file $cigar_dir

    done

    for bam_name in ${novo_list[@]}; do
        version_id=`echo $bam_name | sed -n 's/FLAG_Ro_fclip.dedup.si.Sptan1.unique.parameter_//p' | sed -n 's/\.u.*//p'`
        echo "--$version_id"
        bam_file=$output_dir/partial_sample/novo/$bam_name

        #create individual cigar ratio files
        python cigar_pysam.py novo_$version_id $bam_file $cigar_dir
    done

    #merge files, create counts
    cat $cigar_dir/*.csv > $cigar_dir/final_cigar.csv
    cat $cigar_dir/final_cigar.csv | uniq -c > $cigar_dir/count_cigar.csv 

fi

if [[ "$flag_gap_partial" == "Y" ]]; then
    echo "*** Running GAP analysis***"

    #set dir, counts
    bam_dir=$output_dir/partial_sample/
    cigar_dir=$output_dir/partial_sample/cigar
    gap_min=50

    #mkdir, remove prev counts
    if [[ ! -d "$cigar_dir" ]]; then mkdir -p "$cigar_dir"; fi
    if [[ -f $cigar_dir/cigar_tmp.txt ]]; then rm $cigar_dir/cigar_tmp.txt; fi
    if [[ -f $cigar_dir/final_counts.csv ]]; then rm $cigar_dir/final_counts.csv; fi
    touch $cigar_dir/final_counts.csv

    for i in {3..15} 19 20; do
    
        counter=1
        
        #set verion, bam
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"
        bam_file=$bam_dir/$version_id/${version_id}_Aligned.sortedByCoord.out.bam

        #grep cigar strings with gaps
        samtools view -h -F 4 $bam_file | awk '$6 ~ /[0-9]*M[0-9]*N[0-9]*M/' | awk '{ print $6 }' | \
            sort -n | uniq -c | sed 's/  //g' | sed 's/ /\t/g' > $cigar_dir/cigar_tmp.txt

        #read gap file
        IFS=$'\n' read  -d '' -r -a cigar_array < $cigar_dir/cigar_tmp.txt
        
        for cigar_line in ${cigar_array[@]}; do
            #first line is counts
            if [[ "$counter" -eq 1 ]]; then
                line_count=$cigar_line
                counter=2
            else
                #pull the first match len, gap len, second match len
                match1_len=`echo $cigar_line | awk -F'M' '{ print $1 }'`
                gap_len=`echo $cigar_line | grep -o "[0-9]*N" | sed 's/N//'`
                match2_len=`echo $cigar_line | awk -F'M' '{ print $2}' | awk -F'N' '{ print $2}'`

                if [[ "$gap_len" -gt "$gap_min" ]]; then
                    #set num (min) and denom (max)
                    if [[ "$match1_len" -gt "$match2_len" ]]; then
                        min_val=$match2_len
                        max_val=$match1_len
                    else
                        min_val=$match1_len
                        max_val=$match2_len
                    fi

                    #determine match length ratio
                    ratio=`echo "scale=3; $min_val / $max_val" | bc -l`

                    #echo "line: $line_count | match1: $match1_len | gap1: $gap_len | match2: $match2_len | ratio: $ratio"
                    echo "$line_count,$ratio,$gap_len,$min_val,$max_val,$version_id" >> $cigar_dir/final_counts.csv
                fi
                counter=1
            fi
        done
    done
fi

if [[ "$flag_align_complete" == "Y" ]]; then
    echo "*** Running Variations Partial***"
    
    #set group variables
    round_of_test="complete_sample"
    input_fq="FLAG_Ro_fclip_filtered.fastq" #/data/RBL_NCI/Wolin/mES_fclip_1_YL_012122/01_preprocess/

    # version_id="2e"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="0"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     echo $cmd > $log_dir/complete_$version_id.sh
    #     swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
        
    # version_id="2h"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="Extend5pOfRead1"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="20 10 10 10 10"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     echo $cmd > $log_dir/complete_$version_id.sh
    #     swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
        
    # version_id="2i"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="EndToEnd"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     echo $cmd > $log_dir/complete_$version_id.sh
    #     swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="2j"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="999"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="10 5 5 5"
    #     alignEndsType="EndToEnd"
    #     outSAMunmapped="None"
    #     outSAMattributes="Standard"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="5"
    #     sjdbScore="0"
    #     outSJfilterOverhangMin="15 6 6 6 6"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="10"
    #     outFilterMatchNminOverLread="0.6"
    #     alignSJDBoverhangMin="3"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     echo $cmd > $log_dir/complete_$version_id.sh
    #     swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="ClipSeqTools_v1"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="20"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="3 1 1 1"
    #     alignEndsType="Local"
    #     outSAMunmapped="None"
    #     outSAMattributes="All"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="0"
    #     sjdbScore="2"
    #     outSJfilterOverhangMin="30 12 12 12"
    #     alignIntronMax="50000"
    #     outFilterMatchNmin="15"
    #     outFilterMatchNminOverLread="0.9"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
    #     cmd=`star_run`
    #     echo -e "--$version_id:\n"
    #     echo $cmd > $log_dir/complete_$version_id.sh
    #     swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="ClipSeqTools_v2"
    #     outFilterMismatchNoverReadLmax="0.04"
    #     outFilterMismatchNmax="999"
    #     outFilterMultimapNmax="20"
    #     outSJfilterReads="All"
    #     outSJfilterCountTotalMin="3 1 1 1"
    #     alignEndsType="Local"
    #     outSAMunmapped="None"
    #     outSAMattributes="All"
    #     outFilterType="Normal"
    #     outFilterScoreMin="0"
    #     outFilterMultimapScoreRange="0"
    #     sjdbScore="2"
    #     outSJfilterOverhangMin="30 12 12 12"
    #     alignIntronMax="1000000"
    #     outFilterMatchNmin="15"
    #     outFilterMatchNminOverLread="0.9"
    #     alignSJDBoverhangMin="5"
    #     alignSJoverhangMin="5"
        # cmd=`star_run`
        # echo -e "--$version_id:\n"
        # echo $cmd > $log_dir/complete_$version_id.sh
        # swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    version_id="2e_double"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="10 5 5 5"
            alignEndsType="Extend5pOfRead1"
            outSAMunmapped="None"
            outSAMattributes="Standard"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="0"
            outSJfilterOverhangMin="15 6 6 6 6"
            alignIntronMax="1000000"
            outFilterMatchNmin="10"
            outFilterMatchNminOverLread="0.6"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="20000"
            alignTranscriptsPerReadNmax="20000"
            seedMultimapNmax="20000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/complete_$version_id.sh
            samtools index $output_dir/complete_sample/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
            swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    version_id="clipv1_double"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="20000"
            alignTranscriptsPerReadNmax="20000"
            seedMultimapNmax="20000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/complete_$version_id.sh
            samtools index $output_dir/complete_sample/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
            swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    version_id="2e_triple"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="10 5 5 5"
            alignEndsType="Extend5pOfRead1"
            outSAMunmapped="None"
            outSAMattributes="Standard"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="0"
            outSJfilterOverhangMin="15 6 6 6 6"
            alignIntronMax="1000000"
            outFilterMatchNmin="10"
            outFilterMatchNminOverLread="0.6"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="30000"
            alignTranscriptsPerReadNmax="30000"
            seedMultimapNmax="30000"
            seedNoneLociPerWindow="30"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/complete_$version_id.sh
            samtools index $output_dir/complete_sample/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
            swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
        
    version_id="clipv1_triple"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="30000"
            alignTranscriptsPerReadNmax="30000"
            seedMultimapNmax="30000"
            seedNoneLociPerWindow="30"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/complete_$version_id.sh
            samtools index $output_dir/complete_sample/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
            swarm -f $log_dir/complete_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
fi

if [[ "$flag_gap_complete" == "Y" ]]; then
    echo "*** Running GAP analysis***"

    #set dir, counts
    bam_dir=$output_dir/complete_sample
    cigar_dir=$output_dir/complete_sample/cigar

    #set min gap length to consider
    gap_min=100

    #mkdir
    if [[ ! -d "$cigar_dir" ]]; then mkdir -p "$cigar_dir"; fi
    
    #for each sample 
    for i in 11 14 15 16 19 20 27 28 30; do
        #set verion
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"

        #set bam file
        bam_file=$bam_dir/$version_id/${version_id}_Aligned.sortedByCoord.out.bam

        #grep cigar strings with gaps, print strings
        #remove [0-9]S at the front or back of string, leaving only [0-9]M[0-9]N[0-9]M
        #sort and count unique values, remove large and leading spaces, turn into text file with \t's
        full_cigar_file=$cigar_dir/${version_id}_fullcigars.txt
        if [[ ! -f "$full_cigar_file" ]]; then
            echo "----processing bam file"
            samtools view -h -F 4 $bam_file | awk '$6 ~ /[0-9]*M[0-9]*N[0-9]*M/' | awk '{ print $6 }' | \
            sed 's/[0-9]*S//g' | sort -n | uniq -c | sed 's/  //g' | sed 's/^ *//' | sed 's/ /\t/g' > $full_cigar_file
        fi

        #read cigar file
        IFS=$'\n' read  -d '' -r -a cigar_array < $full_cigar_file
        counter=1
        
        #prep intermed files
        gap_counts_file=$cigar_dir/${version_id}_gap_counts.txt
        gap_complete_file=$cigar_dir/${version_id}_gap_analysis.txt

        if [[ ! -f $gap_counts_file ]] || [[ ! -f $gap_complete_file ]]; then
            echo "----creating gap counts files file"

            #create ouput files
            touch $gap_counts_file
            touch $gap_complete_file

            for cigar_line in ${cigar_array[@]}; do
                #first line is counts
                if [[ "$counter" -eq 1 ]]; then
                    read_count=$cigar_line

                    #ready for read analysis
                    counter=2
                else
                    #determine number of gaps in each read
                    n_gaps=`echo $cigar_line | tr -cd 'N' | wc -c`

                    #set 
                    gaps_gt_min=0

                    # for each gap determine the lengths of alignments before
                    # and after each read. determine the ratio of these alignments
                    # save all data
                    for (( gap_counter=1; gap_counter<=$n_gaps; gap_counter++ )); do
                        
                        # pull the correct alignments section depending on the position of the gap
                        # for example 1M2N3M4N5M6N7M would have
                        # gap_counter1: current_line = 1M2N3M
                        # gap_counter2: current_line = 3M4N5M
                        if [[ "$gap_counter" -eq 1 ]]; then
                            current_line=`echo $cigar_line | grep -o "[0-9]*M[0-9]*N[0-9]*M" | head -n1`
                        elif [[ "$gap_counter" -eq 2 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        elif [[ "$gap_counter" -eq 3 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        elif [[ "$gap_counter" -eq 4 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        else
                            echo "there are a total of $n_gaps gaps in this read ($current_line) - reconfigure"
                        fi

                        #pull the left match len, gap len, right match len
                        match1_len=`echo $current_line | awk -F'M' '{ print $1 }'`
                        gap_len=`echo $current_line | grep -o "[0-9]*N" | sed 's/N//'`
                        match2_len=`echo $current_line | awk -F'M' '{ print $2}' | awk -F'N' '{ print $2}'`

                        #only count gap if it is greater than min
                        if [[ "$gap_len" -gt "$gap_min" ]]; then
                            #determine which (left or right) match is longer
                            #set num (min) and denom (max)
                            if [[ "$match1_len" -gt "$match2_len" ]]; then
                                min_val=$match2_len
                                max_val=$match1_len
                            else
                                min_val=$match1_len
                                max_val=$match2_len
                            fi

                            #determine match length ratio
                            ratio=`echo "scale=3; $min_val / $max_val" | bc -l`

                            echo "$read_count,$ratio,$gap_len,$min_val,$max_val,$version_id" >> $gap_complete_file

                            #mark read as having +1 gap greater than min value set
                            gaps_gt_min=$((gaps_gt_min+1))
                        fi
                    done

                    #save gaps_gt_min to count file
                    echo "$read_count,$n_gaps,$gaps_gt_min,$version_id" >> $gap_counts_file
                    
                    #reset to 1, ready for read count
                    counter=1
                fi
            done
        fi
    done
fi

if [[ "$flag_align_stats" == "Y" ]]; then
    echo "*** Running alignment ***"

    #set output
    bam_dir=$output_dir/complete_sample/
    align_dir=$output_dir/complete_sample/align
    alignment_file=$align_dir/alignment_stats.txt

    #create alignment file
    echo "version_id,total,mapped,unmapped,unique,mm" > $alignment_file

    #for i in 11; do
    for i in 11 14 15 16 19 20 28 30; do
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"
        mapped_bam_file=$bam_dir/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
        unmapped_bam_file=$bam_dir/$version_id/${version_id}_Unmapped.out.mate1

        #gather alignment stats
        echo "----alignment stats"
        total=$(samtools view -c $mapped_bam_file)
        mapped=$(samtools view -c -F 4 $mapped_bam_file)
        unmapped=$(cat $unmapped_bam_file | grep "^@" | wc -l)
        reads_uni=$(samtools view $mapped_bam_file | grep "NH:i:1$" | wc -l)
        reads_mm=$(samtools view $mapped_bam_file | grep -v "NH:i:1$" | wc -l)
        echo "$version_id,$total,$mapped,$unmapped,$reads_uni,$reads_mm,$reads_mm2" >> $alignment_file
    done
fi

if [[ "$flag_subset_rnu" == "Y" ]]; then
    echo "Subsetting for RNU6"
    
    input_fq=${input_dir}/FLAG_Ro_fclip_filtered.fastq
    read_file_star=/home/sevillas2/git/iCLIP/build/STAR_testing/rnu6_readids.txt
    read_file_novo=/home/sevillas2/git/iCLIP/build/STAR_testing/rnu6_readis_piccard.txt
    rnu_fq=${input_dir}/rnu6.fastq
    novo_bam=/data/RBL_NCI/Wolin/mES_fclip_1_YL_012122/alignment_analysis/star/output/complete_sample/novo_complete/FLAG_Ro_fclip.mapq_recalculated.bam
    filt_bam=${output_dir}/rnu_sample/novo/novo_rnu6.bam

    touch $rnu_fq

    if [[ ! -f $rnu_fq ]]; then
        echo "creating FQ"
        IFS=$'\n' read -d '' -r -a read_array < $read_file_star
        
        for f in ${read_array[@]}; do
            echo $f
            grep -A 3 $f $input_fq >> $rnu_fq
        done
    fi

    if [[ ! -f $filt_bam ]]; then
        echo "creating NOVO"
        module load picard
        java -jar $PICARDJARPATH/picard.jar FilterSamReads I=$novo_bam O=$filt_bam READ_LIST_FILE=$read_file_novo FILTER=includeReadList
    fi
fi

if [[ "$flag_align_rnu" == "Y" ]]; then
    echo "*** Running Variations Partial ***"

    #set group variables
    round_of_test="rnu_sample"
    input_fq="rnu6.fastq"
   
    version_id="2e"
        outFilterMismatchNoverReadLmax="0.04"
        outFilterMismatchNmax="999"
        outFilterMultimapNmax="999"
        outSJfilterReads="All"
        outSJfilterCountTotalMin="10 5 5 5"
        alignEndsType="Extend5pOfRead1"
        outSAMunmapped="None"
        outSAMattributes="Standard"
        outFilterType="Normal"
        outFilterScoreMin="0"
        outFilterMultimapScoreRange="0"
        sjdbScore="0"
        outSJfilterOverhangMin="15 6 6 6 6"
        alignIntronMax="1000000"
        outFilterMatchNmin="10"
        outFilterMatchNminOverLread="0.6"
        alignSJDBoverhangMin="5"
        alignSJoverhangMin="5"
        cmd=`star_run`
        echo -e "--$version_id:\n"
        $cmd
        samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam

    version_id="ClipSeqTools_v1"
        outFilterMismatchNoverReadLmax="0.04"
        outFilterMismatchNmax="999"
        outFilterMultimapNmax="20"
        outSJfilterReads="All"
        outSJfilterCountTotalMin="3 1 1 1"
        alignEndsType="Local"
        outSAMunmapped="None"
        outSAMattributes="All"
        outFilterType="Normal"
        outFilterScoreMin="0"
        outFilterMultimapScoreRange="0"
        sjdbScore="2"
        outSJfilterOverhangMin="30 12 12 12"
        alignIntronMax="50000"
        outFilterMatchNmin="15"
        outFilterMatchNminOverLread="0.9"
        alignSJDBoverhangMin="5"
        alignSJoverhangMin="5"
        cmd=`star_run`
        echo -e "--$version_id:\n"
        $cmd
        samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam


    version_id="2e_git"
        outFilterMismatchNoverReadLmax="0.04"
        outFilterMismatchNmax="999"
        outFilterMultimapNmax="10000"
        outSJfilterReads="All"
        outSJfilterCountTotalMin="10 5 5 5"
        alignEndsType="Extend5pOfRead1"
        outSAMunmapped="None"
        outSAMattributes="Standard"
        outFilterType="Normal"
        outFilterScoreMin="0"
        outFilterMultimapScoreRange="0"
        sjdbScore="0"
        outSJfilterOverhangMin="15 6 6 6 6"
        alignIntronMax="1000000"
        outFilterMatchNmin="10"
        outFilterMatchNminOverLread="0.6"
        alignSJDBoverhangMin="5"
        alignSJoverhangMin="5"
        winAnchorMultimapNmax="10000"
        seedPerReadNmax="10000"
        seedPerWindowNmax="10000"
        alignWindowsPerReadNmax="10000"
        alignTranscriptsPerReadNmax="10000"
        seedMultimapNmax="10000"
        seedNoneLociPerWindow="10"
        cmd=`star_run_expanded`
        echo -e "--$version_id:\n"
        $cmd
        samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
        
    version_id="clipv1_git"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="10000"
            alignTranscriptsPerReadNmax="10000"
            seedMultimapNmax="10000"
            seedNoneLociPerWindow="10"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            $cmd
            samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
            
    version_id="2e_double"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="10 5 5 5"
            alignEndsType="Extend5pOfRead1"
            outSAMunmapped="None"
            outSAMattributes="Standard"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="0"
            outSJfilterOverhangMin="15 6 6 6 6"
            alignIntronMax="1000000"
            outFilterMatchNmin="10"
            outFilterMatchNminOverLread="0.6"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="20000"
            alignTranscriptsPerReadNmax="20000"
            seedMultimapNmax="20000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            $cmd
            samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
        
    version_id="clipv1_double"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="20000"
            alignTranscriptsPerReadNmax="20000"
            seedMultimapNmax="20000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            $cmd
            samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
        
    version_id="2e_lowerq"
            outFilterMismatchNoverReadLmax="0.02"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="10 5 5 5"
            alignEndsType="Extend5pOfRead1"
            outSAMunmapped="None"
            outSAMattributes="Standard"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="0"
            outSJfilterOverhangMin="15 6 6 6 6"
            alignIntronMax="1000000"
            outFilterMatchNmin="10"
            outFilterMatchNminOverLread="0.6"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="10000"
            alignTranscriptsPerReadNmax="10000"
            seedMultimapNmax="10000"
            seedNoneLociPerWindow="10"
            cmd=`star_star_run_expandedrun`
            echo -e "--$version_id:\n"
            $cmd
            samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
        
    version_id="clipv1_lowerq"
            outFilterMismatchNoverReadLmax="0.02"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="1000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="5"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="10000"
            alignTranscriptsPerReadNmax="10000"
            seedMultimapNmax="10000"
            seedNoneLociPerWindow="10"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            $cmd
            samtools index $output_dir/$round_of_test/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
fi

if [[ "$flag_align_stats_rnu" == "Y" ]]; then
    echo "*** Running alignment ***"

    #set output
    round_of_test="rnu_sample"
    bam_dir=$output_dir/$round_of_test
    align_dir=$output_dir/$round_of_test/align
    alignment_file=$align_dir/alignment_stats.txt
    log_file=$align_dir/runtime_stats.txt

    #mkdir, clear file
    if [[ ! -d "$align_dir" ]]; then mkdir -p "$align_dir"; fi

    #create alignment,log file
    echo "version_id,total,mapped,unmapped,unique,mm" > $alignment_file
    touch $log_file

    for i in 11 19 21 25 26 27 28 29 30 31; do
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"
        mapped_bam_file=$bam_dir/$version_id/${version_id}_Aligned.sortedByCoord.out.bam
        unmapped_bam_file=$bam_dir/$version_id/${version_id}_Unmapped.out.mate1
        #echo "$mapped_bam_file"
        #echo "--$version_id"

        #gather alignment stats
        echo "----alignment stats"
        total=$(samtools view -c $mapped_bam_file)
        mapped=$(samtools view -c -F 4 $mapped_bam_file)
        unmapped=$(cat $unmapped_bam_file | grep "^@" | wc -l)
        reads_uni=$(samtools view $mapped_bam_file | grep "NH:i:1$" | wc -l)
        reads_mm=$(samtools view $mapped_bam_file | grep -v "NH:i:1$" | wc -l)
        echo "$version_id,$total,$mapped,$unmapped,$reads_uni,$reads_mm,$reads_mm2" >> $alignment_file

        #gather run stats
        echo "----run stats"
        echo "$version_id" >> $log_file
        cat $bam_dir/$version_id/${version_id}_Log.final.out | grep "Started job on" >> $log_file
        cat $bam_dir/$version_id/${version_id}_Log.final.out | grep "Finished on" >> $log_file
        echo >> $log_file
    done
fi

if [[ "$flag_gap_rnu" == "Y" ]]; then
    echo "*** Running GAP analysis***"

    #set dir, counts
    bam_dir=$output_dir/complete_sample
    cigar_dir=$output_dir/complete_sample/cigar

    #set min gap length to consider
    gap_min=100

    #mkdir
    if [[ ! -d "$cigar_dir" ]]; then mkdir -p "$cigar_dir"; fi
    
    #for each sample 
    for i in 11 14 15 16 19 20; do
        #set verion
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"

        #set bam file
        bam_file=$bam_dir/$version_id/${version_id}_Aligned.sortedByCoord.out.bam

        #grep cigar strings with gaps, print strings
        #remove [0-9]S at the front or back of string, leaving only [0-9]M[0-9]N[0-9]M
        #sort and count unique values, remove large and leading spaces, turn into text file with \t's
        full_cigar_file=$cigar_dir/${version_id}_fullcigars.txt
        if [[ ! -f "$full_cigar_file" ]]; then
            echo "----processing bam file"
            samtools view -h -F 4 $bam_file | awk '$6 ~ /[0-9]*M[0-9]*N[0-9]*M/' | awk '{ print $6 }' | \
            sed 's/[0-9]*S//g' | sort -n | uniq -c | sed 's/  //g' | sed 's/^ *//' | sed 's/ /\t/g' > $full_cigar_file
        fi

        #read cigar file
        IFS=$'\n' read  -d '' -r -a cigar_array < $full_cigar_file
        counter=1
        
        #prep intermed files
        gap_counts_file=$cigar_dir/${version_id}_gap_counts.txt
        gap_complete_file=$cigar_dir/${version_id}_gap_analysis.txt

        if [[ ! -f $gap_counts_file ]] || [[ ! -f $gap_complete_file ]]; then
            echo "----creating gap counts files file"

            #create ouput files
            touch $gap_counts_file
            touch $gap_complete_file

            for cigar_line in ${cigar_array[@]}; do
                #first line is counts
                if [[ "$counter" -eq 1 ]]; then
                    read_count=$cigar_line

                    #ready for read analysis
                    counter=2
                else
                    #determine number of gaps in each read
                    n_gaps=`echo $cigar_line | tr -cd 'N' | wc -c`

                    #set 
                    gaps_gt_min=0

                    # for each gap determine the lengths of alignments before
                    # and after each read. determine the ratio of these alignments
                    # save all data
                    for (( gap_counter=1; gap_counter<=$n_gaps; gap_counter++ )); do
                        
                        # pull the correct alignments section depending on the position of the gap
                        # for example 1M2N3M4N5M6N7M would have
                        # gap_counter1: current_line = 1M2N3M
                        # gap_counter2: current_line = 3M4N5M
                        if [[ "$gap_counter" -eq 1 ]]; then
                            current_line=`echo $cigar_line | grep -o "[0-9]*M[0-9]*N[0-9]*M" | head -n1`
                        elif [[ "$gap_counter" -eq 2 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        elif [[ "$gap_counter" -eq 3 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        elif [[ "$gap_counter" -eq 4 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        else
                            echo "there are a total of $n_gaps gaps in this read ($current_line) - reconfigure"
                        fi

                        #pull the left match len, gap len, right match len
                        match1_len=`echo $current_line | awk -F'M' '{ print $1 }'`
                        gap_len=`echo $current_line | grep -o "[0-9]*N" | sed 's/N//'`
                        match2_len=`echo $current_line | awk -F'M' '{ print $2}' | awk -F'N' '{ print $2}'`

                        #only count gap if it is greater than min
                        if [[ "$gap_len" -gt "$gap_min" ]]; then
                            #determine which (left or right) match is longer
                            #set num (min) and denom (max)
                            if [[ "$match1_len" -gt "$match2_len" ]]; then
                                min_val=$match2_len
                                max_val=$match1_len
                            else
                                min_val=$match1_len
                                max_val=$match2_len
                            fi

                            #determine match length ratio
                            ratio=`echo "scale=3; $min_val / $max_val" | bc -l`

                            echo "$read_count,$ratio,$gap_len,$min_val,$max_val,$version_id" >> $gap_complete_file

                            #mark read as having +1 gap greater than min value set
                            gaps_gt_min=$((gaps_gt_min+1))
                        fi
                    done

                    #save gaps_gt_min to count file
                    echo "$read_count,$n_gaps,$gaps_gt_min,$version_id" >> $gap_counts_file
                    
                    #reset to 1, ready for read count
                    counter=1
                fi
            done
        fi
    done
fi

if [[ "$flag_subset_bam_multiple_genes" == "Y" ]]; then
    echo "Subsetting for Rnu6, Sptan1, sympk, Gm24204"
    #chr17:24,361,102-24,364,418
    module load picard

    read_file=/home/sevillas2/git/iCLIP/build/STAR_testing/subset_readids_piccard.txt

    declare -a sample_list=("2e" "2e_double" "ClipSeqTools_v1" "ClipSeqTools_double" "ClipSeqTools_triple" "novo_original")
 
    for sample_id in ${sample_list[@]}; do
        input_bam=${output_dir}/complete_sample/pipeline_v2/02_bam/02_merged/${sample_id}.unaware.merged.si.bam
        filt_bam=${output_dir}/subset_gene_list/${sample_id}_subset.bam

        if [[ ! -f $filt_bam ]]; then
            echo "--creating filtered ${sample_id}"
            java -jar $PICARDJARPATH/picard.jar FilterSamReads I=$input_bam O=$filt_bam READ_LIST_FILE=$read_file FILTER=includeReadList
        fi
    done
fi

if [[ "$flag_subset_fq_multiple_genes" == "Y" ]]; then
    module load samtools

    input_bam=${output_dir}/complete_sample/novo_complete/FLAG_Ro_fclip.mapq_recalculated.bam
    input_fq=${input_dir}/FLAG_Ro_fclip_filtered.fastq
    read_file=/home/sevillas2/git/iCLIP/build/STAR_testing/subset_readids_piccard.txt
    read_updated_file=/home/sevillas2/git/iCLIP/build/STAR_testing/multiplegenes_readids.txt
    mg_fq=${input_dir}/multiplegenes.fastq
    
    if [[ ! -f $read_updated_file ]]; then
        # get all readids for genes, add to list
        echo "Creating read list"
        cp $read_file $read_updated_file
        echo "" >> $read_updated_file
        
        #Gapdh chr6:125161721 - 125165734
        echo "--adding GAPDH"
        samtools view $input_bam | awk '$3 ~ /chr6/' | awk '$4 ~ /12516[1-5][0-9]*./' | awk '{ print $1 }' >> $read_updated_file

        #ACTB chr5:142903115-142906754
        echo "--adding ACTB"
        samtools view $input_bam | awk '$3 ~ /chr5/' | awk '$4 ~ /14290[3-6][0-9]*./' | awk '{ print $1 }' >> $read_updated_file
    fi

    if [[ ! -f $mg_fq ]]; then
        echo "creating FQ"
        sh_file="$log_dir/multiple_genes.sh"
        touch $sh_file
        echo "touch $mg_fq" >> $sh_file

        IFS=$'\n' read -d '' -r -a read_array < $read_updated_file
        for f in ${read_array[@]}; do
            echo "grep -A 3 "@$f" $input_fq >> $mg_fq" >> $sh_file
        done
        swarm -f $sh_file --job-name multiple_genes --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 1-00:00:00
    fi
fi

if [[ "$flag_align_overhang" == "Y" ]]; then
    echo "*** Running Variations Overhang ***"

    #set group variables
    round_of_test="expanded_gene_list"
    input_fq="multiplegenes.fastq"

    version_id="clipfinal_2"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="2"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="10000"
            alignTranscriptsPerReadNmax="10000"
            seedMultimapNmax="10000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/multigenes_$version_id.sh
            swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
    
        version_id="clipfinal_3"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="3"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="10000"
            alignTranscriptsPerReadNmax="10000"
            seedMultimapNmax="10000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/multigenes_$version_id.sh
            swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
    
        version_id="clipfinal_4"
            outFilterMismatchNoverReadLmax="0.04"
            outFilterMismatchNmax="999"
            outFilterMultimapNmax="10000"
            outSJfilterReads="All"
            outSJfilterCountTotalMin="3 1 1 1"
            alignEndsType="Local"
            outSAMunmapped="None"
            outSAMattributes="All"
            outFilterType="Normal"
            outFilterScoreMin="0"
            outFilterMultimapScoreRange="0"
            sjdbScore="2"
            outSJfilterOverhangMin="30 12 12 12"
            alignIntronMax="50000"
            outFilterMatchNmin="15"
            outFilterMatchNminOverLread="0.9"
            alignSJDBoverhangMin="4"
            alignSJoverhangMin="5"
            winAnchorMultimapNmax="10000"
            seedPerReadNmax="10000"
            seedPerWindowNmax="10000"
            alignWindowsPerReadNmax="10000"
            alignTranscriptsPerReadNmax="10000"
            seedMultimapNmax="10000"
            seedNoneLociPerWindow="20"
            cmd=`star_run_expanded`
            echo -e "--$version_id:\n"
            echo $cmd > $log_dir/multigenes_$version_id.sh
            swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
    
    # version_id="clipfinal_5"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="5"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="clipfinal_6"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="6"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
            
    # version_id="clipfinal_7"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="7"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="clipfinal_8"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="8"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="clipfinal_9"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="9"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="clipfinal_10"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="10"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="clipfinal_11"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="11"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

    # version_id="clipfinal_12"
    #         outFilterMismatchNoverReadLmax="0.04"
    #         outFilterMismatchNmax="999"
    #         outFilterMultimapNmax="10000"
    #         outSJfilterReads="All"
    #         outSJfilterCountTotalMin="3 1 1 1"
    #         alignEndsType="Local"
    #         outSAMunmapped="None"
    #         outSAMattributes="All"
    #         outFilterType="Normal"
    #         outFilterScoreMin="0"
    #         outFilterMultimapScoreRange="0"
    #         sjdbScore="2"
    #         outSJfilterOverhangMin="30 12 12 12"
    #         alignIntronMax="50000"
    #         outFilterMatchNmin="15"
    #         outFilterMatchNminOverLread="0.9"
    #         alignSJDBoverhangMin="12"
    #         alignSJoverhangMin="5"
    #         winAnchorMultimapNmax="10000"
    #         seedPerReadNmax="10000"
    #         seedPerWindowNmax="10000"
    #         alignWindowsPerReadNmax="10000"
    #         alignTranscriptsPerReadNmax="10000"
    #         seedMultimapNmax="10000"
    #         seedNoneLociPerWindow="20"
    #         cmd=`star_run_expanded`
    #         echo -e "--$version_id:\n"
    #         echo $cmd > $log_dir/multigenes_$version_id.sh
    #         swarm -f $log_dir/multigenes_$version_id.sh --job-name $version_id --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
fi

if [[ "$flag_gap_overhang" == "Y" ]]; then
    echo "*** Running GAP analysis***"

    #set dir, counts
    bam_dir=$output_dir/expanded_gene_list/bams
    cigar_dir=$output_dir/expanded_gene_list/cigar

    #set min gap length to consider
    gap_min=100

    #mkdir
    if [[ ! -d "$cigar_dir" ]]; then mkdir -p "$cigar_dir"; fi
    
    #for each sample 
    for i in {33..43}; do
        #set verion
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"

        #set bam file
        bam_file=$bam_dir/${version_id}_Aligned.sortedByCoord.out.bam

        #grep cigar strings with gaps, print strings
        #remove [0-9]S at the front or back of string, leaving only [0-9]M[0-9]N[0-9]M
        #sort and count unique values, remove large and leading spaces, turn into text file with \t's
        full_cigar_file=$cigar_dir/${version_id}_fullcigars.txt
        if [[ ! -f "$full_cigar_file" ]]; then
            echo "----processing bam file"
            samtools view -h -F 4 $bam_file | awk '$6 ~ /[0-9]*M[0-9]*N[0-9]*M/' | awk '{ print $6 }' | \
            sed 's/[0-9]*S//g' | sort -n | uniq -c | sed 's/  //g' | sed 's/^ *//' | sed 's/ /\t/g' > $full_cigar_file
        fi

        #read cigar file
        IFS=$'\n' read  -d '' -r -a cigar_array < $full_cigar_file
        counter=1
        
        #prep intermed files
        gap_counts_file=$cigar_dir/${version_id}_gap_counts.txt
        gap_complete_file=$cigar_dir/${version_id}_gap_analysis.txt

        if [[ ! -f $gap_counts_file ]] || [[ ! -f $gap_complete_file ]]; then
            echo "----creating gap counts files file"

            #create ouput files
            touch $gap_counts_file
            touch $gap_complete_file

            for cigar_line in ${cigar_array[@]}; do
                #first line is counts
                if [[ "$counter" -eq 1 ]]; then
                    read_count=$cigar_line

                    #ready for read analysis
                    counter=2
                else
                    #determine number of gaps in each read
                    n_gaps=`echo $cigar_line | tr -cd 'N' | wc -c`

                    #set 
                    gaps_gt_min=0

                    # for each gap determine the lengths of alignments before
                    # and after each read. determine the ratio of these alignments
                    # save all data
                    for (( gap_counter=1; gap_counter<=$n_gaps; gap_counter++ )); do
                        
                        # pull the correct alignments section depending on the position of the gap
                        # for example 1M2N3M4N5M6N7M would have
                        # gap_counter1: current_line = 1M2N3M
                        # gap_counter2: current_line = 3M4N5M
                        if [[ "$gap_counter" -eq 1 ]]; then
                            current_line=`echo $cigar_line | grep -o "[0-9]*M[0-9]*N[0-9]*M" | head -n1`
                        elif [[ "$gap_counter" -eq 2 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        elif [[ "$gap_counter" -eq 3 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        elif [[ "$gap_counter" -eq 4 ]]; then
                            current_line=`echo $cigar_line | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | cut -d "M" -f 2- | cut -d "N" -f 2- | grep -o "[0-9]*M[0-9]*N[0-9]*M"`
                        else
                            echo "there are a total of $n_gaps gaps in this read ($current_line) - reconfigure"
                        fi

                        #pull the left match len, gap len, right match len
                        match1_len=`echo $current_line | awk -F'M' '{ print $1 }'`
                        gap_len=`echo $current_line | grep -o "[0-9]*N" | sed 's/N//'`
                        match2_len=`echo $current_line | awk -F'M' '{ print $2}' | awk -F'N' '{ print $2}'`

                        #only count gap if it is greater than min
                        if [[ "$gap_len" -gt "$gap_min" ]]; then
                            #determine which (left or right) match is longer
                            #set num (min) and denom (max)
                            if [[ "$match1_len" -gt "$match2_len" ]]; then
                                min_val=$match2_len
                                max_val=$match1_len
                            else
                                min_val=$match1_len
                                max_val=$match2_len
                            fi

                            #determine match length ratio
                            ratio=`echo "scale=3; $min_val / $max_val" | bc -l`

                            echo "$read_count,$ratio,$gap_len,$min_val,$max_val,$version_id" >> $gap_complete_file

                            #mark read as having +1 gap greater than min value set
                            gaps_gt_min=$((gaps_gt_min+1))
                        fi
                    done

                    #save gaps_gt_min to count file
                    echo "$read_count,$n_gaps,$gaps_gt_min,$version_id" >> $gap_counts_file
                    
                    #reset to 1, ready for read count
                    counter=1
                fi
            done
        fi
    done
fi

if [[ "$flag_align_stats_overhang" == "Y" ]]; then
    #gather run stats
    echo "----run stats"

    # set output file
    align_stats="$output_dir/expanded_gene_list/align_stats.csv"
    if [[ -f "$align_stats" ]]; then rm $align_stats; fi
    touch $align_stats

    # create array to run greps
    declare -a grep_list=("Started mapping on" "Finished on" "Number of input reads" "Average input read length" "Uniquely mapped reads number" "Uniquely mapped reads %"
                            "Number of splices: Annotated (sjdb)" "Number of splices: Non-canonical" "Number of reads mapped to multiple loci" "% of reads mapped to multiple loci"
                            "% of reads mapped to too many loci" "% of reads unmapped: too many mismatches" "% of reads unmapped: too short" "% of reads unmapped: other")
    
    grep_results=()

    # for each variable in set
    for i in {33..43}; do
        # set verion
        version_id=$(awk -F'\t' -v i=$i 'NR == 1 {print $i}' $variable_file)
        echo "--$version_id"

        grep_results=()
        grep_results+=($version_id)
        
        # grep all data points
        for grep_id in "${grep_list[@]}"; do
            grep_out=`cat $output_dir/expanded_gene_list/$version_id/${version_id}_Log.final.out | grep "${grep_id}" | awk '{split($0, array, "|\t"); print array[2]}'`
            grep_results+=(,$grep_out)
        done

        # echo results
        echo "${grep_results[@]}" >> $align_stats
    done
fi