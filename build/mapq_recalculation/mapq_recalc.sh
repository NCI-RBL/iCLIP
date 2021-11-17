#!/bin/bash

input_dir="/data/RBL_NCI/Wolin/Sam/mapq_recalc/02_bam/04_mapq"
output_dir="/data/RBL_NCI/Wolin/Sam/mapq_recalc/analysis"

if [[ ! -d $output_dir ]]; then mkdir $output_dir; fi
if [[ ! -d $output_dir/counts ]]; then mkdir $output_dir/counts; fi

#use log to determine before / after / diff MAPQ recalculation results
#declare -a sample_list=("KO_fCLIP")
declare -a sample_list=("KO_fCLIP" "KO_NOPFA" "WT_fCLIP" "WT_NOPFA" "Y5KO_fCLIP" "Y5KO_NOPFA")
for sample in ${sample_list[@]}; do
    awk '{print $4}' "${input_dir}/${sample}.mapq_recalculated.log" | sort -n | uniq -c > "${output_dir}/counts/$sample.old_mapq.txt"
    awk '{print $6}' "$input_dir/$sample.mapq_recalculated.log" | sort -n | uniq -c > "${output_dir}/counts/$sample.new_mapq.txt"
    awk '{print $4,$6}' "$input_dir/$sample.mapq_recalculated.log" | sort -n | uniq -c > "${output_dir}/counts/$sample.change_mapq.txt"

    #submit to R to generate reports
    counts_old="$output_dir/counts/$sample.old_mapq.txt"
    counts_new="$output_dir/counts/$sample.new_mapq.txt"
    counts_diff="$output_dir/counts/$sample.change_mapq.txt"

    #module load R

    R_script="Rscript -e 'library(rmarkdown);
    rmarkdown::render(\"/home/sevillas2/git/iCLIP/workflow/scripts/07_mapq_stats.Rmd\",
    output_file = \"$output_dir/${sample}_report.html\",
    params= list(counts_old=\"$counts_old\",
        counts_new=\"$counts_new\",
        counts_diff=\"$counts_diff\",
        sampleid=\"$sample\",
        output_dir=\"$output_dir/counts/\"))'"
    
    echo $R_script

done

#sh /home/sevillas2/git/iCLIP/build/mapq_recalculation/mapq_recalc.sh