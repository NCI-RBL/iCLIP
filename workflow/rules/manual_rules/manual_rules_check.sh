#!/bin/sh

#set dir
log_dir=$1
project_dir=$2

#set args from command line
rulename=$3
sample_name=$4
maxn=$5

#mkdirs if necessary
if [[ ! -d "$log_dir" ]]; then mkdir $log_dir; fi

##############################################################################################################################
# alignment
output_dir="${project_dir}/01_preprocess/02_alignment"
declare -a a_types=("masked" "unmasked" "unaware")

if [ "$rulename" == "alignment" ]; then
    echo "****** $sample_name ***"
    echo "- Alignment Check Starting -"
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi
            
            output_file="$output_dir/${sample_name}.${a_id}.split.$i.sam.gz"
            if [[ ! -f ${output_file} ]]; then
                echo "> MISSING: ${output_file}"
            fi
        done
    done

    echo "- Check Complete -"
    echo
fi

##############################################################################################################################
# cleanup
output_dir="${project_dir}/01_preprocess/03_genomic"
declare -a a_types=("masked" "unmasked" "unaware")

if [ "$rulename" == "cleanup" ]; then
    echo "****** $sample_name ***"
    echo "- Cleanup Check Starting -"
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            #check file exists
            output_file="$output_dir/${sample_name}.${a_id}.split.$i.sam.gz"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done

    echo "- Check Complete -"
    echo
fi

##############################################################################################################################
# unmapped
output_dir="${project_dir}/02_bam/01_unmapped"
a_id="unmasked"
if [ "$rulename" == "unmapped" ]; then
    echo "****** $sample_name ***"
    echo "- Unmapped Check Starting"
    
    #check file exists
    output_file="$output_dir/${sample_name}.${a_id}.complete.bam"          
    if [[ ! -f $output_file ]]; then echo "> MISSING: $output_file"; fi
    
    echo "- Check Complete"
    echo
fi

##############################################################################################################################
# unique_mm
output_dir_u="${project_dir}/01_preprocess/05_unique"
output_dir_m="${project_dir}/01_preprocess/05_mm"

if [ "$rulename" == "unique_mm" ]; then
    echo "****** $sample_name ***"
    echo "- Unique Check Starting -"
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            #check file exists
            output_file="$output_dir_u/${sample_name}.${a_id}.split.${i}.unique.si.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done

    echo "- Check Complete -"
    echo
fi

if [ "$rulename" == "unique_mm" ]; then
    echo "****** $sample_name ***"
    echo "- MM Check Starting -"
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            #check file exists
            output_file="$output_dir_m/${sample_name}.${a_id}.split.${i}.mm.si.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done

    echo "- Check Complete -"
    echo
fi

##############################################################################################################################
# merge_splits
output_dir="${project_dir}/02_bam/02_merged"
if [ "$rulename" == "merge_splits" ]; then
    echo "****** $sample_name ***"
    echo "- Merge Splits Check Starting"
    
    #check file exists
    output_file="$output_dir_u/${sample_name}.${a_id}.merged.unique.si.bam"
    if [[ ! -f $output_file ]]; then echo "> MISSING: $output_file"; fi
    
    #check file exists
    output_file="$output_dir_u/${sample_name}.${a_id}.merged.mm.si.bam"
    if [[ ! -f $output_file ]]; then echo "> MISSING: $output_file"; fi
    

    echo "- Check Complete"
    echo
fi