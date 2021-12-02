#!/bin/sh

#sh /home/sevillas2/git/iCLIP/build/novoalign_testing/2021_10_v2/alignment_check.sh v1 align

#set version
version=$1

if [[ $version == "v1" ]]; then 
    project_dir="/data/RBL_NCI/Wolin/Sam/novoalign"
else
    project_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2"
fi

#set option
option=$2

##############################################################################################################################
# Alignment
output_dir="${project_dir}/aligned"

if [ "$option" == "align" ]; then
    echo
    echo "****** Alignment Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
            #check file exists
            output_file="$output_dir/test${test_id}/Ro_Clip_2.split.$split_id.sam.gz"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done
    echo "****** Alignment Check Complete ******"
    echo
fi

##############################################################################################################################
# cleanup
output_dir="${project_dir}/cleanup"

if [ "$option" == "cleanup" ]; then
    echo
    echo "****** Cleanup Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
            #check file exists
            output_file="$output_dir/test$test_id/Ro_Clip_2.split.$split_id.sam.gz"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done
    echo "****** Cleanup Check Complete ******"
    echo
fi

##############################################################################################################################
# unique_mm
output_dir="${project_dir}/unique_mm"

if [ "$option" == "unique_mm" ]; then
    echo
    echo "****** Unique_mm Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
            #check file exists
            output_file="$output_dir/test$test_id/test$test_id.split.$split_id.unique.si.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

            #check file exists
            output_file="$output_dir/test$test_id/test$test_id.split.$split_id.mm.si.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done
    echo "****** Unique_mm Check Complete ******"
    echo
fi

##############################################################################################################################
# merge_splits
output_dir="${project_dir}/merge_splits"

if [ "$option" == "merge_splits" ]; then
    echo
    echo "****** merge_splits Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
            #check file exists
            output_file="$output_dir/test${test_id}/test${test_id}.merged.unique.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

            #check file exists
            output_file="$output_dir/test${test_id}/test${test_id}.merged.mm.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done
    echo "****** merge_splits Check Complete ******"
    echo
fi

##############################################################################################################################
# merged_mu
output_dir="${project_dir}/merged"

if [ "$option" == "merge_um" ]; then
    echo
    echo "****** merge_um Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
            #check file exists
            output_file="$output_dir/test${test_id}.merged.si.bam"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done
    echo "****** merge_um Check Complete ******"
    echo
fi

##############################################################################################################################
# dedup
output_dir="${project_dir}/dedup"

if [ "$option" == "dedup" ]; then
    echo
    echo "****** dedup Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        #check file exists
        output_file="$output_dir/test${test_id}.dedup.si.bam"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
    done
    echo "****** dedup Check Complete ******"
    echo
fi

##############################################################################################################################
# beds
output_dir="${project_dir}/beds"
bed_dir="$output_dir/01_bed"
saf_dir="$output_dir/02_SAF"

if [ "$option" == "beds" ]; then
    echo
    echo "****** beds Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        #check file exists
        output_file="$bed_dir/test${test_id}_all.bed"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

        #check file exists
        output_file="$bed_dir/test${test_id}_unique.bed"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

        #check file exists
        output_file="$saf_dir/test${test_id}_all.SAF"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

        #check file exists
        output_file="$saf_dir/test${test_id}_unique.SAF"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
    done
    echo "****** beds Check Complete ******"
    echo
fi

##############################################################################################################################
# counts
output_dir="${project_dir}/counts"
all_dir="$output_dir/03_allreadpeaks"
unique_dir="$output_dir/03_uniquereadpeaks"

if [ "$option" == "counts" ]; then
    echo
    echo "****** counts Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        
        #check file exists
        output_file="$all_dir/test${test_id}_uniqueCounts.txt"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

        #check file exists
        output_file="$all_dir/test${test_id}_allFracMMCounts.txt"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

        #check file exists
        output_file="$unique_dir/test${test_id}_uniqueCounts.txt"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

        #check file exists
        output_file="$unique_dir/test${test_id}_allFracMMCounts.txt"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
    done
    echo "****** counts Check Complete ******"
    echo
fi

##############################################################################################################################
# peak_annotations
output_dir="${project_dir}/annotations/02_peaks"
if [ "$option" == "peak_annotations" ]; then
    echo
    echo "****** peak_annotations Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        #check file exists
        output_file="$output_dir/test${test_id}_annotable.bed"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
    done
    echo "****** peak_annotations Check Complete ******"
    echo
fi

##############################################################################################################################
### annotation_report
output_dir="${project_dir}/annotations/03_reports"
if [ "$option" == "annotation_report" ]; then
    echo
    echo "****** annotation_report Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        
        #check file exists
        output_file="$output_dir/test${test_id}_annotation_final_report.html"
        if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
    done
    echo "****** annotation_report Check Complete ******"
    echo
fi