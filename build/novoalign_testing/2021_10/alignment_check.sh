#!/bin/sh

#sh /home/sevillas2/git/iCLIP/build/novoalign_testing/2021_10_v2/alignment_check.sh

option=$1

##############################################################################################################################
# Alignment
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/aligned"

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
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/cleanup"

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
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/unique_mm"

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
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/merge_splits"

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
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/merged"

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
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/dedup"

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
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/beds"
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
    echo "****** dedup Check Complete ******"
    echo
fi

##############################################################################################################################
# counts
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/counts"
all_dir="$output_dir/03_allreadpeaks"
unique_dir="$output_dir/03_uniquereadpeaks/"

if [ "$option" == "counts" ]; then
    echo
    echo "****** beds Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
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
    done
    echo "****** dedup Check Complete ******"
    echo
fi


##############################################################################################################################
# annotation
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/beds"
bed_dir="$output_dir/01_bed"
saf_dir="$output_dir/02_SAF"
if [ "$option" == "annotation" ]; then
    echo
    echo "****** Beds Check Starting ******"
    for test_id in {1..9}; do
        echo "*** test_$test_id ***"
        for split_id in {01..10}; do
            #check file exists
            output_file="$bed_dir/test${test_id}"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi

            #check file exists
            output_file="$saf_dir/test${test_id}"
            if [[ ! -f $output_file ]]; then echo "MISSING: $output_file"; fi
        done
    done
    echo "****** Beds Check Complete ******"
    echo
fi

