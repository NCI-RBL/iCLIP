#!/bin/sh

#sh /home/sevillas2/git/iCLIP/build/novoalign_testing/2021_10_v2/alignment_params_submit.sh v1 align

#set version
version=$1
if [[ $version == "v1" ]]; then 
    log_dir="/data/RBL_NCI/Wolin/Sam/novoalign/log"
    project_dir="/data/RBL_NCI/Wolin/Sam/novoalign"
else
    log_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/log"
    project_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2"
fi

#set option
option=$2

##############################################################################################################################
### alignment
if [ "$option" == "align" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/align_sbatch.sh
    done
fi

### cleanup
if [ "$option" == "cleanup" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/cleanup_sbatch.sh
    done
fi

### unique_mm (OLD)
if [ "$option" == "unique_mm_old" ]; then
    for testid in {1..1}; do \
        sh $log_dir/test${testid}/unique_mm_sbatch.sh
    done
fi

### unique_mm (NEW)
if [ "$option" == "unique_mm" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/unique_mm_split_sbatch.sh
    done
fi

### merge splits
if [ "$option" == "merge_splits" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/merge_splits_sbatch.sh
    done
fi

### merge unique mm
if [ "$option" == "merge_um" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/merge_um_sbatch.sh
    done
fi

### dedup
if [ "$option" == "dedup" ]; then \
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/dedup_sbatch.sh
    done
fi

### beds
if [ "$option" == "beds" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/beds_sbatch.sh
    done
fi

### counts
if [ "$option" == "counts" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/counts_sbatch.sh
    done
fi

### peak_annotations
if [ "$option" == "peak_annotations" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/peak_annotations_sbatch.sh
    done
fi

### annotation_report
if [ "$option" == "annotation_report" ]; then
    for testid in {1..9}; do \
        sh $log_dir/test${testid}/annotation_report_sbatch.sh
    done
fi