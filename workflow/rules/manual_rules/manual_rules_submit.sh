#sh /data/RBL_NCI/Wolin/mES_fclip_1_YL/log/manual_rules/manual_rules_submit.sh

##set dir
log_dir=$1
project_dir=$2
doc="/data/CCBR_Pipeliner/iCLIP/container/USeq_8.9.6/Apps/SamTranscriptomeParser"

#set args from command line
rulename=$3
sample_name=$4
maxn=$5

#mkdirs if necessary
if [[ ! -d "$log_dir" ]]; then mkdir $log_dir; fi

##############################################################################################################################
### cleanup
if [ "$rulename" == "cleanup" ]; then
    sh $log_dir/cleanup_sbatch_${sample_name}.sh
fi

##############################################################################################################################
### unmapped
if [ "$rulename" == "unmapped" ]; then
    sh $log_dir/unmapped_sbatch_${sample_name}.sh
fi

##############################################################################################################################
### unique_mm
if [ "$rulename" == "unique_mm" ]; then
    sh $log_dir/unique_sbatch_${sample_name}.sh
    sh $log_dir/mm_sbatch_${sample_name}.sh
fi

##############################################################################################################################
### merge_splits
if [ "$rulename" == "merge_splits" ]; then
    sh $log_dir/merge_splits_sbatch_${sample_name}.sh
fi

##############################################################################################################################
### merge_um
if [ "$rulename" == "merge_um" ]; then
    sh $log_dir/merge_um_sbatch_${sample_name}.sh
fi