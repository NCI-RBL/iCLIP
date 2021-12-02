#!/bin/sh

# Script controls creation of manual rule SH, SBATCH files and performing
# file check on indiivudal rules, or submitting individual rules to the cluster
# Update the manual manifest with sample anme and max split number in order to run the scripts.

#example  input for sh, sbatch, combo
# sh /home/sevillas2/git/iCLIP/workflow/rules/manual_rules/manual_control.sh \
# /data/RBL_NCI/Wolin/mES_fclip_1_YL/log/manual_rules/log \
# /data/RBL_NCI/Wolin/mES_fclip_1_YL \
# sh

#example  input for check, submit
# sh /home/sevillas2/git/iCLIP/workflow/rules/manual_rules/manual_control.sh \
# /data/RBL_NCI/Wolin/mES_fclip_1_YL/log/manual_rules/log \
# /data/RBL_NCI/Wolin/mES_fclip_1_YL \
# check \
# cleanup


#set source_dir
source_dir="/home/sevillas2/git/iCLIP/workflow/rules/manual_rules"

#set args
log_dir=$1
project_dir=$2
option=$3
    #options 
    # - check runs only check command
    # - sh runs only sh command
    # - sbatch runs only sbatch command
    # - combo runs sh and sbatch command
    # - submit runs submit command, requires rule name
rule_name=$4 #required for check and submit

#set rules list
declare -a rule_list=("cleanup" "unmapped" "unique_mm" "merge_splits" "merge_um")

echo 
echo "Runnning $option"
echo

#check if dir exists, if not create
if [[ ! -d  $log_dir ]]; then mkdir -p $log_dir; fi
if [[ ! -d  $project_dir ]]; then mkdir -p $project_dir; fi

#copy manifest to log
cp ${source_dir}/manual_manifest.txt $log_dir
cp ${source_dir}/manual_rules_check.sh $log_dir
cp ${source_dir}/manual_rules_sh.sh $log_dir
cp ${source_dir}/manual_rules_sbatch.sh $log_dir

#read manifest
cat ${source_dir}/manual_manifest.txt | while read sample_name maxn; do

    # run check only 
    if [ $option == "check" ]; then
        if [ -z $rule_name ]; then
            echo
            echo "******************************"
            echo "Rule name is required for CHECK; Add ARG"
            echo
            exit
        else 
            sh ${source_dir}/manual_rules_check.sh \
            $log_dir \
            $project_dir \
            $rule_name \
            $sample_name \
            $maxn
        fi
    fi

    # run SH creation only
    if [ $option == "sh" ]; then
        for rule_name in ${rule_list[@]}; do
            sh ${source_dir}/manual_rules_sh.sh \
            $log_dir \
            $project_dir \
            $rule_name \
            $sample_name \
            $maxn
        done
    fi

    # run SBATCH creation only
    if [ $option == "sbatch" ]; then
        for rule_name in ${rule_list[@]}; do
            sh ${source_dir}/manual_rules_sbatch.sh \
            $log_dir \
            $project_dir \
            $rule_name \
            $sample_name \
            $maxn
        done
    fi

    # run SH creation AND SBATCH creation
    if [ $option == "combo" ]; then
        for rule_name in ${rule_list[@]}; do
            sh ${source_dir}/manual_rules_sh.sh \
            $log_dir \
            $project_dir \
            $rule_name \
            $sample_name \
            $maxn

            sh ${source_dir}/manual_rules_sbatch.sh \
            $log_dir \
            $project_dir \
            $rule_name \
            $sample_name \
            $maxn
        done
    fi

    # run submit ONLY
    if [ $option == "submit" ]; then
        if [ -z $rule_name ]; then
            echo
            echo "******************************"
            echo "Rule name is required for SUBMIT; Add ARG"
            echo
            exit
        else
            sh ${source_dir}/manual_rules_submit.sh \
            $log_dir \
            $project_dir \
            $rule_name \
            $sample_name \
            $maxn
        fi
    fi
done

