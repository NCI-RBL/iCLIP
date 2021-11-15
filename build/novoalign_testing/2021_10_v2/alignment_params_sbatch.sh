#!/bin/sh
#sh /home/sevillas2/git/iCLIP/build/novoalign_testing/2021_10_v2/alignment_params_sbatch.sh

option=$1

log_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/log"

##############################################################################################################################
#Run alignment
if [ "$option" == "align" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/align_sbatch.sh"
        touch $sbatch_job

        #create error dir
        if [[ ! -d "$log_dir/test${test_id}/err" ]]; then mkdir $log_dir/test${test_id}/err; fi

        for split_id in {01..10}; do \
            sh_job="$log_dir/test${test_id}/align_test${test_id}_${split_id}_sh.sh"
            
            echo "sbatch --job-name=mq.${test_id}.${split_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/align_test${test_id}.${split_id}.err --output=$log_dir/test${test_id}/err/align_test${test_id}.${split_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
        done
    done
fi

##############################################################################################################################
#Run cleanup
if [ "$option" == "cleanup" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/cleanup_sbatch.sh"
        touch $sbatch_job

        #create error dir
        if [[ ! -d "$log_dir/test${test_id}/err" ]]; then mkdir $log_dir/test${test_id}/err; fi

        for split_id in {01..10}; do \
            sh_job="$log_dir/test${test_id}/cleanup_test${test_id}_${split_id}_sh.sh"

            echo "sbatch --job-name=cu.${test_id}.${split_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/cu_test${test_id}.${split_id}.err --output=$log_dir/test${test_id}/err/cu_test${test_id}.${split_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
        done
    done
fi

##############################################################################################################################
#Run create Unique and MM
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/unique_mm"
job_id1="u"
job_id2="m"

if [ "$option" == "unique_mm" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/unique_mm_split_sbatch.sh"
        touch $sbatch_job

        for split_id in {01..10}; do \
            sh_job1="$log_dir/test${test_id}/unique_test${test_id}_${split_id}_sh.sh"
            sh_job2="$log_dir/test${test_id}/mm_test${test_id}_${split_id}_sh.sh"

            if [[ ! -f "$output_dir/test$test_id/test$test_id.split.$split_id.unique.si.bam" ]]; then
                echo "sbatch --job-name=${job_id1}.${test_id}.${split_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/${job_id1}_test${test_id}.${split_id}.err --output=$log_dir/test${test_id}/err/${job_id1}_test${test_id}.${split_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job1;" >> $sbatch_job
            fi

            if [[ ! -f "$output_dir/test$test_id/test$test_id.split.$split_id.mm.si.bam" ]]; then
                echo "sbatch --job-name=${job_id2}.${test_id}.${split_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/${job_id2}_test${test_id}.${split_id}.err --output=$log_dir/test${test_id}/err/${job_id2}_test${test_id}.${split_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job2;" >> $sbatch_job
            fi
        done
    done
fi

##############################################################################################################################
#Run merged unique and mm splits
job_id="splits"
if [ "$option" == "merge_splits" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/merge_splits_sbatch.sh"
        touch $sbatch_job
        
        sh_job="$log_dir/test${test_id}/merge_splits_test${test_id}_sh.sh"
        echo "sbatch --job-name=splits.${test_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/splits_test${test_id}.err --output=$log_dir/test${test_id}/err/splits_test${test_id}.out --mem=50g --gres=lscratch:50 --time 00-06:00:00 $sh_job;" >> $sbatch_job
    done
fi

##############################################################################################################################
#Run merged unique and mm
job_id="merge_um"
if [ "$option" == "$job_id" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/${job_id}_sbatch.sh"
        touch $sbatch_job
        
        sh_job="$log_dir/test${test_id}/${job_id}_test${test_id}_sh.sh"
        echo "sbatch --job-name=$job_id.${test_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/${job_id}_test${test_id}.err --output=$log_dir/test${test_id}/err/${job_id}_test${test_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
    done
fi

##############################################################################################################################
#Run dedup
job_id="dedup"
if [ "$option" == "dedup" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/dedup_sbatch.sh"
        touch $sbatch_job
        
        sh_job="$log_dir/test${test_id}/dedup_test${test_id}_sh.sh"
        echo "sbatch --job-name=de.${test_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/de_test${test_id}.${split_id}.err --output=$log_dir/test${test_id}/err/de_test${test_id}.${split_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
    done
fi

##############################################################################################################################
#Run beds
job_id="bed"
if [ "$option" == "beds" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/beds_sbatch.sh"
        touch $sbatch_job
        
        sh_job="$log_dir/test${test_id}/beds_test${test_id}_sh.sh"
        echo "sbatch --job-name=${job_id}.${test_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/${job_id}_test${test_id}.err --output=$log_dir/test${test_id}/err/${job_id}_test${test_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
    done
fi

##############################################################################################################################
#Run counts
job_id="c"
if [ "$option" == "counts" ]; then \
    for test_id in {1..9}; do \
        sbatch_job="$log_dir/test${test_id}/counts_sbatch.sh"
        touch $sbatch_job
        
        sh_job="$log_dir/test${test_id}/counts_test${test_id}_sh.sh"
        echo "sbatch --job-name=${job_id}.${test_id} --cpus-per-task=32 --verbose --error=$log_dir/test${test_id}/err/${job_id}_test${test_id}.err --output=$log_dir/test${test_id}/err/${job_id}_test${test_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
    done
fi