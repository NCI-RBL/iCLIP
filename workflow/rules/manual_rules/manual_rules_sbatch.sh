#!/bin/sh
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

#set alignment types
declare -a a_types=("masked" "unmasked" "unaware")

##############################################################################################################################
#Run cleanup
job_id="cu"
if [ "$rulename" == "cleanup" ]; then \
    
    #prep sbatch file
    sbatch_job="$log_dir/cleanup_sbatch_${sample_name}.sh"
    if [ -f $sbatch_job ]; then rm $sbatch_job; fi
    touch $sbatch_job

    #create error dir
    if [[ ! -d "$log_dir/err" ]]; then mkdir "$log_dir/err"; fi

    #add to sbatch
    for a_id in ${a_types[@]}; do    
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            sh_job="$log_dir/cleanup_${sample_name}.${a_id}.split.${i}_sh.sh"
            echo "sbatch --job-name=${job_id}.${sample_name}.${a_id}.${i} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
        done
    done
fi

##############################################################################################################################
#Run unmapped
job_id="unmap"
a_id="unmasked"
if [ "$rulename" == "unmapped" ]; then \
    
    #prep sbatch file
    sbatch_job="$log_dir/unmapped_sbatch_${sample_name}.sh"
    if [ -f $sbatch_job ]; then rm $sbatch_job; fi
    touch $sbatch_job

    #create error dir
    if [[ ! -d "$log_dir/err" ]]; then mkdir "$log_dir/err"; fi

    #add to sbatch
    sh_job="$log_dir/unmapped_${sample_name}.${a_id}_sh.sh"
    echo "sbatch --job-name=${job_id}.${sample_name}.${a_id} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
fi

##############################################################################################################################
#Run unique_mm
job_id="u"
if [ "$rulename" == "unique_mm" ]; then \
    
    #prep sbatch file
    sbatch_job="$log_dir/unique_sbatch_${sample_name}.sh"
    if [ -f $sbatch_job ]; then rm $sbatch_job; fi
    touch $sbatch_job

    #create error dir
    if [[ ! -d "$log_dir/err" ]]; then mkdir "$log_dir/err"; fi

    #add to sbatch
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            sh_job="$log_dir/unique_${sample_name}.${a_id}.split.${i}_sh.sh"
            echo "sbatch --job-name=${job_id}.${sample_name}.${a_id}.${i} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
        done
    done
fi

job_id="m"
if [ "$rulename" == "unique_mm" ]; then \
    
    #prep sbatch file
    sbatch_job="$log_dir/mm_sbatch_${sample_name}.sh"
    if [ -f $sbatch_job ]; then rm $sbatch_job; fi
    touch $sbatch_job

    #create error dir
    if [[ ! -d "$log_dir/err" ]]; then mkdir "$log_dir/err"; fi

    #add to sbatch
    for a_id in ${a_types[@]}; do   
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            sh_job="$log_dir/mm_${sample_name}.${a_id}.split.${i}_sh.sh"
            echo "sbatch --job-name=${job_id}.${sample_name}.${a_id}.${i} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
        done
    done
fi

##############################################################################################################################
# merge_splits
job_id="ms"
if [ "$rulename" == "merge_splits" ]; then \
    
    #prep sbatch file
    sbatch_job="$log_dir/merge_splits_sbatch_${sample_name}.sh"
    if [ -f $sbatch_job ]; then rm $sbatch_job; fi
    touch $sbatch_job

    #create error dir
    if [[ ! -d "$log_dir/err" ]]; then mkdir "$log_dir/err"; fi

    #add to sbatch
    for a_id in ${a_types[@]}; do    
        sh_job="$log_dir/merge_splits_u_${sample_name}_${a_id}_sh.sh"
        echo "sbatch --job-name=${job_id}.u.${sample_name}.${a_id} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job

        sh_job="$log_dir/merge_splits_m_${sample_name}_${a_id}_sh.sh"
        echo "sbatch --job-name=${job_id}.m.${sample_name}.${a_id}.${i} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
    done
fi

##############################################################################################################################
# merge_um
job_id="ms"
if [ "$rulename" == "merge_um" ]; then \
    
    #prep sbatch file
    sbatch_job="$log_dir/merge_um_sbatch_${sample_name}_${a_id}.sh"
    if [ -f $sbatch_job ]; then rm $sbatch_job; fi
    touch $sbatch_job

    #create error dir
    if [[ ! -d "$log_dir/err" ]]; then mkdir "$log_dir/err"; fi

    #add to sbatch
    for a_id in ${a_types[@]}; do    
        sh_job="$log_dir/merge_um_${sample_name}_${a_id}_sh.sh"
        echo "sbatch --job-name=${job_id}.${sample_name}.${a_id} --cpus-per-task=32 --verbose --error=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.err --output=$log_dir/err/${job_id}.${sample_name}.${a_id}.${i}.out --mem=50g --gres=lscratch:50 --time 01-00:00:00 $sh_job;" >> $sbatch_job
    done
fi