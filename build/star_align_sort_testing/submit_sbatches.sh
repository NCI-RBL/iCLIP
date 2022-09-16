#!/bin/sh
create="Y"
submit="Y"
run_list=("run1" "run2" "run3" "run4" "run5" "run6" "run7" "run8" "run9" "run10" "run11" "run12" "run13" "run14" "run15")
sample_list=("Y5KO_fCLIP" "WT_fCLIP")

# to run
run_list=("run14")
sample_list=("Y5KO_fCLIP")

if [[ $create == "Y" ]]; then
    for run_id in ${run_list[@]}; do
        for sample_id in ${sample_list[@]}; do
            sbatch_name="marco_${sample_id}_${run_id}.sh"
            cp sbatch_template_${run_id}.sh $sbatch_name
            sed -i "s/FILL_SAMPLE/$sample_id/g" $sbatch_name
            sed -i "s/FILL_RUN/$run_id/g" $sbatch_name
        done
    done
fi
if [[ $submit == "Y" ]]; then
    for run_id in ${run_list[@]}; do
        for sample_id in ${sample_list[@]}; do
            sbatch_name="marco_${sample_id}_${run_id}.sh"
            
            # submit sbatch
            sbatch --cpus-per-task=32 --verbose \
            --output=/data/sevillas2/marco_star/logs/%j.out \
            --mem=200g --gres=lscratch:800 --time 1-00:00:00 \
            --error=/data/sevillas2/marco_star/logs/%j.err \
            $sbatch_name 
        done
    done
fi

#            --mem=80g --gres=lscratch:800 --time 1-00:00:00 \
