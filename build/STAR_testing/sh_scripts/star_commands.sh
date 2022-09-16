#!/bin/bash

function create_sh()
{
    cp /data/sevillas2/star/logs/submit_star.sh $sh_file
    sed -i "s/fill_sample/$sample_id/g" $sh_file
    sed -i "s/fill_project/$project_id/g" $sh_file
    sed -i "s/fill_trial/$trial/g" $sh_file
    sed -i "s/fill_seed/$seedPerWindowNmax/g" $sh_file
    sed -i "s/fill_anchor/$winAnchorMultimapNmax/g" $sh_file
}

function set_variables()
{
    trial="${cpu}_${mem}_${seedPerWindowNmax}_${winAnchorMultimapNmax}"
    output_dir="/data/sevillas2/star/${sample_id}/${trial}"
    sh_file="/data/sevillas2/star/logs/${sample_id}_${trial}.sh"

    if [[ ! -d $output_dir ]]; then mkdir -p $output_dir; fi
}

function submit_sh()
{
    sbatch --cpus-per-task=${cpu} --verbose --output=/data/sevillas2/star/logs/%j.out --mem=${mem} --gres=lscratch:800 --time 10:00:00 --error=/data/sevillas2/star/logs/%j.err $sh_file
}

function stats()
{
    val=`cat $output_dir/$sample_id.out | grep "Uniquely mapped reads %"`
    descrip="unique%"
    cleanup
    val=`cat $output_dir/$sample_id.out | grep "Number of splices: Total"`
    descrip="splice_num"
    cleanup
    val=`cat $output_dir/$sample_id.out | grep "Number of splices: Annotated (sjdb)"`
    descrip="splice_anno_num"
    cleanup
    val=`cat $output_dir/$sample_id.out | grep "% of reads mapped to multiple loci"`
    descrip="mm%"
    cleanup
    val=`cat $output_dir/$sample_id.out | grep "% of reads mapped to too many loci"`
    descrip="unmapped_toomany_%"
    cleanup
    val=`cat $output_dir/$sample_id.out | grep "% of reads unmapped: too short"`
    descrip="unmapped_tooshort_%"
    cleanup
    val=`cat $output_dir/$sample_id.out | grep "% of reads unmapped: other"`
    descrip="unmapped_other_%"
    cleanup
}

function cleanup ()
{
    local cleanup=`echo $val | cut -f2 -d"|" | sed -s "s/\t//g" | sed -s "s/%//g" | sed -s "s/ //g"`
    echo "$cleanup,$descrip,$sample_id,$cpu,$mem,$seedPerWindowNmax,$winAnchorMultimapNmax" >> stats.csv

}

# create output
touch stats.csv
echo "val,descrip,sample_id,cpu,mem,seed,anchor" > stats.csv
################################################################################################
project_id="mESC_clip_4_v2.0"
################################################################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="5000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="100"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="50"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="10"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="5"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
# stats - failed
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="100"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="Control1hr_Clip"
seedPerWindowNmax="50"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
stats
################################################
sample_id="YKO7hr_Clip"
seedPerWindowNmax="5000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="YKO7hr_Clip"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="YKO7hr_Clip"
seedPerWindowNmax="50"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="YKO7hr_Clip"
seedPerWindowNmax="10"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Ro7hr2_Clip"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Ro7hr2_Clip"
seedPerWindowNmax="5000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Ro7hr2_Clip"
seedPerWindowNmax="50"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Ro7hr2_Clip"
seedPerWindowNmax="10"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Ro1hrNuc_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro1hrCyt_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro7hrNuc_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro7hrCyt_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro1hrNuc_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro1hrCyt_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro7hrNuc_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################
sample_id="Ro7hrCyt_Clip"
seedPerWindowNmax="500"
winAnchorMultimapNmax="50"
cpu="32"
mem="75g"

set_variables
create_sh
# submit_sh
################################################################################################
project_id="mESC_clip_2_v2.0"
################################################################################################
sample_id="Y_Clip_2"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Y_Clip_2"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################################################################
project_id="8-09-21-HaCaT_fCLIP_v2.0"
################################################################################################
sample_id="Y5KO_fCLIP"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Y5KO_fCLIP"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
#create_sh - fixed index for hg38
#submit_sh
################################################################################################
project_id="mESC_clip_3_clash_v2.0"
################################################################################################
sample_id="Ro_clash_uc"
seedPerWindowNmax="1000"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
################################################
sample_id="Ro_clash_uc"
seedPerWindowNmax="500"
winAnchorMultimapNmax="10000"
cpu="32"
mem="75g"

set_variables
create_sh
#submit_sh
