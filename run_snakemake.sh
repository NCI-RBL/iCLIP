#!/bin/bash

##############################################################################
#handle arguments

helpFunction()
{
   echo ""
   echo "Usage: $0 -p pipeline"
   echo -e "\t-p options: initialize, cluster, local, dry-run, unlock, report"
   echo "Usage: $1 -o output_dir"
   echo -e "\t-o path to output directory"
   exit 1 # Exit script after printing help
}

while getopts "p:o:" opt
do
   case "$opt" in
      p ) pipeline="$OPTARG" ;;
      o ) output_dir="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$pipeline" ] || [ -z "$output_dir" ]; then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

#####
#handle yaml file
parse_yaml() {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

##########
#Run check that output dir in config matches output dir on cmd
check_initialization(){
  if [[ ! -d $output_dir ]] || [[ ! -d "${output_dir}/log" ]]; then 
    echo "ERROR: You must initalize the dir before beginning pipeline"
    exit 1
  fi
}
check_output_dir(){
  eval $(parse_yaml ${output_dir}/snakemake_config.yaml "config_")
  config_output_dir=$(echo $config_output_dir | sed 's:/*$::')

  if [[ ! ${config_output_dir} == $output_dir ]]; then 
    echo "ERROR: Output dir provided: $output_dir does not match snakemake_config: ${config_output_dir}. Update and re-run"
    exit 1
  fi
}    
  
# set timestamp
log_time=`date +"%Y%m%d_%H%M"`
s_time=`date +"%Y%m%d_%H%M%S"`

#remove trailing / on directories
output_dir=$(echo $output_dir | sed 's:/*$::')

#Run initialization step
if [[ $pipeline = "initialize" ]]; then
  echo
  echo "Initializing pipeline"

  #check output dir, log dir
  if [ -d "${output_dir}" ]; then
    echo
    echo "Output dir: ${output_dir}"

    if [ ! -d "${output_dir}/log" ]; then
      mkdir "${output_dir}/log/" 
    fi
  else
    mkdir "${output_dir}"
    mkdir "${output_dir}/log"
    echo
    echo "Creating output dir with configs: ${output_dir}"
  fi
  
  # copy config inputs to edit
  files_save=('config/snakemake_config.yaml' 'config/cluster_config.yml' 'config/index_config.yaml')

  for f in ${files_save[@]}; do
    IFS='/' read -r -a strarr <<< "$f"
    cp $f "${output_dir}/${strarr[-1]}"
  done
#Run pipeline on cluster or locally
elif [[ $pipeline = "cluster" ]] || [[ $pipeline = "local" ]]; then
  echo
  echo "Running pipeline"

  #run checks
  check_initialization
  check_output_dir

  #parse config
  eval $(parse_yaml ${output_dir}/snakemake_config.yaml "config_")

  #remove trailing / on directories
  source_dir=$(echo $config_source_dir | sed 's:/*$::')

  #create run log dir
  mkdir "${output_dir}/log/${log_time}"
  
  # copy config inputs for ref
  files_save=('${output_dir}/snakemake_config.yaml' '${output_dir}/cluster_config.yml' '${output_dir}/index_config.yaml' ${config_multiplex_manifest} ${config_sample_manifest} 'workflow/Snakefile' 'workflow/scripts/create_error_report.sh')

  for f in ${files_save[@]}; do
    IFS='/' read -r -a strarr <<< "$f"
    cp $f "${output_dir}/log/${log_time}/00_${strarr[-1]}"
  done

  #submit jobs to cluster
  if [[ $pipeline = "cluster" ]]; then
    echo
    echo "Pipeline jobid:"
    
    sbatch \
    --job-name="iCLIP" \
    --gres=lscratch:200 \
    --time=24:00:00 \
    --output=${output_dir}/log/${log_time}/00_%j_%x.out \
    --mail-type=BEGIN,END,FAIL \
    snakemake \
    --use-envmodules \
    --latency-wait 120 \
    -s ${output_dir}/log/${log_time}/00_Snakefile \
    --configfile ${output_dir}/log/${log_time}/00_snakemake_config.yaml \
    --printshellcmds \
    --cluster-config ${output_dir}/log/${log_time}/00_cluster_config.yml \
    --keep-going \
    --restart-times 1 \
    -j 500 \
    --rerun-incomplete \
    --stats ${output_dir}/log/${log_time}/snakemake.stats \
    --cluster \
    "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} \
    -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} \
    --job-name={params.rname} --output=${output_dir}/log/${log_time}/{params.rname}.out"
  #submit jobs locally
  else
    #remove iCount dir if it already exist - will cause error in demux
    if [ -d "/tmp/iCount" ]; then rm -r /tmp/iCount/; fi
    
    snakemake \
    -s ${output_dir}/log/${log_time}/00_Snakefile \
    --use-envmodules \
    --configfile ${output_dir}/log/${log_time}/00_snakemake_config.yaml \
    --printshellcmds \
    --cluster-config ${output_dir}/log/${log_time}/00_cluster_config.yml \
    --cores 8 \
    --stats ${output_dir}/log/${log_time}/snakemake.stats
  fi
#Unlock pipeline
elif [[ $pipeline = "unlock" ]]; then
  echo
  echo "Unlocking pipeline"
  snakemake \
  -s workflow/Snakefile \
  --use-envmodules \
  --unlock \
  --cores 8 \
  --configfile ${output_dir}/snakemake_config.yaml
#Run github actions
elif [[ $pipeline = "test" ]]; then
  snakemake \
  -s workflow/Snakefile \
  --configfile .tests/snakemake_config.yaml \
  --printshellcmds \
  --cluster-config ${output_dir}/cluster_config.yml \
  -npr
#Create DAG
elif [[ $pipeline = "DAG" ]]; then
  snakemake \
  -s workflow/Snakefile \
  --configfile .tests/snakemake_config.yaml \
  --rulegraph | dot -Tpdf > ${output_dir}/dag.pdf
elif [[ $pipeline = "report" ]]; then
  snakemake -s workflow/Snakfile \
  --report ${output_dir}/runlocal_snakemake_report.html \
  --directory $output_dir \
  --configfile ${output_dir}/snakemake_config.yaml 
#Dry-run pipeline
else
  echo
  echo "Starting dry-run"
  
  #run check
  check_initialization
  check_output_dir
  
  snakemake -s workflow/Snakefile \
  --configfile ${output_dir}/snakemake_config.yaml \
  --printshellcmds \
  --cluster-config ${output_dir}/cluster_config.yml \
  -npr
fi