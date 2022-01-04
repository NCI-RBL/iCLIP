#!/usr/bin/env bash

#########################################################
# Arguments
#########################################################

helpFunction()
{
   echo "#########################################################" 
   echo "Usage: bash $0 -p <PIPELINEMODE> -o <OUTPUTDIR>"
   echo "#########################################################" 
   echo "Acceptable inputs:"
   echo -e "\t<PIPELINEMODE> options: initialize, check, dry, cluster, local, git, unlock, DAG, report"
   echo -e "\t<OUTPUTDIR> : absolute path to output folder required"
   echo "#########################################################" 
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

#########################################################
# functions
#########################################################

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

check_initialization(){
  if [[ ! -d $output_dir ]] || [[ ! -d "${output_dir}/log" ]]; then 
    echo "ERROR: You must initalize the dir before beginning pipeline"
    exit 1
  fi
}

check_output_dir(){
  eval $(parse_yaml ${output_dir}/config/snakemake_config.yaml "config_")
  config_outputDir=$(echo $config_outputDir | sed 's:/*$::')

  if [[ ! ${config_outputDir} == $output_dir ]]; then 
    echo "ERROR: Output dir provided: $output_dir does not match snakemake_config: ${config_outputDir}. Update and re-run."
    exit 1
  fi
}    

check_existence(){
 if [[ ! -f "$1" ]]; then
    echo "ERROR: File does not exist: $1"
    exit 1
  fi
}

check_readaccess(){
  check_existence()
  
  if [[ -r "$1" ]]; then
    echo "ERROR: User does does not have read access to: $1"
    exit 1
  fi
}
  
check_writeaccess(){
  if [[ -w "$1" ]]; then
    echo "ERROR: User does does not have write access to: $1"
    exit 1
  fi
}

#########################################################  
# Formatting
#########################################################
log_time=`date +"%Y%m%d_%H%M"`
s_time=`date +"%Y%m%d_%H%M%S"`

#remove trailing / on directories
output_dir=$(echo $output_dir | sed 's:/*$::')

# ## setting PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))
echo "Pipeline Dir: $PIPELINE_HOME"

#########################################################
# Pipeline options
#########################################################
####################### INITIALIZE #######################
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

    if [ ! -d "${output_dir}/config" ]; then
      mkdir "${output_dir}/config/" 
    fi

    if [ ! -d "${output_dir}/manifest" ]; then
      mkdir "${output_dir}/manifest/" 
    fi

  else
    # Create any initialized parent directories
    mkdir -p "${output_dir}"
    mkdir "${output_dir}/config"
    mkdir "${output_dir}/manifest"
    mkdir "${output_dir}/log"
    echo
    echo "Creating output dir with configs: ${output_dir}"
  fi
  
  # copy config inputs to edit
  files_save=('config/snakemake_config.yaml' 'config/cluster_config.yaml' 'config/index_config.yaml' 'config/annotation_config.txt' 'workflow/Snakefile')

  for f in ${files_save[@]}; do
  # converting $f path to absolute path so that
  # run_snakemake.sh can be run from any location on the file system
  # not just from the PIPELINE_HOME folder
    f="${PIPELINE_HOME}/$f"
    IFS='/' read -r -a strarr <<< "$f"
    sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/OUTPUT_DIR/${output_dir//\//\\/}/g" $f > "${output_dir}/config/${strarr[-1]}"
  done

  # copy example manifests
  files_save=('manifests/contrasts_example.tsv' 'manifests/samples_example.tsv' 'manifests/multiplex_example.tsv')

  for f in ${files_save[@]}; do
  # converting $f path to absolute path so that
  # run_snakemake.sh can be run from any location on the file system
  # not just from the PIPELINE_HOME folder
    f="${PIPELINE_HOME}/$f"
    IFS='/' read -r -a strarr <<< "$f"
    cp $f "${output_dir}/manifest/${strarr[-1]}"
  done
####################### CHECK, CLUSTER, LOCAL #######################
#Run check of pipeline OR run pipeline locally/cluster
elif [[ $pipeline = "check" ]] || [[ $pipeline = "cluster" ]] || [[ $pipeline = "local" ]]; then
  echo
  echo "Running pipeline"

  ####################### Preparation
  #parse config
  eval $(parse_yaml ${output_dir}/config/snakemake_config.yaml "config_")  
  eval $(parse_yaml ${output_dir}/config/index_config.yaml "yaml_")
  
  #save source dir
  source_dir=$(echo $config_sourceDir | sed 's:/*$::')

  #run checks
  check_initialization
  check_output_dir

  if [[ $config_reference == "hg38" ]]; then
    check_readaccess "${yaml_hg38_std}"
    check_readaccess "${yaml_hg38_spliceawareunmasked_50bp}"
    check_readaccess "${yaml_hg38_spliceawareunmasked_75bp}"
    check_readaccess "${yaml_hg38_spliceawaremasked_50bp}"
    check_readaccess "${yaml_hg38_spliceawaremasked_75bp}"
    check_readaccess "${yaml_hg38_gencodepath}"
    check_readaccess "${yaml_hg38_refseqpath}"
    check_readaccess "${yaml_hg38_canonicalpath}"
    check_readaccess "${yaml_hg38_intronpath}"
    check_readaccess "${yaml_hg38_rmskpath}"
    check_readaccess "${yaml_hg38_sypath}"
    check_readaccess "${yaml_hg38_aliaspath}"
  else
    check_readaccess "${yaml_mm10_std}"
    check_readaccess "${yaml_mm10_spliceawareunmasked_50bp}"
    check_readaccess "${yaml_mm10_spliceawareunmasked_75bp}"
    check_readaccess "${yaml_mm10_spliceawaremasked_50bp}"
    check_readaccess "${yaml_mm10_spliceawaremasked_75bp}"
    check_readaccess "${yaml_mm10_gencodepath}"
    check_readaccess "${yaml_mm10_refseqpath}"
    check_readaccess "${yaml_mm10_canonicalpath}"
    check_readaccess "${yaml_mm10_intronpath}"
    check_readaccess "${yaml_mm10_rmskpath}"
    check_readaccess "${yaml_mm10_sypath}"
    check_readaccess "${yaml_mm10_aliaspath}"
  fi

  #create run log dir
  mkdir "${output_dir}/log/${log_time}"
  
  # copy config inputs for ref
  files_save=("${output_dir}/config/snakemake_config.yaml" "${output_dir}/config/cluster_config.yaml" "${output_dir}/config/index_config.yaml" "${output_dir}/config/Snakefile" ${config_multiplexManifest} ${config_sampleManifest} "${PIPELINE_HOME}/workflow/scripts/create_error_report.sh")

  for f in ${files_save[@]}; do
    IFS='/' read -r -a strarr <<< "$f"
    cp $f "${output_dir}/log/${log_time}/00_${strarr[-1]}"
  done

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  ####################### Selection
  # run check - only includes check_manifest rule - run locally
  if [[ $pipeline = "check" ]]; then
    echo
    echo "Running manifest check"

    module load snakemake
    snakemake \
    -s ${output_dir}/log/${log_time}/00_Snakefile \
    --use-envmodules \
    --configfile ${output_dir}/log/${log_time}/00_snakemake_config.yaml \
    --printshellcmds \
    --cluster-config ${output_dir}/log/${log_time}/00_cluster_config.yaml \
    --cores 8 \
    --until check_manifest

    echo "If running differential expression, review the qc/manifest_check.txt file to confirm your settings"

  # run cluster - includes all rules - run on cluster
  elif [[ $pipeline = "cluster" ]]; then
    echo
    echo "Output dir: ${output_dir}"
    echo "Pipeline jobid:"
    module load snakemake

    sbatch \
    --job-name="iCLIP" \
    --gres=lscratch:200 \
    --time=10-00:00:00 \
    --output=${output_dir}/log/${log_time}/00_%j_%x.out \
    --mail-type=BEGIN,END,FAIL \
    snakemake \
    --use-envmodules \
    --latency-wait 120 \
    -s ${output_dir}/log/${log_time}/00_Snakefile \
    --configfile ${output_dir}/log/${log_time}/00_snakemake_config.yaml \
    --printshellcmds \
    --cluster-config ${output_dir}/log/${log_time}/00_cluster_config.yaml \
    --keep-going \
    --restart-times 1 \
    -j 500 \
    --rerun-incomplete \
    --stats ${output_dir}/log/${log_time}/snakemake.stats \
    --cluster \
    "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} \
    -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} \
    --job-name={params.rname} --output=${output_dir}/log/${log_time}/{params.rname}{cluster.output} --error=${output_dir}/log/${log_time}/{params.rname}{cluster.error}" \
    |tee ${output_dir}/log/${log_time}/snakemake.log

  # run local - includes all rules - run locally
  else
    #remove iCount dir if it already exist - will cause error in demux
    if [ -d "/tmp/iCount" ]; then rm -r "/tmp/iCount/"; fi
    module load snakemake

    snakemake \
    -s ${output_dir}/log/${log_time}/00_Snakefile \
    --use-envmodules \
    --configfile ${output_dir}/log/${log_time}/00_snakemake_config.yaml \
    --printshellcmds \
    --cluster-config ${output_dir}/log/${log_time}/00_cluster_config.yaml \
    --cores 8 \
    --stats ${output_dir}/log/${log_time}/snakemake.stats
  fi

####################### UNLOCK #######################
elif [[ $pipeline = "unlock" ]]; then
  echo
  echo "Unlocking pipeline"
  module load snakemake

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  snakemake \
  -s "${output_dir}/config/Snakefile" \
  --use-envmodules \
  --unlock \
  --cores 8 \
  --configfile ${output_dir}/config/snakemake_config.yaml
######################## GIT #######################
#Run github actions
elif [[ $pipeline = "git" ]]; then
  echo
  echo "Starting Git Tests"

  snakemake \
  -s "${output_dir}/config/Snakefile" \
  --configfile .tests/snakemake_config.yaml \
  --printshellcmds \
  --cluster-config .tests/cluster_config.yaml \
  -npr
######################## DAG #######################
elif [[ $pipeline = "DAG" ]]; then
  module load snakemake

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  snakemake \
  -s "${output_dir}/config/Snakefile" \
  --configfile ${output_dir}/config/snakemake_config.yaml \
  --rulegraph | dot -Tpdf > ${output_dir}/log/dag.pdf
######################## Report #######################
elif [[ $pipeline = "report" ]]; then
  echo
  echo "Generating Snakemake Report Command"
  echo

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  report_cmd="module load R; Rscript -e 'library(rmarkdown);
  rmarkdown::render(\"${PIPELINE_HOME}/workflow/scripts/create_snakemake_report.RMD\",
  output_file = \"${output_dir}/snakemake_report.html\", 
  params= (log_dir=\"${output_dir}/\"))'"
  echo $report_cmd

######################## cleanup #######################
elif [[ $pipeline = "cleanup" ]]; then
  echo
  echo "Starting cleanup"
  
  module load snakemake

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  snakemake \
  --cleanup-metadata $2 \
  -s ${output_dir}/config/Snakefile \
  --configfile ${output_dir}/config/snakemake_config.yaml \
  --cores 1 \
  --cleanup-shadow

######################## Dry #######################
else
  echo
  echo "Starting dry-run"
  
  #run check
  check_initialization
  check_output_dir
  
  module load snakemake

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  snakemake -s ${output_dir}/config/Snakefile \
  --configfile ${output_dir}/config/snakemake_config.yaml \
  --printshellcmds \
  --cluster-config ${output_dir}/config/cluster_config.yaml \
  -npr | tee ${output_dir}/dryrun.${log_time}.log
fi
echo
