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
   echo -e "\t<PIPELINEMODE> options: initialize, dry, cluster, local, git, unlock, DAG, report, check"
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
load_modules(){
  if [[ $1 =~ "python" ]]; then module load python/3.7; fi
  if [[ $1 =~ "snakemake" ]]; then module load snakemake/7.19.1; fi
  if [[ $1 =~ "graphviz" ]]; then module load graphviz/2.40; fi
}

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

# Read in multiplex manifest, sample manifest, contrast manifest.
# Use python script to check file matching, sample matching, and invalid characters.
# Will only check contrast manifest if "MANORM" is selected
# If files are correct, outputs: 
# 1) an empty "no_error" file
# 2) barcode-to-sample file needed for barcode/adaptor removal
# If there are errors, outputs:
# 1) a detailed error file
check_manifests(){

  load_modules "python"
  
  # set args
  py_script_in="$1"
  manifest_file_prefix_in="$2"

  # if the third argument is missing, it is a standard run
  # if the third argument is present, but the fifth is missing it's a test run - NOT for DE
  # if the third and fifth arguments are given, then its a test run - for DE
  # for manifest testing
  if [[ -z "$3" ]]; then
    multiplexManifest=$config_multiplexManifest
    sampleManifest=$config_sampleManifest
    DEmethod=`echo "${config_DEmethod}" | cut -f1 -d"#"`
    contrastManifest=$config_contrastManifest

    check_existence $multiplexManifest
    check_existence $sampleManifest
    check_existence $contrastManifest
  elif [[ -z "$5" ]]; then
    multiplexManifest=$3
    sampleManifest=$4
    DEmethod="NONE"
    contrastManifest=""

    check_existence $multiplexManifest
    check_existence $sampleManifest
  else
    multiplexManifest=$3
    sampleManifest=$4
    DEmethod=$5
    contrastManifest=$6

    check_existence $multiplexManifest
    check_existence $sampleManifest
    check_existence $contrastManifest
  fi
  
  # Run manifest check
  python $py_script_in \
    $manifest_file_prefix_in \
    ${multiplexManifest} \
    ${sampleManifest} \
    ${DEmethod} \
    ${contrastManifest}
}

# create ultraplex barcode manifest
create_barcode_manifest(){
    # create list of all multiplexed samples
    mp_id_list=`sed -n '1d;p' ${config_multiplexManifest} | cut -f2 -d"," | cut -f1 -d"." | uniq`
    # for each id create a single barcode manifest
    for mp_id in ${mp_id_list[@]}; do
      bc_file="${output_dir}/manifests/${mp_id}_barcode_manifest.txt"
      
      # if the file exists, skip
      if [[ ! -f ${output_dir}/manifests/${mp_id}_barcode_manifest.txt ]]; then
        cat ${config_sampleManifest} | grep "$mp_id" | awk -F"," '{ print $4":"$2 }' > $bc_file
      fi
    done
}

check_manifest_qc(){
  #check for errors in manifests (sample,multiplex,contrasts)
  if [[ ! -f "${manifest_file}no_errors.txt" ]]; then
    echo "The manifest check FAILED. Check "${output_dir}/qc/" for more information."
    exit 1
  else
    echo "-- user manifest check completed successfully"
  fi

  # check for errors in barcode manifest, if multiplex flag is turned on
  mp_flag=`echo "${config_multiplexflag}" | cut -f1 -d" "`
  if [[ ${mp_flag} == "Y" ]]; then
    # search for all barcodes
    mp_id_list=`sed -n '1d;p' ${config_multiplexManifest} | cut -f2 -d"," | cut -f1 -d"." | uniq`
    
    for mp_id in ${mp_id_list[@]}; do
        # check length of file
        bc_file="${output_dir}/manifests/${mp_id}_barcode_manifest.txt"
        len_bc=`cat $bc_file | wc -l`

        #if 0 then print error
        if [[ "$len_bc" -lt 1 ]]; then
          echo "The barcode manifest check FAILED. Check "${bc_file}" for more information."
          exit 1
        else
          echo "-- barcode manifest check completed successfully (sample $mp_id)"
        fi
    done
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

#########################################################
# Pipeline options
#########################################################
####################### INITIALIZE #######################
if [[ $pipeline = "initialize" ]]; then
  echo "------------------------------------------------------------------------"
	echo "*** STARTING Initialization ***"

  # create dirs
  if [[ ! -d "${output_dir}" ]]; then mkdir ${output_dir}; fi
  dir_list=(log config manifests qc)
  for pd in "${dir_list[@]}"; do if [[ ! -d $output_dir/$pd ]]; then mkdir -p $output_dir/$pd; fi; done
  
  # copy config inputs to edit
  files_save=('config/snakemake_config.yaml' 'config/cluster_config.yaml' 'config/index_config.yaml' 'config/annotation_config.txt' 'workflow/Snakefile')
  for f in ${files_save[@]}; do
    f="${PIPELINE_HOME}/$f"
    IFS='/' read -r -a strarr <<< "$f"
    sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/OUTPUT_DIR/${output_dir//\//\\/}/g" $f > "${output_dir}/config/${strarr[-1]}"
  done

  # copy example manifests
  files_save=('manifests/contrasts_example_diffbind.tsv' 'manifests/contrasts_example_manorm.tsv' 'manifests/samples_example.tsv' 'manifests/multiplex_example.tsv')
  for f in ${files_save[@]}; do
    f="${PIPELINE_HOME}/$f"
    IFS='/' read -r -a strarr <<< "$f"
    cp $f "${output_dir}/manifests/${strarr[-1]}"
  done

  echo "*** COMPLETE Initialization ***"
  echo "------------------------------------------------------------------------"

####################### CHECK, CLUSTER, LOCAL #######################
# Run unit checks
elif [[ $pipeline = "check" ]]; then
    echo "------------------------------------------------------------------------"
	  echo "*** STARTING CHECKS ***"

    # Run manifest checks
    manifest_file="${output_dir}/manifest_"
    py_script="${PIPELINE_HOME}/workflow/scripts/01_check_manifest.py"

    #set files
    date_stamp=`date +'%Y%m%d'`
    manifest_log=${output_dir}/manifest_log.txt
    pass_qc="${output_dir}/manifest_no_errors.txt"
    manorm_qc="${output_dir}/manifest_MANORM_comparisons.txt"
    diffbind_qc="${output_dir}/manifest_DIFFBIND_comparisons.txt"
    fail_qc="${output_dir}/manifest_contains_errors_${date_stamp}.txt"

    # Run manifest check 
    ## Test 1-2 expected no errors
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_two.tsv" "MANORM" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/contrasts_manorm.tsv"
    if [[ -f $pass_qc ]] && [[ -f $manorm_qc ]]; then rm $pass_qc; rm $manorm_qc; echo "Test1 pass" > $manifest_log; else echo "Test1 fail" > $manifest_log; fi

    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_two.tsv" "DIFFBIND" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/contrasts_diffbind.tsv"
    if [[ -f $pass_qc ]] && [[ -f $diffbind_qc ]]; then rm $pass_qc; rm $diffbind_qc; echo "Test2 pass" >> $manifest_log; else echo "Test2 fail" > $manifest_log; fi

    ## Test 3 - expect error with duplicate values
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_sampledups.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "sample names must be unique" | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test3 pass" >> $manifest_log; else echo "Test3 fail" >> $manifest_log; fi
    else
      echo "Test3 fail" >> $manifest_log
    fi

    ## Test 4 - barcode has duplicates
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_barcodedups.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "Barcodes must be unique by sample" | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test4 pass" >> $manifest_log; else echo "Test4 fail" >> $manifest_log; fi
    else
      echo "Test4 fail" >> $manifest_log
    fi

    ## Test 5 - FASTQ files are duplicated
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_filedups.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_two.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "File names must be unique" | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test5 pass" >> $manifest_log; else echo "Test5 fail" >> $manifest_log; fi
    else
      echo "Test5 fail" >> $manifest_log
    fi

    ## Test 6 - FASTQ files wrong extension
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_fileext.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_two.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "All items in file_name column must end in fastq.gz" | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test6 pass" >> $manifest_log; else echo "Test6 fail" >> $manifest_log; fi
    else
      echo "Test6 fail" >> $manifest_log
    fi

    ## Test 7 - sample names don't match between multiplex and sample
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_extrasamples.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "Multiplex ID's must be consistent between both files." | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test7 pass" >> $manifest_log; else echo "Test7 fail" >> $manifest_log; fi
    else
      echo "Test7 fail" >> $manifest_log
    fi

    ## Test 8-9 - DE method requests do not match sample names / group names
    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_two.tsv" "MANORM" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/contrasts_manormfail.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "was/were not found in the sample_manifest tsv but were found in the contrasts manifest" | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test8 pass" >> $manifest_log; else echo "Test8 fail" >> $manifest_log; fi
    else
      echo "Test8 fail" >> $manifest_log
    fi

    check_manifests $py_script $manifest_file \
      "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/multiplex_two.tsv" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/samples_two.tsv" "DIFFBIND" "/data/CCBR_Pipeliner/iCLIP/test/pipeline_checks/manifests/contrasts_diffbindfail.tsv"
    if [[ -f "${fail_qc}" ]]; then 
      check=`cat ${fail_qc} | grep "was/were not found in the sample_manifest tsv but were found in the contrasts manifest" | wc -l`
      if [[ $check > 0 ]]; then rm $fail_qc; echo "Test9 pass" >> $manifest_log; else echo "Test9 fail" >> $manifest_log; fi
    else
      echo "Test9 fail" >> $manifest_log
    fi
    
    cat $manifest_log

 
# Run pipeline locally/cluster
elif [[ $pipeline = "cluster" ]] || [[ $pipeline = "local" ]]; then
  echo "------------------------------------------------------------------------"
	echo "*** STARTING PIPELINE ***"

  ####################### Preparation
  #parse config
  eval $(parse_yaml ${output_dir}/config/snakemake_config.yaml "config_")  
  eval $(parse_yaml ${output_dir}/config/index_config.yaml "yaml_")
  
  #save source dir
  source_dir=$(echo $config_sourceDir | sed 's:/*$::')

  # prep manifest check
  manifest_file="${output_dir}/qc/manifest_"
  py_script="${config_sourceDir}/workflow/scripts/01_check_manifest.py"

  #run checks
  check_initialization
  check_output_dir
  check_manifests $py_script $manifest_file
  create_barcode_manifest
  check_manifest_qc

  if [[ $config_reference == "hg38" ]]; then
    check_readaccess "${yaml_hg38_stargtf}"
    check_readaccess "${yaml_hg38_stardir}"
    check_readaccess "${yaml_hg38_gencodepath}"
    check_readaccess "${yaml_hg38_refseqpath}"
    check_readaccess "${yaml_hg38_canonicalpath}"
    check_readaccess "${yaml_hg38_intronpath}"
    check_readaccess "${yaml_hg38_rmskpath}"
    check_readaccess "${yaml_hg38_sypath}"
    check_readaccess "${yaml_hg38_aliaspath}"
  else
    check_readaccess "${yaml_mm10_stargtf}"
    check_readaccess "${yaml_mm10_stardir}"
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
  files_save=("${output_dir}/config/snakemake_config.yaml" "${output_dir}/config/cluster_config.yaml" "${output_dir}/config/index_config.yaml" "${output_dir}/config/Snakefile" ${config_multiplexManifest} ${config_sampleManifest}  ${config_contrastManifest} "${PIPELINE_HOME}/workflow/scripts/create_error_report.sh")

  for f in ${files_save[@]}; do
    IFS='/' read -r -a strarr <<< "$f"
    cp $f "${output_dir}/log/${log_time}/00_${strarr[-1]}"
  done

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"

  ####################### Selection
  if [[ $pipeline = "cluster" ]]; then
    echo "-- submitting to cluster"
    echo "---- output dir: ${output_dir}"
    echo "---- pipeline jobid:"
    
  cat > ${output_dir}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name=iCLIP_v2.2
#SBATCH --mem=40g
#SBATCH --gres=lscratch:200
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=${output_dir}/log/${log_time}/00_%j_%x.out \
#SBATCH --mail-type=BEGIN,END,FAIL
module load snakemake
module load graphviz
cd \$SLURM_SUBMIT_DIR
    snakemake \
    --scheduler greedy \
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
    --job-name={params.rname} --output=${output_dir}/log/${log_time}/{params.rname}{cluster.output} \
    --error=${output_dir}/log/${log_time}/{params.rname}{cluster.error}" \
    2>&1|tee ${output_dir}/log/${log_time}/snakemake.log

if [ "\$?" -eq "0" ];then
  snakemake -s ${output_dir}/log/${log_time}/00_Snakefile \
  --directory $output_dir \
  --report ${output_dir}/log/${log_time}/runslurm_snakemake_report.html \
  --configfile ${output_dir}/log/${log_time}/00_snakemake_config.yaml 
fi

#bash <(curl https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/gather_cluster_stats_biowulf.sh 2>/dev/null) ${output_dir}/log/${log_time}/snakemake.log > ${output_dir}/log/${log_time}/snakemake.log.HPC_summary.txt
#${PIPELINE_HOME}/workflow/scripts/jobinfo -s ${output_dir}/log/${log_time}/snakemake.log -o ${output_dir}/log/${log_time}/snakemake.log.HPC_summary.txt

EOF

sbatch ${output_dir}/submit_script.sbatch

  # run local - includes all rules - run locally
  else
    #remove iCount dir if it already exist - will cause error in demux
    if [ -d "/tmp/iCount" ]; then rm -r "/tmp/iCount/"; fi
    
    load_modules "snakemake"

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
  echo "------------------------------------------------------------------------"
	echo "*** STARTING Unlock ***"
  load_modules "snakemake"

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
  echo "------------------------------------------------------------------------"
	echo "*** STARTING GITTests ***"

  snakemake \
  -s "${output_dir}/config/Snakefile" \
  --configfile .tests/snakemake_config.yaml \
  --printshellcmds \
  --cluster-config .tests/cluster_config.yaml \
  -npr
######################## DAG #######################
elif [[ $pipeline = "DAG" ]]; then
  echo "------------------------------------------------------------------------"
	echo "*** STARTING DAG ***"
  load_modules "snakemake graphviz"

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  snakemake \
  -s "${output_dir}/config/Snakefile" \
  --configfile ${output_dir}/config/snakemake_config.yaml \
  --rulegraph | dot -Tpdf > ${output_dir}/log/dag.pdf

  echo "DAG is available at ${output_dir}/log/dag.pdf"
######################## Report #######################
elif [[ $pipeline = "report" ]]; then
  echo "------------------------------------------------------------------------"
	echo "*** STARTING Snakemake Report ***"

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
  echo "------------------------------------------------------------------------"
	echo "*** STARTING Cleanup ***"
  
  load_modules "snakemake"

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
  echo "------------------------------------------------------------------------"
	echo "*** STARTING DryRun ***"
  
  load_modules "snakemake"

  ####################### Preparation
  #parse config
  eval $(parse_yaml ${output_dir}/config/snakemake_config.yaml "config_")  

  # prep manifest check
  manifest_file="${output_dir}/qc/manifest_"
  py_script="${config_sourceDir}/workflow/scripts/01_check_manifest.py"

  #run check
  check_initialization
  check_output_dir
  check_manifests $py_script $manifest_file
  create_barcode_manifest
  check_manifest_qc

  # Change directories before running any important
  # snakemake commands, ensures the .snakemake 
  # directory get created in the output directory
  cd "${output_dir}"
  if [[ ! -d "${output_dir}/log/dryrun" ]]; then mkdir "${output_dir}/log/dryrun"; fi

  snakemake -s ${output_dir}/config/Snakefile \
  --configfile ${output_dir}/config/snakemake_config.yaml \
  --printshellcmds \
  --verbose \
  --rerun-incomplete \
  --scheduler greedy  \
  --rerun-triggers mtime \
  --cluster-config ${output_dir}/config/cluster_config.yaml \
  --cluster \
    "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} \
    -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} \
    --job-name={params.rname} --output=${output_dir}/log/${log_time}/{params.rname}{cluster.output} --error=${output_dir}/log/${log_time}/{params.rname}{cluster.error}" \
  --jobs 100 \
  -npr | tee ${output_dir}/log/dryrun/dryrun.${log_time}.log
fi
echo
