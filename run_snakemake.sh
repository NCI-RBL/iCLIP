test=$1

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

# read yaml file
eval $(parse_yaml config/snakemake_config.yaml "config_")

if [ -d "${config_out_dir}/log" ]
then
  echo "Pipeline re-run"
else
  mkdir "${config_out_dir}/log"
fi

#cp "config/snakemake_config.yaml ${config_out_dir}/log/snakemake_config.yaml"
#cp "config/cluster_config.yml ${config_out_dir}/log/cluster_config.yml"
#cp "${config_sample_manifest} ${config_out_dir}/log/sample_manifest.tsv"
#cp "${config_multiplex_manifest} ${config_out_dir}/log/multiplex_manifest.tsv"

if [[ $test = "run" ]]; then
  sbatch --job-name="iCLIP" --gres=lscratch:200 --time=120:00:00 --mail-type=BEGIN,END,FAIL \
  snakemake --latency-wait 120  -s workflow/Snakefile --printshellcmds --cluster-config config/cluster_config.yml --keep-going \
  --restart-times 1 --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} \
  --mem {cluster.mem} --cores {cluster.cores} --job-name={params.rname} --output=${config_out_dir}/log/{params.rname}.out" -j 500 --rerun-incomplete

elif [[ $test = "unlock" ]]; then
  snakemake -s workflow/Snakefile --unlock --cores=8
else
  #run snakemake
  snakemake -s workflow/Snakefile -npr
fi
