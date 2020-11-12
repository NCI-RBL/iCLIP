test=$1

if [[ $test = "run" ]]; then
  snakemake -s workflow/Snakefile --cores=8 --latency-wait=60 --printshellcmds
elif [[ $test = "unlock" ]]; then
  snakemake -s workflow/Snakefile --unlock --cores=8
else
  #run snakemake
  snakemake -s workflow/Snakefile -npr
fi
