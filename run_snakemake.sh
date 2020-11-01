test=$1

if [ $test = "run" ]; then
  snakemake -s workflow/Snakefile --cores=8 --latency-wait=60 --printshellcmds
elif [[ $test = "unlock" ]]; then
  snakemake -s workflow/Snakefile --unlock --cores=8
else
  #cleanup output for testing
  rm -r /data/sevillas2/iCLIP/output/*.txt
  rm /home/sevillas2/git/iCLIP/iCount.log

  #run snakemake
  snakemake -s workflow/Snakefile -npr
fi
