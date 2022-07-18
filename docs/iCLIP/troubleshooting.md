# Troubleshooting
Recommended steps to troubleshoot the pipeline.

## 1.1 Email
Check your email for an email regarding pipeline failure. You will receive an email from slurm@biowulf.nih.gov with the subject: Slurm Job_id=[#] Name=iCLIP Failed, Run time [time], FAILED, ExitCode 1

## 1.2 Error Report
Run the error report script
```
cd /[output_dir]/log/[time_of_run]
sh 00_create_error_report.sh
cat error.log
```

Review the report for the rules that erred, and the sample information. An example report is listed below:
```
The following error(s) were found in rules:
*********************************************
Error in rule rule1:
Error in rule rule2:
Error in rule rule3:

The following samples are affected by memory and must be deleted:
rule1.[sbatchid].sp=[sample_name].err:[E::hts_open_format] Disc quota exceeded

The following samples are affected by missing input files/output dir and should be reviewed:
rule2.[sbatchid].sp=[sample_name].err:[E::hts_open_format] Failed to open file "[file_name]" : No such file or directory

The following samples are affected by other error_rules and should be reviewed:
rule3.[sbatchid].sp=[sample_name].err:[E::hts_open_format] TIMEOUT
```

## 1.3 Restart the run
After addressing the issue, unlock the output directory, perform another dry-run and check the status of the pipeline, then resubmit to the cluster.
```
#unlock dir
sh run_snakemake.sh -p unlock -o /path/to/output/dir

#perform dry-run
sh run_snakemake.sh -p dry -o /path/to/output/dir

#submit to cluster
sh run_snakemake.sh -p cluster -o /path/to/output/dir
```