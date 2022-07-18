# 5. Pipeline Tutorial
Welcome to the iCLIP Pipeline Tutorial!

## 5.1 Getting Started
Review the information on the [Getting Started](https://rbl-nci.github.io/iCLIP/iCLIP/getting-started/) for a complete overview the pipeline. The tutorial below will use test data available on NIH Biowulf HPC only.

A. Change working directory to the iCLIP repository
```
cd /data/RBL_NCI/Pipelines/iCLIP/[version number]
```

B. Initialize Pipeline
```
sh run_snakemake.sh -p initialize -o /path/to/output/dir
```

## 5.2 Prepare the test set

A. Four different test data sets are available, depending on the need. These include:

- test_1: Single test (multiplex_flag="N", splice_aware="N", DE_method="none")
- test_2: Multiplexed test (multiplex_flag="Y", splice_aware="Y", DE_method="none")
- test_3: MANORM test (multiplex_flag="N", splice_aware="Y", DE_method="MANORM")
- test_4: DIFFBIND test (multiplex_flag="N", splice_aware="Y", DE_method="DIFFBIND")

B. Pull the test data to your output directory
```
# example
sh /data/CCBR_Pipeliner/iCLIP/test/run_test.sh -t TESTNAME -o /path/to/output/dir

# example running test_3:
sh /data/CCBR_Pipeliner/iCLIP/test/run_test.sh -t test_3 -o /path/to/output/dir -s /data/RBL_NCI/Pipelines/iCLIP/[version number]
```

## 5.3 Complete dry-run

A. Complete a dry-run and review output
```
sh run_snakemake.sh -p dry -o /path/to/output/dir/
```

Ensure that an expected output is displayed. An expected output for test_3 is as follows:
```
job                       count    min threads    max threads
----------------------  -------  -------------  -------------
MANORM_RMD                    2              2              2
MANORM_analysis               4              4              4
MANORM_beds                   4              4              4
MANORM_post_processing        2              2              2
all                           1              1              1
annotation_report             4              1              1
bgzip_beds                    4              4              4
create_beds_safs              4              8              8
dedup                         4              8              8
feature_counts                4              8              8
index_stats                   4              8              8
multiqc                       1              1              1
nondemux                      4              1              1
peak_ExonIntron               4             32             32
peak_RMSK                     8             32             32
peak_Transcripts              8             32             32
peak_junctions                4             32             32
peak_process                  4             32             32
project_annotations           1              1              1
qc_fastq                      4              1              1
qc_screen_validator           4             32             32
qc_troubleshoot               1              1              1
rename_fastqs                 1              1              1
star                          4             32             32
total                        85              1             32
```

## 5.4 Run the pipeline
Execute pipeline on the cluster OR locally
```
#submit to the cluster (recommended)
sh run_snakemake.sh -p cluster -o /path/to/output/dir/
```

## 5.5 Review outputs
Review the expected outputs on the [Output](https://rbl-nci.github.io/iCLIP/iCLIP/output/) page. If there are errors, review and performing stesp described on the [Troubleshooting](https://rbl-nci.github.io/iCLIP/iCLIP/troubleshooting/) page as needed.