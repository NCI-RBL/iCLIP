# 5. Pipeline Tutorial
Welcome to the iCLIP Pipeline Tutorial!

## 5.1 Getting Started
Review the information on the [Getting Started](https://rbl-nci.github.io/iCLIP/iCLIP/getting-started/) for a complete overview the pipeline. The tutorial below will use test data available on NIH Biowulf HPC only. All example code will assume you are running v2.2 of the pipeline, from the shared RBL_NCI storage directory, using test_1 data.

A. Change working directory to the iCLIP repository
```
# general format
cd /data/RBL_NCI/Pipelines/iCLIP/[version number]

# example
cd /data/RBL_NCI/Pipelines/iCLIP/v2.2
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
NOTE: Test data is currently available for v1.8, v2.0, v2.1, v2.2. Please contact samantha.sevilla@nih.gov to create other test data.

```
# general format
sh /data/CCBR_Pipeliner/iCLIP/test/run_test.sh \
    -t test_number \
    -v version_id \
    -s /path/to/source/dir
    -o /path/to/output/dir

# example running test_1, v2.2:
sh /data/CCBR_Pipeliner/iCLIP/test/run_test.sh \
    -t test_1 \
    -v v2.2 \
    -s /data/RBL_NCI/Pipelines/iCLIP/v2.2 \
    -o /path/to/output/dir 
```

## 5.3 Complete dry-run

A. Complete a dry-run and review output
```
sh run_snakemake.sh -p dry -o /path/to/output/dir/
```

Ensure that an expected output is displayed. 
- An expected output for test_1 is as follows:
```
job                    count    min threads    max threads
-------------------  -------  -------------  -------------
all                        1              1              1
annotation_report          1              1              1
bgzip_beds                 1              4              4
create_beds_safs           1              8              8
dedup                      1              8              8
feature_counts             1              8              8
index_stats                1              8              8
multiqc                    1              1              1
nondemux                   1              1              1
peak_ExonIntron            1             32             32
peak_RMSK                  2             32             32
peak_Transcripts           2             32             32
peak_junctions             1             32             32
peak_process               1             32             32
project_annotations        1              1              1
qc_fastq                   1              1              1
qc_screen_validator        1             32             32
qc_troubleshoot            1              1              1
rename_fastqs              1              1              1
star                       1             32             32
total                     22              1             32
```

- An expected output for test_3 is as follows:
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
Execute pipeline on the cluster
```
#submit to the cluster
sh run_snakemake.sh -p cluster -o /path/to/output/dir/
```

## 5.5 Review outputs
Review the expected outputs on the [Output](https://rbl-nci.github.io/iCLIP/iCLIP/output/) page. If there are errors, review and performing stesp described on the [Troubleshooting](https://rbl-nci.github.io/iCLIP/iCLIP/troubleshooting/) page as needed.