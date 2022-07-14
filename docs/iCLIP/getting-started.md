# Getting Started
The iCLIP github repository is stored locally, and will be used for project deployment. Multiple projects can be deployed from this one point simultaneously, without concern.

```
#Change working directory to the iCLIP repository

cd /data/RBL_NCI/Pipelines/iCLIP/[version number]
```

1. Getting Started
1.1 Introduction
The iCLIP Pipelie beings with raw FASTQ files and performs demultiplexing, adaptor removal, and trimming, as designated by the user. Alignment is then performed using [STAR](), followed by deduplication. Peaks are then called, where both unique and multimapped read counts are identified. Peaks are then collapsed, annotated, and summarized into reports. If designated either MANORM or DIFFBIND is performed on peak findings. QC reports are also generated with each project.

The following are sub-commands used within iCLIP:

- initialize: initalize the pipeline
- dryrun: predict the binding of peptides to any MHC molecule
- cluster: execute the pipeline on the Biowulf HPC
- local: execute a local, interactive, session
- git: execute GitHub actions
- unlock: unlock directory
- DAG: create DAG report
- report: create SNAKEMAKE report

1.2 Setup Dependencies¶
iCLIP has several dependencies listed below. These dependencies can be installed by a sysadmin. All dependencies will be automatically loaded if running from Biowulf.

- bedtools: "bedtools/2.29.2"
- bowtie2: "bowtie/2-2.3.4"
- fastq_screen: "fastq_screen/0.14.0"
- fastqc: "fastqc/0.11.9"
- manorm: "manorm/1.1.4"
- multiqc: "multiqc/1.9"
- perl: "perl/5.24.3"
- python: "python/3.8"
- R: "R/4.0"
- samtools: "samtools/1.11"
- star: "STAR/2.7.8a"
- subread: "subread/2.0.1"
- ultraplex: "ultraplex/1.2.5"
- umitools: "umitools/1.1.1"

1.3 Login to the cluster¶
iCLIP has been exclusively tested on Biowulf HPC. Login to the cluster's head node and move into the pipeline location.
```
# ssh into cluster's head node
ssh -Y $USER@biowulf.nih.gov

# move the iCLIP dir
cd /
```

1.4 Load an interactive session 
An interactive session should be started before performing any of the pipeline sub-commands, even if the pipeline is to be executed on the cluster.

```
# Grab an interactive node
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
```