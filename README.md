# iCLIP

![Build](https://github.com/RBL-NCI/iCLIP/workflows/Tests/badge.svg)  [![GitHub issues](https://img.shields.io/github/issues/RBL-NCI/iCLIP)](https://github.com/RBL-NCI/iCLIP/issues)  [![GitHub license](https://img.shields.io/github/license/RBL-NCI/iCLIP)](https://github.com/RBL-NCI/iCLIP/blob/main/LICENSE)

An RNA Biology pipeline to characterize protein-RNA interactions.

### 1. Getting Started

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (tba later!).

##### 1.1 Download the workflow
Please [clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/RBL-NCI/iCLIP.git
# Change your working directory to the iCLIP repo
cd iCLIP/
```
##### 1.2 Add snakemake to PATH
Please make sure that `snakemake>=5.19` is in your `$PATH`. If you are in Biowulf, please load the following environment module:
```bash
# Recommend running snakemake>=5.19
module load snakemake/5.24.1
```

##### 1.3 Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `snakemake_config.yaml` to configure the workflow execution, `cluster_config.yml` to configure the cluster settings, `multiplex_manifest.tsv` and `samples_manifest.tsv` to specify your sample setup.

##### 1.4 Dry-run the workflow

Please run the following either command to dry-run the snakemake pipeline:
```bash
snakemake -npr -s workflow/Snakefile

sh run_snakemake.sh dry-run
```

### 2. Usage

Submit master job to the cluster:
```bash
sh run_snakemake.sh cluster
```
Submit master job locally:
```bash
sh run_snakemake.sh local
```

### 3. Contribute

This section is for new developers working with the iCLIP pipeline. If you have added new features or adding new changes, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or org account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local filesystem.
3. Copy the modified files to the cloned fork.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) to this repository.
