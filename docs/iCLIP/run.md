# 3. Running the Pipeline

## 3.1 Pipeline Overview
The Snakemake workflow has a multiple options:
```
Usage: /data/RBL_NCI/Pipelines/iCLIP/[version number]/run_snakemake.sh -p pipeline
	-p options: initialize, dry-run, cluster, local, unlock, git, DAG, report, check
Usage:  -o output_dir
	-o path to output directory
```

## 3.2 Commands explained
The following explains each of the command options:

- Preparation Commands
  - initialize (REQUIRED): This must be performed before any Snakemake run (dry, local, cluster) can be performed. This will copy the necessary config, manifest and Snakefiles needed to run the pipeline to the provided output directory.
  - dry-run (OPTIONAL): This is an optional step, to be performed before any Snakemake run (local, cluster). This will check for errors within the pipeline, and ensure that you have read/write access to the files needed to run the full pipeline.
- Processing Commands
  - local - This will run the pipeline on a local node. NOTE: This should only be performed on an interactive node.
  - cluster - This will submit a master job to the cluster, and subsequent sub-jobs as needed to complete the workflow. An email will be sent when the pipeline begins, if there are any errors, and when it completes.
- Other Commands (All optional)
  - unlock:  This will unlock the pipeline if an error caused it to stop in the middle of a run.
  - git:  This is only utilized for GITHUB Actions testing.
  - DAG: This will produce a DAG of the workflow and dependencies, saved to the /output/dir/log directory
  - report:  This will produce a report generated from the snakemake statistics produced by your pipeline, saved to the /output/dir/log directory.
  - check: This will check for errors in the rules surrounding input manifests. If there are errors they will be printed to the command line and to the output error file.

To run any of these commands, follow the the syntax:
```
sh run_snakemake.sh -p COMMAND -o /path/to/output/dir
```

## 3.3 Typical Workflow
A typical command workflow, running on the cluser, is as follows:

- initialize
- dry-run
- cluster