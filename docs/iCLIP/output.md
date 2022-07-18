# 4. Expected Outputs
The following directories are created under the output_directory:

- 01_preprocess: this directory includes FASTQ and alignment files
    - 01_fastq: FASTQ files, demultiplexed
    - 02_alignment: aligned, indexed, and sorted, alignment files
- 02_bam: this directory includes processed bam files
    - 01_merged: unique and multi-mapped reads BAMS, sorted and indexed
    - 03_dedup: all files in the 02_merged directory, deduplicated
- 03_peaks: this directory includes the bed and SAF files for the pipeline, sorted by:
    - 01_bed: bed files sorted by all reads or unique reads
    - 02_SAF: SAF files sorted by all reads or unique reads
    - 03_counts: peak coutns for all and unique reads, split by unique and MM peaks
- 04_annotation: this directory includes the annotation files at a project and sample level, sorted by:
    - 01_project: includes project level annotation information
    - 02_peaks: includes annotation bed files, complete annotated peak text files
    - final annotation report (HTML) and table (TXT)
- 05_demethod: this directory is only produced when MANORM or DIFFBIND is selected from DE_METHOD
    - 01_input: this includes bed files for any samples being compared
    - 02_analysis: this includes raw DE files (excel MANORM, text DIFFBIND) by comparison
    - 03_report: this includes the final reports (HTML) by comparison
- qc: this directory includes the qc reports, sorted by:
    - multiqc_report: this includes the fastqc results, as well as fastq screen results of each sample before and after filtering
    - qc_report: this includes barcode and alignment information of each sample before and after filtering
- log: this includes log files
    - [date of run]: the slurm output files of the pipeline sorted by pipeline start time; copies of config and manifest files used in this specific pipeline run; error reporting script
    - STAR: star-related log output files