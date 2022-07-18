# 2. Preparing Files
The pipeline is controlled through editing configuration and manifest files. Defaults are found in the /output/dir/config and /output/dir/manifest directories, after initialization.

## 2.1 Configs
The configuration files control parameters and software of the pipeline. These files are listed below:

- annotation_config.txt
- cluster_config.yml
- fqscreen_rrna_config.conf
- fqscreen_species_config.conf
- index_config.yaml
- multiqc_config.yaml
- snakemake_config.yaml

### 2.1.1 Annotaiton Config (annotation_config.txt)
The annotation config dictates the type and location of annotaiton files to use for each species (hg38 and mm10). The file must include the following headers:

- type
- AnnoType
- hg38_options
- hg38_selection
- hg38_Custom
- mm10_options
- mm10_selection
- mm10_Custom
- description
- rmsk_flag
- notes
- Status

### 2.1.2 Cluster Config (cluster_config.yml)
The cluster configuration file dictates the resouces to be used during submission to Biowulf HPC. There are two differnt ways to control these parameters - first, to control the default settings, and second, to create or edit individual rules. These parameters should be edited with caution, after significant testing.

### 2.1.3 FASTQ Screen Config (fqscreen_rrna_config.conf and fqscreen_species_config.conf)
The FASTQ sceren configuration files dictates the parameters used for FASTQ Screen.

### 2.1.4 Index Config (index_config.yaml)
The index config is used to dictate what versions of reference index is used for each species (mm10, hg38). The structure of annotation is as follows:
 - organism:
   - std: '/path/to/index/'
   - spliceaware:
     - valuebp1: '/path/to/index1/'
     - valuebp2: '/path/to/index2/'

### 2.1.5 MultiQC Config (multiqc_config.yaml)
The MultiQC screen configuration files control the parameters used for MultiQC.

### 2.1.6 Snakemake Config (snakemake_config.yaml)
There are several groups of parameters that are editable for the user to control the various aspects of the pipeline. These are :

- Folders and Paths
  - These parameters will include the input and ouput files of the pipeline, as well as list all manifest names.
- User parameters
  - These parameters will control the pipeline features. These include thresholds and whether to perform processes.
- STAR parameters
  - These parameters will control the STAR parameters for alignment.
- Modules, container parameters
  - These parameters will control the version of tools used in the pipeline.

## 2.2 Preparing Manifests
There are three manifests, two of which are required for all pipeliens and one that is only required if running a differential expression method. These files describe information on the samples and desired contrasts. The paths of these files are defined in the snakemake_config.yaml file. These files are:

- multiplexManifest
- sampleManifest
- contrastManifest

### 2.2.1 Multiplex Manifest (REQUIRED)
This manifest will include information to map fastq files to their multiple sample ID. It includes the following column headers:

- file_name: the full file name of the multiplexed sample, which must be unique; example: 'test_1.fastq.gz'
- multiplex: the multiplexID associated the fastq file, which must be unique. These names must match the multiplex column of the sampleManifest. example: 'test_1'

An example multplex_manifest.tsv file:
```
file_name,multiplex
test_1.fastq.gz,test_1
test_2.fastq.gz,test_2
```
### 2.2.2 Samples Manifest (REQUIRED)
This manifest will include information to sample level information. It includes the following column headers:

- multiplex: the multiplexID associated with the fasta file, and will not be unique. These names must match the multiplex column of the multiplex_manifest.tsv file. example: 'SIM_iCLIP_S1'
- sample: the final sample name; this column must be unique. example: 'Ro_Clip'
- barcode: the barcode to identify multiplexed sample; this must be unique per each multiplex sample name but can repeat between multiplexid's. example: 'NNNTGGCNN'
- adaptor: the adaptor sequence, to be removed from sample; this may or may not be unique. example: 'AGATCGGAAGAGCGGTTCAG'
- group: groupings for samples, may or may not be unique values. example: 'CNTRL'

An example sampleManifest file with multiplexing of one sample. Notice that the multiplexID test_1 is repeated, as Ro_Clip and Control_Clip are both found in the same fastq file, whereas test_2 is not multiplexed:

```
multiplex,sample,group,barcode,adaptor
test_1,Ro_Clip,CLIP,NNNTGGCNN,AGATCGGAAGAGCGGTTCAG
test_1,Control_Clip,CNTRL,NNNCGGANN,AGATCGGAAGAGCGGTTCAG
test_2,Ro_Clip2,CLIP,NNNCGTANN,AGATCGGAAGAGCGGTTCAG
```

### 2.2.3 Contrast Manifest (REQUIRED with DE_Method of MANORM or DIFFBIND)
This manifest will  include sample or group information to performed differential expresison comparisons (MANORM or DIFFBIND). The column requirements differ by DE method.
- if MANORM:
  - sample: the sample name, identified in the samplesManifest [sample] column, of the sample to compare. example: 'Ro_Clip'
  - background: the background sample name, identified in the samplesManifest [sample] column, of the background to remove. example: 'Control_Clip'
- if DIFFBIND:
  - sample: the sample group, identified in the samplesManifest [group] column, of the sample group to compare. example: 'CLIP' will include samples 'Ro_Clip' and 'Ro_Clip2'
  - background: the background group name, identified in the samplesManifest [group] column, of the background group to remove. example: 'CNTRL' will include sample 'Control_Clip'

An example contrastManifest file for MANORM:
```
sample,background
Ro_Clip,Control_Clip
```

An example contrastManifest file for DIFFBIND:
```
group,background
CLIP,CNTRL
```