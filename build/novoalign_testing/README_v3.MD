# Novoalign Testing

http://www.novocraft.com/documentation/novoalign-2/novoalign-user-guide/novoalign-command-options/

####################################################################################
# OVERVIEW
####################################################################################
## variable summary
-- k Enables quality calibration. Default is off.
-- t Sets the maximum alignment score acceptable for the best alignment. Default 'automatic'
-- l Sets the minimum number of good quality bases for a read. Alignment will not be attempted for reads with less bases. Default log4(Ng) + 5 where Ng is the length of the reference genome. Default is log4(3E9) + 5 = ~21.
-- x Sets the gap extend penalty. Default 6
-- g Sets the gap opening penalty. Default 40
-- s Turns on read trimming and sets step size. Reads that fail to align will be progressively shortened by specified amount until they either align or length reduces to less that the length set by the -l. Default is 2
-- R Sets score difference between best and second best alignment for calling a repeat. Default 5.

## iCLIP default settings
- k
- t 15,3
- l 20
- x 4 
- g 20
- s 1
- R 0

####################################################################################
# TESTING INFO
####################################################################################
## Testing summary
## Tests have been parsed down from v1 and v2 parameters

- test1: ICLIP
  - iclip default settings (see above)
- test4: STRICTER VALUES
  - t set to 20,3;  g set to 30; R set to 5
- test5: NOVO DEFAULT
  - t set to 20,3; l removed; x set to 6; g set to 40; R set to 5; S set to 2

## Testing details
- input files are taken from two studies
    - Dir 1: /data/RBL_NCI/Wolin/mESC_clip_2/fastq
      - demux from 
          - 1_PARiCLIP_mRO1_S1_R1_001.fastq.gz	clip1
      - sample 
          - clip1	Ro_Clip_1	CLIP	NNNNNGCCTAANNNN	AGATCGGAAGAGCGGTTCAG
      - index file was the splice aware, bp = 75, mm10
        - /data/CCBR_Pipeliner/iCLIP/index/active/phil/mm10/mm10_splice75bp_unmasked.nix
    - Dir 2: 
      - 1 sample from  /data/RBL_NCI/Wolin/8-0reprocess/_fCLIP/01_preprocess
        - WT_fCLIP_filtered.fastq.gz  WT_fCLIP
      - index file was the splice aware, bp = 75, hg38
        - /data/CCBR_Pipeliner/iCLIP/index/active/phil/hg38/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc_75.nix

## Directories
input_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/splits"
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/aligned"
log_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/log"

####################################################################################
# Workflow
####################################################################################
Run complete iCLIP pipeline with adjustments for alignment params only

sh run_snakemake.sh -p cluster -o /data/RBL_NCI/Wolin/Sam/novoalign_v3/hg38/test1/ 27643847
sh run_snakemake.sh -p cluster -o /data/RBL_NCI/Wolin/Sam/novoalign_v3/hg38/test4/
sh run_snakemake.sh -p cluster -o /data/RBL_NCI/Wolin/Sam/novoalign_v3/hg38/test5/

sh run_snakemake.sh -p cluster -o /data/RBL_NCI/Wolin/Sam/novoalign_v3/mm10/test1/ 27643853
sh run_snakemake.sh -p cluster -o /data/RBL_NCI/Wolin/Sam/novoalign_v3/mm10/test4/
sh run_snakemake.sh -p cluster -o /data/RBL_NCI/Wolin/Sam/novoalign_v3/mm10/test5/