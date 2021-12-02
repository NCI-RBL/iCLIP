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

- test1: ICLIP
  - iclip default settings (see above)
- test2: R FLAG
  - R flag set to 5
- test3: R FLAG
  - R flag set to 20
- test4: STRICTER VALUES
  - t set to 20,3;  g set to 30; R set to 5
- test5: NOVO DEFAULT
  - t set to 20,3; l removed; x set to 6; g set to 40; R set to 5; S set to 2
- test6: S FLAG
  - s flag set to 0
- test7: S FLAG
  - s flag set to 2
- test8: S FLAG
  - no s flag
- test9: K FLAG
  - no k flag

## Testing details
- input files are taken from the sample 
- the first 10 of the splits (out of 98) were used for the analysis
  - taken from ls /data/RBL_NCI/Wolin/10212021_mESC_clip_2/01_preprocess/01_splits/
- index file was the splice aware, bp = 75
  - /data/CCBR_Pipeliner/iCLIP/index/active/phil/hg38/ gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc_75.nix

## Directories
input_dir="/data/RBL_NCI/Wolin/Sam/novoalign/splits"
output_dir="/data/RBL_NCI/Wolin/Sam/novoalign/aligned"
log_dir="/data/RBL_NCI/Wolin/Sam/novoalign/log"

####################################################################################
# Workflow
####################################################################################
############## Commands
align
cleanup
unique_mm
merge_splits
merge_um
dedup
beds
counts
project_annotations (sh only, will submit one per project)
peak_annotations
annotation_report

############## FILE CHECK
### check previous files exist
sh /home/sevillas2/git/iCLIP/build/novoalign_testing/alignment_check.sh v1 

############## SH FILES
### create individual sh files, one file for each test (N=9) and each split (N=10); total of 90 files
sh /home/sevillas2/git/iCLIP/build/novoalign_testing/alignment_params_sh.sh v1 

############## SBATCH FILES
### create sbatch file which includes each test's sh files; total of 9 files
sh /home/sevillas2/git/iCLIP/build/novoalign_testing/alignment_params_sbatch.sh v1 

############## SUBMIT TO CLUSTER
### submit sbatch jobs to cluster
sh /home/sevillas2/git/iCLIP/build/novoalign_testing/alignment_params_submit.sh v1 
