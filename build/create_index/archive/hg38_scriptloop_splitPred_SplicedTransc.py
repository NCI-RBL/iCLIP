import os
import sys
import glob
import os.path
save=os.getcwd()
from subprocess import call

#Before running this script make sure to set the following values:

#Location of directory of files to be worked on. This must be an Absolute path.
input_dir = "/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/genome"
out_dir = "/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc"
 
#Pattern for matching files to work on.  For example: "*.fastq"
match_pattern = "*.fa"

#Name of task being processed.  For example: "fastqc"
task_name = "hg38_scriptloop_splitPred_SplicedTransc"
   

#Find files in specified directory and iterate through the list, creating associated shell scripts
script_name = "hg38_splitPred_SplicedTransc.sh"
completeName = os.path.join(save, script_name)

os.chdir(input_dir)

script_f = open(completeName, "w+")

script = ""
script += "#swarm -f "+script_name+" --verbose 1 -g 50 -t 32 -b 1 --time 24:00:00 --job-name "+script_name+" --module samtools,novocraft"
script += "\n"
script += "\n"
#script += "samtools merge "+out_dir+"/Ro_Clip_iCountcutadpt_all.unique.NH.mm.bam "

for file in glob.glob(match_pattern):
    fname_toks = file.split(".fa")

    script += "grep "+fname_toks[0]+" /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.phred2 > "+out_dir+"/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc."+fname_toks[0]+".phred2 \n"
    
script += "\n"
script_f.write(script)
script_f.close
 
    #Now that the shell script has been created - submit it...
    #call(["msub", script_name])

