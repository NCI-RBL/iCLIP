import os
import sys
import glob
import os.path
save=os.getcwd()
from subprocess import call

#Before running this script make sure to set the following values:

#Location of directory of files to be worked on. This must be an Absolute path.
rad="46"
transType="SplicedTransc"
input_phred = "/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_"+transType
input_fasta = "/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/genome"
out_dir = "/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_"+transType+"_"+rad+"R"


#Pattern for matching files to work on.  For example: "*.fastq"
match_pattern = "*.fa"

#Name of task being processed.  For example: "fastqc"
task_name = rad+"R_"+transType+"_MakeTransctome_mm10_rRNA"
   

#Find files in specified directory and iterate through the list, creating associated shell scripts
script_name = "mm10_rRNA_scriptloop_MakeTransctome_"+transType+"_"+rad+"R.sh"
completeName = os.path.join(save, script_name)

os.chdir(input_fasta)

script_f = open(completeName, "w+")

script = ""
script += "#swarm -f "+script_name+" --verbose 1 -g 50 -t 32 -b 1 --time 10-00:00:00 --job-name "+task_name+" --module samtools,novocraft"
script += "\n"
script += "\n"
#script += "samtools merge "+out_dir+"/Ro_Clip_iCountcutadpt_all.unique.NH.mm.bam "

for file in glob.glob(match_pattern):
    fname_toks = file.split(".fa")
    script += "java -Xmx100G -jar /data/RBL_NCI/Phil/Tools/USeq_8.9.6/Apps/MakeTranscriptome \\\n"
    script += "-f "+input_fasta+"/ \\\n"
    script += "-u "+input_phred+"/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc."+fname_toks[0]+".phred2 \\\n"
    script += "-r "+rad+" -s \n"
    script += "\n"
    script += "\n"


script += "\n"
script += "\n"
script += "# mv "+input_phred+"/*.gz "+out_dir+" \n"
script += "# gunzip "+out_dir+"/*.gz \n"
script_f.write(script)
script_f.close
 
    #Now that the shell script has been created - submit it...
    #call(["msub", script_name])

