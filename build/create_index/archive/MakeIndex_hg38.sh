#####################
### 1. Create genome folder
####################
module load ucsc
faSplit byname /data/RBL_NCI/Phil/Reference/fasta/hg38/Gencode_V32/GRCh38.p13.genome.fa /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/genome/


#############################################################
### 2. Create exon Annotation file
############################################################

gtfToGenePred -genePredExt /data/RBL_NCI/Phil/Reference/gtf/hg38/Gencode_V32/fromGencode/gencode.v32.chr_patch_hapl_scaff.annotation.gtf /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred

awk 'BEGIN {FS="\t"; OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9,$10, $11 }' /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred > /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred2

#####################
###select only Spliced Transcripts

awk '{if($9>1){print}}' /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred2 > /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.phred2


#############################################################
### 3. split by chromosome
#############################################################
# Edit input_dir + out_dir
# script will create a shell script to submit at swarm job on biowulf.
## Splits phred file into individual chromosomes so that MakeTranscriptome step can run. Job will timeout if not split by chromosome

hg38_scriptloop_splitPred_SplicedTransc.py

hg38_splitPred_SplicedTransc.sh
#####################################################################################################################################
################################################################
### 4. make Exon maksed genome -SplicedTransc
###################################################################

java -Xmx50G -jar /data/RBL_NCI/Phil/Tools/USeq_8.9.6/Apps/MaskExonsInFastaFiles \
-f /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/genome/ \
-u /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.phred2 \
-s /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/genome_maskedexon_SplicedTransc/

awk 'FNR==0{print ""}1' /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/genome_maskedexon_SplicedTransc/*.fasta > /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/GRCh38.p13.genome.gencode32.maskedexon.SplicedTransc.fa





#####################################################################################################################################
############################################
### 5. create Transcriptome - SplicedTransc
######################################
hg38_scriptloop_MakeTranscriptome_SplicedTransc_46R.py
hg38_scriptloop_MakeTranscriptome_SplicedTransc_71R.py
hg38_scriptloop_MakeTranscriptome_SplicedTransc_146R.py
## creates bash script for next step. bash script can be submitted as swarm job

hg38_scriptloop_MakeTransctome_SplicedTransc_46R.sh
hg38_scriptloop_MakeTransctome_SplicedTransc_71R.sh
hg38_scriptloop_MakeTransctome_SplicedTransc_146R.sh
## makes transcriptome for each chromsome

##### Example from hg38_scriptloop_MakeTransctome_SplicedTransc_71R.sh
## for chr1

java -Xmx100G -jar /data/RBL_NCI/Phil/Tools/USeq_8.9.6/Apps/MakeTranscriptome \
-f /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/genome/ \
-u /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.chr1.phred2 \
-r 71 -s
## Example for read length of 75, -r was addjusted for 50 and 150 (readlength -4)


####################################
## Combine chr transcriptome fasta files
######## 50nt #########
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad46Num100kMin10Splices.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/Chr_split
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad46Num100kMin10Transcripts.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/Chr_split

gunzip /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/Chr_split/*.gz
python Comb_TranscFasta_scriptloop.py
## creates bash script for next step. bash script can be submitted as swarm job

bash Comb_TranscFasta_46R.sh
## combines transcritome files to sinlge fasta for index creation

## Clean up junctions, script removes potential duplicate junctions
### edit pearl script change $file to the combined Combine Transcriptome fasta from Comb_TranscFasta_46R.sh
perl rmv_duplicate_splice_junctions.pl
rm gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Splices.fasta
rm gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Transcripts.fasta

############################################
######## 75nt #########
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad71Num100kMin10Splices.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/Chr_split
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad71Num100kMin10Transcripts.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/Chr_split

gunzip /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/Chr_split/*.gz
python Comb_TranscFasta_scriptloop.py
## creates bash script for next step. bash script can be submitted as swarm job

bash Comb_TranscFasta_71R.sh
## combines transcritome files to sinlge fasta for index creation

## Clean up junctions, script removes potential duplicate junctions
### edit pearl script change $file to the combined Combine Transcriptome fasta from Comb_TranscFasta_46R.sh
perl rmv_duplicate_splice_junctions.pl
rm gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad71Num100kMin10Splices.fasta
rm gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad71Num100kMin10Transcripts.fasta

############################################
######## 125nt #########
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad146Num100kMin10Splices.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/Chr_split
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad146Num100kMin10Transcripts.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/Chr_split

gunzip /data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/Chr_split/*.gz
python Comb_TranscFasta_scriptloop.py
## creates bash script for next step. bash script can be submitted as swarm job

bash Comb_TranscFasta_146R.sh
## combines transcritome files to sinlge fasta for index creation

## Clean up junctions, script removes potential duplicate junctions
### edit pearl script change $file to the combined Combine Transcriptome fasta from Comb_TranscFasta_46R.sh
perl rmv_duplicate_splice_junctions.pl
rm gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad146Num100kMin10Splices.fasta
rm gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad146Num100kMin10Transcripts.fasta


############################################
### 6. make Index
######################################

#####
## No Splicing
## masked Exon + Transcriptome
novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/GRCh38.p13.genome.gencode32_BK000964.3.nix \
/data/RBL_NCI/Phil/Reference/fasta/hg38/Gencode_V32/GRCh38.p13.genome.fa


#####
## masked Exon + Transcriptome
novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.maskedexon_46.nix \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/GRCh38.p13.genome.gencode32.maskedexon.SplicedTransc.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Transcripts.fasta.uniq


novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.maskedexon_71.nix \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/GRCh38.p13.genome.gencode32.maskedexon.SplicedTransc.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad71Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad71Num100kMin10Transcripts.fasta.uniq


novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.maskedexon_146.nix \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/GRCh38.p13.genome.gencode32.maskedexon.SplicedTransc.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad146Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad146Num100kMin10Transcripts.fasta.uniq


#####
## Full genome + Transcriptome

novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc_46.nix \
/data/RBL_NCI/Phil/Reference/fasta/hg38/Gencode_V32/GRCh38.p13.genome.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad46Num100kMin10Transcripts.fasta.uniq


novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc_71.nix \
/data/RBL_NCI/Phil/Reference/fasta/hg38/Gencode_V32/GRCh38.p13.genome.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad71Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad71Num100kMin10Transcripts.fasta.uniq


novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/final_index/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc_146.nix \
/data/RBL_NCI/Phil/Reference/fasta/hg38/Gencode_V32/GRCh38.p13.genome.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad146Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/hg38/hg38_gencode_v32/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.SplicedTransc.Rad146Num100kMin10Transcripts.fasta.uniq
