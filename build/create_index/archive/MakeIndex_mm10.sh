

#####################
### 1. Create genome folder
####################

#### comb fasta
cat /data/RBL_NCI/Phil/Reference/fasta/mm10/GRCm38.p6.genome.gencodeM23.fa \
/data/RBL_NCI/Phil/Reference/fasta/mm10/preRibosomal/BK000964.3_TPA_exp.fa \
> /data/RBL_NCI/Phil/Reference/fasta/Combined_FASTA/GRCm38.p6.genome.gencodeM23_prerRNArepeat.fa

## create Fasta file directory
module load ucsc
faSplit byname /data/RBL_NCI/Phil/Reference/fasta/Combined_FASTA/GRCm38.p6.genome.gencodeM23_prerRNArepeat.fa /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/genome/


#############################################################
### 2. Create exon Annotation file
############################################################

#####################
### BK000964.3.bed (rRNA)
bedToGenePred /data/RBL_NCI/Phil/Reference/gtf/mouse/mm10/BK000964.3.bed /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/BK000964.3.bed.phred

### format file correctly
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $2, $3, $4, $5, $6, $7, $8, $9,$10, $5*0 }' /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/BK000964.3.bed.phred > /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/BK000964.3.bed.phred2


#####################
### gencode.vM23.chr_patch_hapl_scaff.annotation.gtf
gtfToGenePred -genePredExt /data/RBL_NCI/Phil/Reference/gtf/mouse/mm10/Gencode_VM23/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred

### format file correctly
awk 'BEGIN {FS="\t"; OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9,$10, $11 }' /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred > /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred2


#####################
### combine gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred2 + BK000964.3.bed.phred2

cat /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.phred2 /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/BK000964.3.bed.phred2 > /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.phred2

#####################
###select only Spliced Transcripts

awk '{if($9>1){print}}' /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.phred2 > /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.phred2


#############################################################
### 3. split by chromosome
#############################################################
# Edit input_dir + out_dir
# script will create a shell script to submit at swarm job on biowulf.
## Splits phred file into individual chromosomes so that MakeTranscriptome step can run. Job will timeout if not split by chromosome

mm10_rRNA_scriptloop_splitPred_SplicedTransc.py
mm10_rRNA_scriptloop_splitPred_FullTransc.py


#####################################################################################################################################
################################################################
### 4. make Exon maksed genome -SplicedTransc
###################################################################

java -Xmx50G -jar /data/RBL_NCI/Phil/Tools/USeq_8.9.6/Apps/MaskExonsInFastaFiles \
-f /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/genome/ \
-u /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.phred2 \
-s /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/genome_maskedexon_SplicedTransc/

awk 'FNR==0{print ""}1' /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/genome_maskedexon_SplicedTransc/*.fasta > /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon.SplicedTransc.fa


#####################################################################################################################################
############################################
### 5. create Transcriptome - SplicedTransc
######################################
python mm10_rRNA_scriptloop_MakeTransctome_SplicedTransc_46R.py
python mm10_rRNA_scriptloop_MakeTransctome_SplicedTransc_71R.py
python mm10_rRNA_scriptloop_MakeTransctome_SplicedTransc_146R.py
## creates bash script for next step. bash script can be submitted as swarm job

mm10_rRNA_scriptloop_MakeTransctome_SplicedTransc_46R.sh
mm10_rRNA_scriptloop_MakeTransctome_SplicedTransc_71R.sh
mm10_rRNA_scriptloop_MakeTransctome_SplicedTransc_146R.sh
## makes transcriptome for each chromsome

## Example for chr1
#java -Xmx100G -jar /data/RBL_NCI/Phil/Tools/USeq_8.9.6/Apps/MakeTranscriptome \
-f /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/genome/ \
-u /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.chr1.phred2 \
-r 71 -s
## Example for read length of 75, -r was addjusted for 50 and 150 (readlength -4)


####################################
## Combine chr transcriptome fasta files
######## 50nt #########
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad46Num100kMin10Splices.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/Chr_split
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad46Num100kMin10Transcripts.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/Chr_split

gunzip /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/Chr_split/*.gz
python Comb_TranscFasta_scriptloop.py
## creates bash script for next step. bash script can be submitted as swarm job

bash Comb_TranscFasta_46R.sh
## combines transcritome files to sinlge fasta for index creation

## Clean up junctions, script removes potential duplicate junctions
### edit pearl script change $file to the combined Combine Transcriptome fasta from Comb_TranscFasta_46R.sh
perl rmv_duplicate_splice_junctions.pl
rm gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad46Num100kMin10Splices.fasta
rm gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad46Num100kMin10Transcripts.fasta

######## 75nt #########
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad71Num100kMin10Splices.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/Chr_split
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad71Num100kMin10Transcripts.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/Chr_split

gunzip /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/Chr_split/*.gz
python Comb_TranscFasta_scriptloop.py
## creates bash script for next step. bash script can be submitted as swarm job

bash Comb_TranscFasta_71R.sh
## combines transcritome files to sinlge fasta for index creation

## Clean up junctions, script removes potential duplicate junctions
### edit pearl script change $file to the combined Combine Transcriptome fasta from Comb_TranscFasta_71R.sh
perl rmv_duplicate_splice_junctions.pl
rm gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad71Num100kMin10Splices.fasta
rm gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad71Num100kMin10Transcripts.fasta


######## 150nt #########
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad146Num100kMin10Splices.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/Chr_split
mv /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/splitphred_SplicedTransc/*Rad146Num100kMin10Transcripts.fasta.gz /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/Chr_split

gunzip /data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/Chr_split/*.gz
python Comb_TranscFasta_scriptloop.py
## creates bash script for next step. bash script can be submitted as swarm job

bash Comb_TranscFasta_146R.sh
## combines transcritome files to sinlge fasta for index creation

## Clean up junctions, script removes potential duplicate junctions
### edit pearl script change $file to the combined Combine Transcriptome fasta from Comb_TranscFasta_146R.sh
perl rmv_duplicate_splice_junctions.pl
rm gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad146Num100kMin10Splices.fasta
rm gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad146Num100kMin10Transcripts.fasta



############################################
### 6. make Index
######################################

#####
## No Splicing
## masked Exon + Transcriptome
novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.nix \
/data/RBL_NCI/Phil/Reference/fasta/Combined_FASTA/GRCm38.p6.genome.gencodeM23_prerRNArepeat.fa



#####
## masked Exon + Transcriptome
novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon_46.nix \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon.SplicedTransc.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad46Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad46Num100kMin10Transcripts.fasta.uniq


novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon_71.nix \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon.SplicedTransc.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad71Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad71Num100kMin10Transcripts.fasta.uniq


novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon_146.nix \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/GRCm38.p6.genome.gencodeM23_BK000964.3.maskedexon.SplicedTransc.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad146Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad146Num100kMin10Transcripts.fasta.uniq

###############################################################################################################################################################################
###############################################################################################################################################################################

#####
## Full genome + Transcriptome
novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.FullGenome_Trnstome.SplicedTransc_46.nix \
/data/RBL_NCI/Phil/Reference/fasta/Combined_FASTA/GRCm38.p6.genome.gencodeM23_prerRNArepeat.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad46Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_46R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad46Num100kMin10Transcripts.fasta.uniq
    

novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.FullGenome_Trnstome.SplicedTransc_71.nix \
/data/RBL_NCI/Phil/Reference/fasta/Combined_FASTA/GRCm38.p6.genome.gencodeM23_prerRNArepeat.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad71Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_71R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad71Num100kMin10Transcripts.fasta.uniq
    
    
novoindex \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/final_index/GRCm38.p6.genome.gencodeM23_BK000964.3.FullGenome_Trnstome.SplicedTransc_146.nix \
/data/RBL_NCI/Phil/Reference/fasta/Combined_FASTA/GRCm38.p6.genome.gencodeM23_prerRNArepeat.fa \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad146Num100kMin10Splices.fasta.uniq \
/data/RBL_NCI/Phil/Reference/Index/novoalign/mm10_rRNA/MakeTranscriptome_genome/MakeTranscriptome_SplicedTransc_146R/gencode.v32.chr_patch_hapl_scaff.annotation_BK000964.3.bed.SplicedTransc.Rad146Num100kMin10Transcripts.fasta.uniq
