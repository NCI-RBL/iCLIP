#!/usr/bin/env bash
set -exo pipefail

#full path to jcounts and countstxt files
if [[ "$#" != "5" ]];then
	echo "bash $0 <jcountsfile> <countstxtfile> <sampleName> <outdir> <pyscriptpath>"
	echo "Please use absolute paths for all inputs!"
	exit
fi
#module load bedtools

keep_files=1
delete_unrequired_intermediate_files=1

# get arguments
jcounts="$1"
countstxt="$2"
sample_name="$3"
outdir="$4"
pyscriptpath="$5"

#jcounts="/data/RBL_NCI/Wolin/11-30-21_7849iCLIP_010822/03_peaks/03_allreadpeaks/KO_allFracMMCounts.txt.jcounts"
#countstxt="/data/RBL_NCI/Wolin/11-30-21_7849iCLIP_010822/03_peaks/03_allreadpeaks/KO_allFracMMCounts.txt"
#outdir="/data/RBL_NCI/Wolin/11-30-21_7849iCLIP_010822/03_peaks/03_allreadpeaks"

jcounts_bn=$(basename $jcounts)
jcounts_filtered="${outdir}/${jcounts_bn}.filtered"

# ensure that site1 and site2 are on the same chr
tail -n+2 ${jcounts} | awk '$3==$6' > ${jcounts_filtered}

# get strand
cut -f1 ${jcounts_filtered} | awk -F"_" '{print $2}' > ${outdir}/${sample_name}_strand
# get site2 coordinates and save to bed file
awk '{OFS="\t";print $6,$7,$7,".","."}' ${jcounts_filtered} > ${outdir}/${sample_name}_site2.bed.tmp
paste ${outdir}/${sample_name}_site2.bed.tmp ${outdir}/${sample_name}_strand | sort -k1,1V -k2,2n | uniq > ${outdir}/${sample_name}_site2.uniq.bed
if [[ "$delete_unrequired_intermediate_files" == "1" ]];then rm -f ${outdir}/site2.bed.tmp;fi

# similarly get site1 bed file
awk '{OFS="\t";print $3,$4,$4,".","."}' ${jcounts_filtered} > ${outdir}/${sample_name}_site1.bed.tmp
paste ${outdir}/${sample_name}_site1.bed.tmp ${outdir}/${sample_name}_strand | sort -k1,1V -k2,2n | uniq > ${outdir}/${sample_name}_site1.uniq.bed
if [[ "$delete_unrequired_intermediate_files" == "1" ]];then rm -f ${outdir}/site1.bed.tmp ${outdir}/strand;fi

# get all peaks
tail -n+3 $countstxt | cut -f1 | awk -F"_" '{print $1}' | sed "s/:/\t/g" | sed "s/-/\t/g" > ${outdir}/${sample_name}_allpeaks.bed.tmp
tail -n+3 $countstxt | cut -f1 > ${outdir}/${sample_name}_name
tail -n+3 $countstxt | cut -f1 | awk -F"_" '{print $2}' > ${outdir}/${sample_name}_strand
paste ${outdir}/${sample_name}_allpeaks.bed.tmp ${outdir}/${sample_name}_name ${outdir}/${sample_name}_strand | awk '{OFS="\t";print $1,$2,$3,$4,".",$5}' | sort -k1,1V -k2,2n | uniq > ${outdir}/${sample_name}_allpeaks.uniq.bed
if [[ "$delete_unrequired_intermediate_files" == "1" ]];then rm -f ${outdir}/allpeaks.bed.tmp ${outdir}/name ${outdir}/strand;fi

# intersect site1 and site2 with peaks to get site-to-peak annotations
cat ${outdir}/${sample_name}_site1.uniq.bed ${outdir}/${sample_name}_site2.uniq.bed | sort -k1,1V -k2,2n | uniq > ${outdir}/${sample_name}_sites.uniq.bed
if [[ "$delete_unrequired_intermediate_files" == "1" ]];then rm -f ${outdir}/site1.uniq.bed ${outdir}/site2.uniq.bed;fi
bedtools intersect -nonamecheck -s -wa -wb -a ${outdir}/${sample_name}_sites.uniq.bed -b ${outdir}/${sample_name}_allpeaks.uniq.bed | awk '{OFS="\t";print $1"_"$2"_"$6,$(NF-2)}' > ${outdir}/${sample_name}_site2peak.lookup.txt

# if bedtools is followed by grep and the bedtools output is zero lines then grep exitcode is non-zero ... causing the pipeline to fail...commenting these lines for now
#bedtools intersect -nonamecheck -loj -s -wa -wb -a ${outdir}/${sample_name}_sites.uniq.bed -b ${outdir}/${sample_name}_allpeaks.uniq.bed | grep "\-1\b" | cut -f1-6 > ${outdir}/${sample_name}_sites_w_no_same_strand_peak.bed
#bedtools intersect -nonamecheck -wa -wb -a ${outdir}/${sample_name}_sites_w_no_same_strand_peak.bed -b ${outdir}/${sample_name}_allpeaks.uniq.bed | awk '{OFS="\t";print $1"_"$2"_"$6,$(NF-2)}' > ${outdir}/${sample_name}_site_to_opposite_strand_peak.lookup.txt
#bedtools intersect -nonamecheck -loj -wa -wb -a ${outdir}/${sample_name}_sites_w_no_same_strand_peak.bed -b ${outdir}/${sample_name}_allpeaks.uniq.bed | grep "\-1\b" | cut -f1-6 > ${outdir}/${sample_name}_sites_w_no_peak.bed

ls -alrth ${outdir}

# switching to python
echo "running python3 script now"
echo "python3 $pyscriptpath ${outdir}/${sample_name}_site2peak.lookup.txt ${jcounts} > ${outdir}/${sample_name}_connected_peaks.txt"
python3 --version
python3 $pyscriptpath ${outdir}/${sample_name}_site2peak.lookup.txt ${jcounts} > ${outdir}/${sample_name}_connected_peaks.txt

# if keep_files == 1 then gzip files and keep them
# else delete
if [[ "$keep_files" == "1" ]];then
    gzip -n -f ${outdir}/${sample_name}_allpeaks.uniq.bed
    gzip -n -f ${outdir}/${sample_name}_sites.uniq.bed
    gzip -n -f ${outdir}/${sample_name}_sites_w_no_same_strand_peak.bed
    gzip -n -f ${outdir}/${sample_name}_site2peak.lookup.txt
    gzip -n -f ${outdir}/${sample_name}_site_to_opposite_strand_peak.lookup.txt
    gzip -n -f ${outdir}/${sample_name}_sites_w_no_peak.bed
else
    rm -f ${outdir}/${sample_name}_allpeaks.uniq.bed
    rm -f ${outdir}/${sample_name}_sites.uniq.bed
    rm -f ${outdir}/${sample_name}_sites_w_no_same_strand_peak.bed
    rm -f ${outdir}/${sample_name}_site2peak.lookup.txt
    rm -f ${outdir}/${sample_name}_site_to_opposite_strand_peak.lookup.txt
    rm -f ${outdir}/${sample_name}_sites_w_no_peak.bed
fi
