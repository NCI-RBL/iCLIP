# SUPPLEMENTARY DATA 1
# Busch et al, Methods, 2019

# This file provides the bash code for basic read processing, conversion into crosslink events and peak calling,
# as described in Chapters 3 - 5.1.

# Details on the input files, preset variables and external tools required to run this code are listed in Chapters 2 and 3.1.

output_dir="/data/sevillas2/test_2/geos"
tmp=$output_dir/tmp

sample_id="SRR5646572"
raw_fq="$output_dir/fq/${sample_id}.fastq.gz"

#block1
filtered_IDS="$output_dir/block_1/1_${sample_id}.filteredIDS.list"
filtered_fq="$output_dir/block_1/2_${sample_id}.filtered.fastq.gz"
barcodes_found="$output_dir/block_1/3_${sample_id}.expbarcodes.list"
barcodes_fasta="$output_dir/block_2/4_${sample_id}.barcodes.fasta"

demux_fq="$output_dir/block_2/5_SRR5646572.fastq.gz"
demux_fa="$output_dir/block_2/6_SRR5646572.fasta.gz"
read_png="$output_dir/block_2/7_SRR5646572.readlength.png"

barcodeLength=9
minBaseQuality=15
umi_len_1=3
bc_len=4
umi_len_2=2
adapter_seq="AGATCGGAAGAGCGGTTCAG"
min_read=15
flexbar_out="$output_dir/block_1/SRR5646572"

dir_list=($tmp)
for pd in "${dir_list[@]}"; do if [[ ! -d $pd ]]; then mkdir -p $pd; fi; done

block_1="Y"
block_2="N"

##====================
## Quality control
##====================

###--------------------------
### General quality check
###--------------------------
if [[ $block_1 == "Y" ]]; then
    module load fastqc fastxtoolkit seqtk
    
    #### FastQC on the full data set
    fastqc --extract --nogroup --outdir $output_dir $raw_fq

    ###------------------------------------------
    ### Quality filter on the barcode regions
    ###------------------------------------------

    #### List of read IDs of reads with high quality barcode regions (using FASTX-Toolkit)
    zcat $raw_fq | fastx_trimmer -l $barcodeLength | fastq_quality_filter -q $minBaseQuality -p 100 | awk 'FNR%4==1 { print $1 }' | sed 's/@//' > $filtered_IDS

    #### Extract reads of given read IDs (using seqtk) and remove problematic characters and whitespaces from read IDs#S
    seqtk subseq $raw_fq $filtered_IDS | sed 's/ /#/g; s/\\//#/g' | gzip > $filtered_fq

    ###------------------------
    ### Barcode frequencies 
    ###------------------------

    #### Extract all detected experimental barcodes and their frequencies (x = length of UMI1, y = length of the experimental barcodes)
    zcat $filtered_fq | awk -v umi1_len=$umi_len_1 -v exp_bc_len=$bc_len '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > $barcodes_found
fi

if [[ $block_2 == "Y" ]]; then
    module load flexbar fastxtoolkit
    ##=================================================
    ## Demultiplexing, adapter and barcode trimming
    ##=================================================

    #### Demultiplexing, adapter and barcode trimming using Flexbar
    flexbar -r $filtered_fq --zip-output GZ --barcodes $barcodes_fasta --barcode-unassigned --barcode-trim-end LTAIL \
    --barcode-error-rate 0 --adapter-seq $adapter_seq --adapter-trim-end RIGHT --adapter-error-rate 0.1 --adapter-min-overlap 1 \
    --min-read-length $min_read --umi-tags -t $flexbar_out


    #### Plot reads length distribution using FASTX-Toolkit
    ##### fastq to fasta
    zcat $demux_fq | fastq_to_fasta -n -r | gzip > $demux_fa

    ##### create the plot
    fasta_clipping_histogram.pl $demux_fa $read_png
fi

if [[ $block_3 == "Y"]]; then
    ##====================
    ## Genomic mapping
    ##====================

    #### Genomic mapping using STAR

    #STAR --runMode alignReads --genomeDir genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile <annotation.gtf> --sjdbOverhang maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn <sampleX.fastq.gz>
    sh 
fi


# ##======================================
# ## Duplicate removal (deduplication)
# ##======================================

# #### Duplicate removal (deduplication) using UMI-tools

# umi_tools dedup -I <sampleX.bam> -L <sampleX.duprm.log> -S <sampleX.duprm.bam> --extract-umi-method read_id --method unique



# ##==========================================
# ## Extraction of crosslinked nucleotides
# ##==========================================

# #### Create file of chromosome length using SAMtools (produces a file genome.fasta.fai)

# samtools faidx <genome.fasta>


# #### Convert all read locations to intervals in bed file format using BEDTools
 
# bedtools bamtobed -i <sampleX.duprm.bam> > <sampleX.bed>
 

# #### Shift intervals depending on the strand by 1 bp upstream using BEDTools
 
# bedtools shift -m 1 -p -1 -i <sampleX.bed> -g <genome.fasta.fai> > <sampleX.shifted.bed>
 

# #### Extract the 5' end of the shifted intervals and pile up into coverage track in bedgraph file format (separately for each strand) using BEDTools (in case of RPM-normalised coverage tracks, use additional parameter -scale with 1,000,000/#mappedReads)
 
# bedtools genomecov -bg -strand + -5 -i <sampleX.shifted.bed> -g <genome.fasta.fai> > <sampleX.plus.bedgraph>
# bedtools genomecov -bg -strand - -5 -i <sampleX.shifted.bed> -g <genome.fasta.fai> > <sampleX.minus.bedgraph>


# #### Optional convertion of bedgraph files to bw file format files using bedGraphToBigWig of the kentUtils suite
 
# bedGraphToBigWig <sampleX.strand.bedgraph> <genome.fasta.fai> <sampleX.strand.bw>
 

# #### Depending on the system and the version of bedGraphToBigWig, it might be necessary to sort the bedgraph files before converting them to bw files:
 
# export LC_COLLATE=C
# sort -k1,1 -k2,2n <sampleX.strand.bedgraph> > <sampleX.strand.sorted.bedgraph>



# ##========================================================
# ## Diagnostic plots and measures of library complexity
# ##========================================================

# ###------------------------------
# ### Duplicate removal summary
# ###------------------------------

# #### Number of crosslink events, i.e. reads after duplicate removal:
 
# cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalcount=0 }{ totalcount += (($3 - $2) * $4) }END{ print totalcount }'
 

# #### Number of crosslinked nucleotides, i.e. positions harbouring crosslinked nucleotides (if both strands are covered, count as 2):
 
# cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalpos=0 }{ totalpos += ($3 - $2) }END{ print totalpos }'


# #### Number of stacked crosslink events, i.e. crosslink events on positions with >1 crosslink events:
 
# cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalstackedcount=0 }{ if($4 > 1) totalstackedcount += (( $3 - $2) * $4) }END{ print totalstackedcount }'


# #### Number of nucleotides with stacked crosslink events, i.e. positions with >1 crosslink events: 
 
# cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalstackedpos=0 }{ if($4 > 1) totalstackedpos += ($3 - $2) }END{ print totalstackedpos }' 


# ###----------------------------------------
# ### Reads with insertions and deletions
# ###----------------------------------------

# #### Convert bam to sam file format (using SAMtools)

# samtools view <sampleX.duprm.bam> -o <sampleX.duprm.sam>


# #### Number of reads mapped with deletions:

# cut -f6 <sampleX.duprm.sam> | grep D | wc -l


# #### Number of reads mapped with insertions:
 
# cut -f6 <sampleX.duprm.sam> | grep I | wc -l


# ###---------------------
# ### iCLIPro analysis
# ###---------------------

# #### iCLIPro analysis and plots

# iCLIPro -r 50 -b 300 -f 50 -g "L15:15,L16:16,L17:17,L18:18,L19:19,L20:20,L21:21,L22:22,L23:23,L24:24,L25:25,L26:26,L27:27,L28:28,L29:29,L30:30,L31:31,L32:32,L33:33,L34:34,L35:35,L36:36,L37:37,L38:38,L39:39,L40:40,R:41" -p "L15-R,L16-R,L17-R,L18-R,L19-R,L20-R,L21-R,L22-R,L23-R,L24-R,L25-R,L26-R,L27-R,L28-R,L29-R,L30-R,L31-R,L32-R,L33-R,L34-R,L35-R,L36-R,L37-R,L38-R,L39-R,L40-R" -o $output_dir <sampleX.duprm.bam>



# ##===============================
# ## Peak calling with PureCLIP
# ##===============================

# #### Merge BAM files (sampleX.duprm.bam)

# samtools merge -f <merged.bam> -b <list_of_bam_files>

# samtools index <merged.bam>


# #### Run PureCLIP

# pureclip -i <merged.bam> -bai <merged.bam.bai> -g <genome.fasta> -ld -nt 8 -o PureCLIP.crosslink_sites.bed -or PureCLIP.crosslink_regions.bed


# #### Remove 7th column of PureCLIP output file

# cat PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6 > PureCLIP.crosslink_sites_short.bed
