#!/bin/bash
parent_dir="/data/sevillas2/test_2/geos"
fq_dir="$parent_dir/fq"
pureclip_dir="/data/CCBR_Pipeliner/bin/PureCLIP/PureCLIP/build/pureclip"
fa_file="/data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/hg38/ref/GRCh38.p13.genome.fa"

#samples
sample_list=("SRR5646572" "SRR5646573" "SRR5646574" "SRR5646575")

# clean ref
clean_reference="/data/RBL_NCI/Wolin/Sam/pureclip/GRCh38.p13.genome.fa"
if [[ ! -f $clean_reference ]]; then
    echo "** Creating REF **"
    cat $reference_file | sed '/^\>/! s/[RYKMSWBVHD]/N/g' > $clean_reference
fi


##############################
# control
##############################
# controls
#run_trim="Y"
#run_demux="Y"
#run_star="Y"
run_dedup="Y"
run_crosslink="N"
run_summary="N"
run_pro="N"
run_pure="N"

##############################
# setup
##############################
##########
if [[ $run_ == "Y" ]]; then
    echo "** Running  **"

    sub_dir=$parent_dir/trim
    if [[ ! -f $sub_dir/tmp ]]; then mkdir -p $sub_dir/tmp; fi
fi

##############################
# function
##############################
function star_run()
{
    local star_command="#!/bin/bash
    
            module load STAR
            
            STAR \
            --runMode alignReads \
            --genomeDir /data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/hg38/index \
            --sjdbGTFfile /data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/hg38/ref/gencode.v32.annotation.gtf \
            --readFilesCommand zcat \
            --readFilesIn $input_file \
            --outFileNamePrefix $sub_dir/tmp/${sample_id}_ \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --alignEndsType Local \
            --alignIntronMax 50000 \
            --alignSJDBoverhangMin 3 \
            --alignSJoverhangMin 5 \
            --alignTranscriptsPerReadNmax=10000 \
            --alignWindowsPerReadNmax=10000 \
            --outFilterMatchNmin 15 \
            --outFilterMatchNminOverLread 0.9 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterMultimapNmax 10000 \
            --outFilterMultimapScoreRange 0 \
            --outFilterScoreMin 0 \
            --outFilterType Normal \
            --outSAMattributes All \
            --outSAMunmapped None \
            --outSJfilterCountTotalMin 3 1 1 1 \
            --outSJfilterOverhangMin 30 12 12 12 \
            --outSJfilterReads All \
            --seedMultimapNmax=10000 \
            --seedNoneLociPerWindow=20 \
            --seedPerReadNmax=10000 \
            --seedPerWindowNmax=10000 \
            --sjdbScore 2 \
            --winAnchorMultimapNmax=10000"
    echo "$star_command"
}
function umi_command()
{
    local umi_cmd="#!/bin/bash
            module load umitools

            umi_tools dedup \
            -I $aligned_file \
            -L $log_file \
            -S $deduped_file \
            --extract-umi-method read_id \
            --method unique"
    echo "$umi_cmd"
}
function pureclip_command()
{
    local pc_cmd="$pureclip_dir \
        -i $pro_file \
        -bai ${pro_file}.bai \
        -g $clean_reference \
        -ld -nt 8 \
        -o ${sub_dir}/${sample_id}_sites.bed \
        -or ${sub_dir}/${sample_id}_regions.bed"
    echo "$pc_cmd"
}

##############################
# run
##############################
# trimming
if [[ $run_trim == "Y" ]]; then
    echo "** CHECKING TRIMMING **"
    #https://github.com/agordon/fastx_toolkit/blob/master/scripts/fasta_clipping_histogram.pl
    module load fastxtoolkit seqtk fastqc

    sub_dir=$parent_dir/01_trim
    if [[ ! -f $sub_dir/tmp ]]; then mkdir -p $sub_dir/tmp; fi

    for sample_id in "${sample_list[@]}"; do
        input_file=$fq_dir/$sample_id.fastq.gz
        tmp_file=$sub_dir/tmp/$sample_id.qualFilteredIDS.list
        filtered_file=$sub_dir/$sample_id.filtered.fastq.gz
        fastqc_file=$sub_dir/$sample_id.filtered_fastqc.html
        barcode_file=$sub_dir/$sample_id.exp_barcodes

        if [[ ! -f $fastqc_file ]]; then
            echo "--$sample_id"
            #### List of read IDs of reads with high quality barcode regions (using FASTX-Toolkit)
            #zcat <data.fastq.gz> | fastx_trimmer -l barcodeLength | fastq_quality_filter -q minBaseQuality -p 100 | awk 'FNR%4==1 { print $1 }' | sed 's/@//' > <tmp/data.qualFilteredIDs.list>
            zcat $input_file | fastx_trimmer -l 9 | fastq_quality_filter -q 10 -p 100 | awk 'FNR%4==1 { print $1 }' | sed 's/@//' > $tmp_file

            #### Extract reads of given read IDs (using seqtk) and remove problematic characters and whitespaces from read IDs
            #seqtk subseq <data.fastz.gz> <tmp/data.qualFilteredIDs.list> | sed 's/ /#/g; s/\\//#/g' | gzip > <data.filtered.fastq.gz>
            seqtk subseq $input_file $tmp_file | sed 's/ /#/g; s/\\//#/g' | gzip > $filtered_file

            #### run FASTQ to ensure adaptors were already removed in experiment
            fastqc $filtered_file

            #unncessary step since barcodes and adaptors were already removed
            #### Extract all detected experimental barcodes and their frequencies (x = length of UMI1, y = length of the experimental barcodes)
            #zcat <data.filtered.fastq.gz> | awk -v umi1_len=x -v exp_bc_len=y '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > <exp_barcodes.detected>
            #zcat $filtered_file | awk -v umi1_len=3 -v exp_bc_len=4 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > $barcode_file
        fi
    done

    # remove tmps
    rm -r $sub_dir/tmp
fi

# demultiplex
if [[ $run_demux == "Y" ]]; then
    echo "** Checking DEMUX **"
    module load fastxtoolkit flexbar
    
    sub_dir=$parent_dir/02_demux
    if [[ ! -f $sub_dir ]]; then mkdir -p $sub_dir; fi

    for sample_id in "${sample_list[@]}"; do
        filtered_file=$parent_dir/01_trim/$sample_id.filtered.fastq.gz
        barcode_file=$parent_dir/01_trim/$sample_id.barcodes.fasta
        demux_fastq=$sub_dir/$sample_id.fastq.gz
        demux_fasta=$sub_dir/$sample_id.fasta.gz
        png_file=$sub_dir/$sample_id.readlength.png

        if [[ ! -f $png_file ]]; then
            echo "--$sample_id"
            #### Demultiplexing, adapter and barcode trimming using Flexbar
            #https://github-wiki-see.page/m/seqan/flexbar/wiki/Manual
            #flexbar -r <data.filtered.fastq.gz> --zip-output GZ --barcodes barcodes.fasta --barcode-unassigned --barcode-trim-end LTAIL --barcode-error-rate 0 --adapter-seq adapter.seq --adapter-trim-end RIGHT --adapter-error-rate 0.1 --adapter-min-overlap 1 --min-read-length minReadLength --umi-tags
            cmd=`echo "flexbar -r $filtered_file \
            --zip-output GZ \
            --barcodes $barcode_file \
            --barcode-unassigned \
            --barcode-trim-end LTAIL \
            --barcode-error-rate 0 \
            --adapter-seq AGATCGGAAGAGCGGTTCAG \
            --adapter-trim-end RIGHT \
            --adapter-error-rate 0.1 \
            --adapter-min-overlap 1 \
            --min-read-length 15 \
            --target ${sub_dir}/$sample_id \
            --umi-tags \
            --length-dist \
            -t ${sub_dir}/${sample_id}"`
            echo $cmd

            #### Plot reads length distribution using FASTX-Toolkit
            ##### fastq to fasta
            #zcat <sampleX.fastq.gz> | fastq_to_fasta -n -r | gzip > <sampleX.fasta.gz>
            zcat $demux_fastq | fastq_to_fasta -n -r | gzip > $demux_fasta

            ##### create the plot
            fasta_clipping_histogram.pl $demux_fasta $png_file
        fi
    done
fi

# STAR
if [[ $run_star == "Y" ]]; then
    echo "** Checking STAR **"

    sub_dir=$parent_dir/03_STAR
    if [[ ! -f $sub_dir/log ]]; then mkdir -p $sub_dir/log; fi
    if [[ ! -f $sub_dir/tmp ]]; then mkdir -p $sub_dir/tmp; fi

    for sample_id in "${sample_list[@]}"; do
        input_file=$parent_dir/01_trim/$sample_id.filtered.fastq.gz
        aligned_file=$sub_dir/$sample_id.aligned.bam
        log_file=$sub_dir/log/$sample_id.log
        sh_file=$sub_dir/log/$sample_id.sh

        if [[ ! -f $aligned_file ]]; then
            echo "--$sample_id"
            ### Genomic mapping using STAR
            #STAR --runMode alignReads --genomeDir genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 
            #--sjdbGTFfile <annotation.gtf> --sjdbOverhang maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate 
            #--readFilesIn <sampleX.fastq.gz>
            cmd=`star_run`
            # echo $cmd > $sh_file
            # echo "mv $sub_dir/tmp/${sample_id}_Aligned.sortedByCoord.out.bam $aligned_file" >> $sh_file
            # echo "mv $sub_dir/tmp/${sample_id}_Log.final.out $log_file" >> $sh_file
            
            sh /home/sevillas2/general_code/submit.sbatch.swarm sbatch $sub_dir $sh_file
        fi
    done
fi

#DEDUP
if [[ $run_dedup == "Y" ]]; then
    echo "** Running  **"

    sub_dir=$parent_dir/04_dedup
    if [[ ! -f $sub_dir/log ]]; then mkdir -p $sub_dir/log; fi

    for sample_id in "${sample_list[@]}"; do
        aligned_file=$parent_dir/STAR/$sample_id.aligned.bam
        log_file=$sub_dir/log/$sample_id.log
        deduped_file=$sub_dir/$sample_id.deduped.bam
        sh_file=$sub_dir/log/$sample_id.sh

        if [[ ! -f $aligned_file ]]; then
            echo "--$sample_id"
            cmd=`umi_command`
            #echo $cmd > $sh_file
            
            #sh /home/sevillas2/general_code/submit.sbatch.swarm sbatch $sub_dir $sh_file
        fi
    done
fi

# CROSSLINKED EVENTS
if [[ $run_crosslink == "Y" ]]; then
    sub_dir=$parent_dir/05_crosslink
    if [[ ! -f $sub_dir ]]; then mkdir -p $sub_dir; fi

    for sample_id in "${sample_list[@]}"; do
        deduped_file=$parent_dir/04_dedup/$sample_id.deduped.bam
        bed_file=$sub_dir/$sample_id.bed
        shiftedbed_file=$sub_dir/${sample_id}_shifted.bed
        bedgraph_p_file=$sub_dir/${sample_id}_p.bedgraph
        bedgraph_n_file=$sub_dir/${sample_id}_n.bedgraph

        if [[ ! -f $bedgraph_n_file ]]; then
            echo "--$sample_id"
            module load samtools bedtools

            #### Create file of chromosome length using SAMtools (produces a file genome.fasta.fai)
            # samtools faidx <genome.fasta>
            samtools faidx $fa_file

            #### Convert all read locations to intervals in bed file format using BEDTools
            # bedtools bamtobed -i <sampleX.duprm.bam> > <sampleX.bed>
            bedtools bamtobed -i $deduped_file > $bed_file


            #### Shift intervals depending on the strand by 1 bp upstream using BEDTools
            # bedtools shift -m 1 -p -1 -i <sampleX.bed> -g <genome.fasta.fai> > <sampleX.shifted.bed>
            bedtools shift -m 1 -p -1 -i $bed_file -g $fa_file.fai > $shiftedbed_file
                

            #### Extract the 5' end of the shifted intervals and pile up into coverage track in bedgraph file format (separately for each strand) using BEDTools (in case of RPM-normalised coverage tracks, use additional parameter -scale with 1,000,000/#mappedReads)
            # bedtools genomecov -bg -strand + -5 -i <sampleX.shifted.bed> -g <genome.fasta.fai> > <sampleX.plus.bedgraph>
            # bedtools genomecov -bg -strand - -5 -i <sampleX.shifted.bed> -g <genome.fasta.fai> > <sampleX.minus.bedgraph>
            bedtools genomecov -bg -strand + -5 -i $shiftedbed_file -g $fa_file.fai > $bedgraph_p_file
            bedtools genomecov -bg -strand - -5 -i $shiftedbed_file -g $fa_file.fai > $bedgraph_n_file

            #### Optional convertion of bedgraph files to bw file format files using bedGraphToBigWig of the kentUtils suite
            # bedGraphToBigWig <sampleX.strand.bedgraph> <genome.fasta.fai> <sampleX.strand.bw>                

            #### Depending on the system and the version of bedGraphToBigWig, it might be necessary to sort the bedgraph files before converting them to bw files:
            # export LC_COLLATE=C
            # sort -k1,1 -k2,2n <sampleX.strand.bedgraph> > <sampleX.strand.sorted.bedgraph>

        fi
    done
fi

# run stats summary
if [[ $run_summary == "Y" ]]; then
    sub_dir=$parent_dir/06_summary
    if [[ ! -f $sub_dir ]]; then mkdir -p $sub_dir; fi

    for sample_id in "${sample_list[@]}"; do
        bedgraph_p_file=$parent_dir/05_crosslink/${sample_id}_p.bedgraph
        bedgraph_n_file=$parent_dir/05_crosslink/${sample_id}_n.bedgraph
        deduped_file=$parent_dir/04_dedup/$sample_id.deduped.bam
        deduped_sam_file=$sub_dir/$sample_id.deduped.sam

        if [[ ! -f $deduped_sam_file ]]; then
            echo "--$sample_id"
            module load samtools

            #### Number of crosslink events, i.e. reads after duplicate removal:
            # cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalcount=0 }{ totalcount += (($3 - $2) * $4) }END{ print totalcount }'
            cat $bedgraph_p_file $bedgraph_n_file | awk 'BEGIN{ totalcount=0 }{ totalcount += (($3 - $2) * $4) }END{ print totalcount }'
            
            #### Number of crosslinked nucleotides, i.e. positions harbouring crosslinked nucleotides (if both strands are covered, count as 2):
            # cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalpos=0 }{ totalpos += ($3 - $2) }END{ print totalpos }'
            cat $bedgraph_p_file $bedgraph_n_file  | awk 'BEGIN{ totalpos=0 }{ totalpos += ($3 - $2) }END{ print totalpos }'

            #### Number of stacked crosslink events, i.e. crosslink events on positions with >1 crosslink events:
            # cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalstackedcount=0 }{ if($4 > 1) totalstackedcount += (( $3 - $2) * $4) }END{ print totalstackedcount }'
            cat $bedgraph_p_file $bedgraph_n_file | awk 'BEGIN{ totalstackedcount=0 }{ if($4 > 1) totalstackedcount += (( $3 - $2) * $4) }END{ print totalstackedcount }'

            #### Number of nucleotides with stacked crosslink events, i.e. positions with >1 crosslink events: 
            # cat <sampleX.plus.bedgraph> <sampleX.minus.bedgraph> | awk 'BEGIN{ totalstackedpos=0 }{ if($4 > 1) totalstackedpos += ($3 - $2) }END{ print totalstackedpos }' 
            cat $bedgraph_p_file $bedgraph_n_file  | awk 'BEGIN{ totalstackedpos=0 }{ if($4 > 1) totalstackedpos += ($3 - $2) }END{ print totalstackedpos }' 

            #### Convert bam to sam file format (using SAMtools)
            # samtools view <sampleX.duprm.bam> -o <sampleX.duprm.sam>
            samtools view $deduped_file -o $deduped_sam_file

            #### Number of reads mapped with deletions:
            # cut -f6 <sampleX.duprm.sam> | grep D | wc -l
            cut -f6 $deduped_sam_file | grep D | wc -l

            #### Number of reads mapped with insertions:
            # cut -f6 <sampleX.duprm.sam> | grep I | wc -l
            cut -f6 $deduped_sam_file | grep I | wc -l
        fi
    done
fi

# run iCLIPro
if [[ $run_pro == "Y" ]]; then
    sub_dir=$parent_dir/07_pro
    if [[ ! -f $sub_dir ]]; then mkdir -p $sub_dir; fi

    for sample_id in "${sample_list[@]}"; do
        deduped_file=$parent_dir/04_dedup/$sample_id.deduped.bam
        pro_file=$sub_dir/$sample_id.something.bam
            
        if [[ ! -f $pro_file ]]; then
            echo "--$sample_id"
            
            #### iCLIPro analysis and plots
            # iCLIPro -r 50 -b 300 -f 50 \
            # -g "L15:15,L16:16,L17:17,L18:18,L19:19,L20:20,L21:21,L22:22,L23:23,L24:24,L25:25,L26:26,L27:27,L28:28,L29:29,L30:30,L31:31,L32:32,L33:33,L34:34,L35:35,L36:36,L37:37,L38:38,L39:39,L40:40,R:41" \
            # -p "L15-R,L16-R,L17-R,L18-R,L19-R,L20-R,L21-R,L22-R,L23-R,L24-R,L25-R,L26-R,L27-R,L28-R,L29-R,L30-R,L31-R,L32-R,L33-R,L34-R,L35-R,L36-R,L37-R,L38-R,L39-R,L40-R" \
            # -o <outdir> <sampleX.duprm.bam>
            iCLIPro -r 50 -b 300 -f 50 \
            -g "L15:15,L16:16,L17:17,L18:18,L19:19,L20:20,L21:21,L22:22,L23:23,L24:24,L25:25,L26:26,L27:27,L28:28,L29:29,L30:30,L31:31,L32:32,L33:33,L34:34,L35:35,L36:36,L37:37,L38:38,L39:39,L40:40,R:41" \
            -p "L15-R,L16-R,L17-R,L18-R,L19-R,L20-R,L21-R,L22-R,L23-R,L24-R,L25-R,L26-R,L27-R,L28-R,L29-R,L30-R,L31-R,L32-R,L33-R,L34-R,L35-R,L36-R,L37-R,L38-R,L39-R,L40-R" \
            -o $sub_dir $deduped_file
        fi
    done
fi

# run PureCLIP
if [[ $run_pure == "Y" ]]; then
    sub_dir=$parent_dir/08_pure
    if [[ ! -f $sub_dir ]]; then mkdir -p $sub_dir; fi

    for sample_id in "${sample_list[@]}"; do
        pro_file=$parent_dir/07_pro/$sample_id.something.bam
            
        if [[ ! -f $pro_file ]]; then
            echo "--$sample_id"
            module load samtools pureclip
            #### Merge BAM files (sampleX.duprm.bam)
            #samtools merge -f <merged.bam> -b <list_of_bam_files>
            # samtools index <merged.bam>
            samtools index $pro_file

            #### Run PureCLIP
            # pureclip -i $pro_file -bai $pro_file.bai -g <genome.fasta> -ld -nt 8 -o PureCLIP.crosslink_sites.bed -or PureCLIP.crosslink_regions.bed
            cmd=`pureclip_command`
            echo $cmd
        fi
    done
fi