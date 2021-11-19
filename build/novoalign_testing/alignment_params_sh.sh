#!/bin/sh

#sh /home/sevillas2/git/iCLIP/build/novoalign_testing/2021_10_v2/alignment_params_sh.sh

#set version
version=$1
if [[ $version == "v1" ]]; then 
    log_dir="/data/RBL_NCI/Wolin/Sam/novoalign/log"
    project_dir="/data/RBL_NCI/Wolin/Sam/novoalign"
else
    log_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2/log"
    project_dir="/data/RBL_NCI/Wolin/Sam/novoalign_v2"
fi

#set option
option=$2

index_file="/data/CCBR_Pipeliner/iCLIP/index/active/phil/mm10/mm10_splice75bp_unmasked.nix"
doc="/data/CCBR_Pipeliner/iCLIP/container/USeq_8.9.6/Apps/SamTranscriptomeParser"
alias_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/hg38.chromAlias.txt"
gencode_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/Gencode_V32/fromGencode/gencode.v32.annotation.gtf.txt"
refseq_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/NCBI_RefSeq/GCF_000001405.39_GRCh38.p13_genomic.gtf.txt"
canonical_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/Gencode_V32/fromUCSC/KnownCanonical/KnownCanonical_GencodeM32_GRCh38.txt"
intron_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/Gencode_V32/fromUCSC/KnownGene/KnownGene_GencodeV32_GRCh38_introns.bed"
rmsk_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/repeatmasker/rmsk_GRCh38.txt"
custom_path="/data/CCBR_Pipeliner/iCLIP/ref/annotations/hg38/additional_anno/"

#mkdirs if necessary
if [[ ! -d "$log_dir" ]]; then mkdir $log_dir; fi
for test_id in {1..9}; do \
    if [[ ! -d "$log_dir/test$test_id/" ]]; then mkdir "$log_dir/test$test_id/"; fi
done

##############################################################################################################################
#Run splits
if [ "$option" == "splits" ]; then \
    zcat /data/RBL_NCI/Wolin/Phil/novoalign/01_preprocess/01_NOyRNA/Ro_Clip_2_filtered.NOyRNA.fastq.gz | \
    split --additional-suffix .fastq -l 1906896 --numeric-suffixes=1 --filter='gzip > $FILE.gz' - ${project_dir}/splits/Ro_Clip_2.split.
fi

##############################################################################################################################
#Run alignment
input_dir="${project_dir}/splits"
output_dir="${project_dir}/aligned"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
for test_id in {1..9}; do \
    if [[ ! -d "$output_dir/test$test_id/" ]]; then mkdir "$output_dir/test$test_id/"; fi
done

if [ "$option" == "align" ]; then \
    # test1: iCLIP default
    test_id="test1"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -k \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -s 1 \
            -o SAM \
            -R 0 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"
            " > $sh_file
    done

    # test2: R flag
    test_id="test2"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -k \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -s 1 \
            -o SAM \
            -R 5 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test3: R flag
    test_id="test3"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -k \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -s 1 \
            -o SAM \
            -R 20 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test4: STRICTER VALUES
    test_id="test4"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
        novoalign -d $index_file \
        -k \
        -f - \
        -F STDFQ \
        -c 32 \
        -t 20,3 \
        -l 20 \
        -x 4 \
        -g 30 \
        -s 1 \
        -o SAM \
        -R 5 \
        -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test5: NOVOALIGN default
    test_id="test5"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
        novoalign -d $index_file \
        -k \
        -f - \
        -F STDFQ \
        -c 32 \
        -t 20,3 \
        -x 6 \
        -g 40 \
        -s 2 \
        -o SAM \
        -R 5 \
        -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test 6 S FLAG
    test_id="test6"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -k \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -s 0 \
            -o SAM \
            -R 0 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test 7 S FLAG
    test_id="test7"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -k \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -s 2 \
            -o SAM \
            -R 0 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test 8 S FLAG
    test_id="test8"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -k \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -o SAM \
            -R 0 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

    # test 9 K FLAG
    test_id="test9"
    for i in {01..10}; do \
        input_file="$input_dir/Ro_Clip_2.split.$i.fastq.gz"
        split_base="$output_dir/$test_id/Ro_Clip_2.split.$i"
        sh_file="$log_dir/${test_id}/align_${test_id}_${i}_sh.sh"
        
        echo "#!/bin/sh
        module load novocraft

        zcat $input_file | \
            novoalign -d $index_file \
            -f - \
            -F STDFQ \
            -c 32 \
            -t 15,3 \
            -l 20 \
            -x 4 \
            -g 20 \
            -s 1 \
            -o SAM \
            -R 0 \
            -r EXHAUSTIVE 999 | gzip -c > "${split_base}.sam.gz"" > $sh_file
    done

fi

##############################################################################################################################
#Run cleanup
input_dir="${project_dir}/aligned"
output_dir="${project_dir}/cleanup"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
for test_id in {1..9}; do \
    if [[ ! -d "$output_dir/test$test_id/" ]]; then mkdir "$output_dir/test$test_id/"; fi
done
if [ "$option" == "cleanup" ]; then \
    for test_id in {1..9}; do
        for split_id in {01..10}; do
            sh_file="$log_dir/test${test_id}/cleanup_test${test_id}_${split_id}_sh.sh"
            if [ -f $sh_file ]; then rm $sh_file; fi

            input_file="$input_dir/test$test_id/Ro_Clip_2.split.$split_id.sam.gz"
            base="$test_id.split.$split_id"
            output_file="$output_dir/test$test_id/Ro_Clip_2.split.$split_id.sam"
            
            echo "#!/bin/sh
            module load samtools java;
                        
            #set / create tmp dir
            if [ -d "/lscratch/\${SLURM_JOB_ID}" ]; then
                tmpdir="/lscratch/\${SLURM_JOB_ID}"
            fi
   
            #run cleanup
            zcat $input_file | samtools view | awk '{ if ((\$4 == 1 && \$6!~/^[0-9]I/ && \$1~/:/ )||(\$4 > 1 && \$1~/:/ )) { print } }' > \${tmpdir}/$base.tmp.sam
            zcat $input_file | samtools view -H | cat - \${tmpdir}/$base.tmp.sam > \${tmpdir}/$base.tmp.mapped.sam

            #run genomic conversion
            java -Djava.io.tmpdir=\${tmpdir} -jar -Xmx100G $doc \
                -f \${tmpdir}/$base.tmp.mapped.sam \
                -a 50000 -n 25 -u -s $output_file
          " > $sh_file
        done
    done
fi

##############################################################################################################################
#Run create Unique and MM
input_dir="${project_dir}/cleanup"
output_dir="${project_dir}/unique_mm"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
for test_id in {1..9}; do \
    if [[ ! -d "$output_dir/test$test_id/" ]]; then mkdir "$output_dir/test$test_id/"; fi
done
if [ "$option" == "unique_mm" ]; then \

    for test_id in {1..9}; do
        for split_id in {01..10}; do
            
            f1="$input_dir/test$test_id/Ro_Clip_2.split.$split_id.sam.gz"
            p_base="test$test_id.split.$split_id.tmp"
            sort_u="$output_dir/test$test_id/test$test_id.split.$split_id.unique.si.bam"
            sh_file="$log_dir/test${test_id}/unique_test${test_id}_${split_id}_sh.sh"

            if [[ ! -f "$sort_u" ]]; then
                echo "#!/bin/sh
                module load samtools;

                ################################################
                #set / create tmp dir
                ################################################
                tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"
                tmpdir_m=\"\${tmpdir}/m\"
                tmpdir_u=\"\${tmpdir}/u\"

                mkdir \$tmpdir_m
                mkdir \$tmpdir_u
                ################################################
                #create mm, unique, header files
                ################################################
                set +e
                gunzip -c $f1 | samtools view | grep -v 'IH:i' > \"\${tmpdir_u}/${p_base}.unique.txt\";
                gunzip -c $f1 | samtools view -H > \"\${tmpdir}/${p_base}.header.txt\"
                
                ################################################
                #create unique bam, correcting NH header
                ################################################
                #if the file exists, remove it
                if [ -f \"\${tmpdir_u}/${p_base}.unique.sam\" ]; then rm \"\${tmpdir_u}/${p_base}.unique.sam\"; fi

                #create the new files
                touch \"\${tmpdir_u}/${p_base}.unique.sam\"

                #read the sam file
                while read rd; do 
                    #set tab variable
                    T=\$(printf \"\t\");
                    
                    #if the line doesnt have NH then add NH:i:1
                    if [[ ! \$rd == \"NH:i\" ]]; then
                        echo \"\$rd\" | awk '{ \$(NF+1) = \"NH:i:1\"; print }' | sed \"s/[[:blank:]]\+/\$T/g\" >> \"\${tmpdir_u}/${p_base}.unique.sam\"

                    #if it does have NH then replace it with NH:i:1
                    else 
                        replacement=\"NH:i:1\"
                        fixed=\${rd/NH:i:[0-9]*/\$replacement}
                        echo \"\$fixed\" >> \"\${tmpdir_u}/${p_base}.unique.sam\"
                    fi
                done < \"\${tmpdir_u}/${p_base}.unique.txt\"

                #place header back on file, sort, convert to bam
                cat \"\${tmpdir}/${p_base}.header.txt\" \"\${tmpdir_u}/${p_base}.unique.sam\" | samtools sort -T \${tmpdir_u} | samtools view -Sb > \"\${tmpdir_u}/${p_base}.unique.bam\";

                ################################################
                #sort unique and mm bam
                ################################################
                samtools sort -T \${tmpdir_u} \"\${tmpdir_u}/${p_base}.unique.bam\" -o $sort_u;
                samtools index $sort_u;

                exitcode=\$?
                if [ \$exitcode -eq 1 ]; then
                    exit 0
                else
                    exit 0
                fi
                " > $sh_file
            fi
        done
    done
fi

if [ "$option" == "unique_mm" ]; then \
    for test_id in {1..9}; do
        for split_id in {01..10}; do
            
            f1="$input_dir/test$test_id/Ro_Clip_2.split.$split_id.sam.gz"
            p_base="test$test_id.split.$split_id.tmp"
            sort_m="$output_dir/test$test_id/test$test_id.split.$split_id.mm.si.bam"

            sh_file="$log_dir/test${test_id}/mm_test${test_id}_${split_id}_sh.sh"

            if [[ ! -f "$sort_m" ]]; then
                echo "#!/bin/sh
                module load samtools;

                ################################################
                #set / create tmp dir
                ################################################
                tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"
                tmpdir_m=\"\${tmpdir}/m\"
                tmpdir_u=\"\${tmpdir}/u\"

                mkdir \$tmpdir_m
                mkdir \$tmpdir_u
                ################################################
                #create mm, unique, header files
                ################################################
                set +e
                gunzip -c $f1 | samtools view | grep 'IH:i' > \"\${tmpdir_m}/${p_base}.mm.txt\";
                gunzip -c $f1 | samtools view -H > \"\${tmpdir}/${p_base}.header.txt\"
                
                ################################################
                #create mm bam, correcting NH header
                ################################################
                #if the file exists, remove it
                if [ -f \"\${tmpdir_m}/${p_base}.mm.sam\" ]; then rm \"\${tmpdir_m}/${p_base}.mm.sam\"; fi

                #create the new files
                touch \"\${tmpdir_m}/${p_base}.mm.sam\"

                #read the sam file
                while read rd; do 
                    #replace the NH value with the IH value
                    IH=\`echo \$rd | grep -o \"IH:i:[0-9]*\"\`
                    IH_num=\`echo \${IH/IH:i:/}\`
                    
                    NH=\`echo \$rd | grep -o \"NH:i:.*\S\"\`
                    NH_num=\`echo \${NH/NH:i:/}\`
                    
                    NH_complete=\`echo \${NH/\$NH_num/\$IH_num}\`
                    
                    echo \"\${rd/\$NH/\$NH_complete}\" >> \"\${tmpdir_m}/${p_base}.mm.sam\"
                done < \"\${tmpdir_m}/${p_base}.mm.txt\"

                #place header back on file, sort, convert to bam
                cat \"\${tmpdir}/${p_base}.header.txt\" \"\${tmpdir_m}/${p_base}.mm.sam\" | samtools sort -T \${tmpdir_m} | samtools view -Sb > \"\${tmpdir_m}/${p_base}.mm.bam\";

                ################################################
                #sort unique and mm bam
                ################################################
                samtools sort -T \${tmpdir_m} \"\${tmpdir_m}/${p_base}.mm.bam\" -o $sort_m;
                samtools index $sort_m;

                exitcode=\$?
                if [ \$exitcode -eq 1 ]; then
                    exit 0
                else
                    exit 0
                fi
                " > $sh_file
            fi
        done
    done
fi

##############################################################################################################################
#Run merged splits
input_dir="${project_dir}/unique_mm"
output_dir="${project_dir}/merge_splits"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
for test_id in {1..9}; do \
    if [[ ! -d "$output_dir/test$test_id/" ]]; then mkdir "$output_dir/test$test_id/"; fi
done

if [ "$option" == "merge_splits" ]; then \
    for test_id in {1..9}; do
        sh_file="$log_dir/test${test_id}/merge_splits_test${test_id}_sh.sh"
        if [ -f $sh_file ]; then rm $sh_file; fi
        touch $sh_file

        echo "#!/bin/bash
        module load samtools;" >> $sh_file    
        
        #run unique
        echo "samtools merge -f $output_dir/test${test_id}/test${test_id}.merged.unique.bam \\" >> $sh_file
        for split_id in {01..10}; do
            if [ "$split_id" == "10" ]; then 
                echo "$input_dir/test$test_id/test$test_id.split.$split_id.unique.si.bam;" >> $sh_file
            else
                echo "$input_dir/test$test_id/test$test_id.split.$split_id.unique.si.bam \\" >> $sh_file
            fi
        done
        echo >> $sh_file     

        #run mm
        echo "samtools merge -f $output_dir/test${test_id}/test${test_id}.merged.mm.bam \\" >> $sh_file
        for split_id in {01..10}; do
            if [ "$split_id" == "10" ]; then 
                echo "$input_dir/test$test_id/test$test_id.split.$split_id.mm.si.bam;" >> $sh_file
            else
                echo "$input_dir/test$test_id/test$test_id.split.$split_id.mm.si.bam \\" >> $sh_file
            fi
        done
        echo >> $sh_file
    done
fi


##############################################################################################################################
#Run merged unique and mm
input_dir="${project_dir}/merge_splits"
output_dir="${project_dir}/merged"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi

if [ "$option" == "merge_um" ]; then \
    for test_id in {1..9}; do \
        sh_file="$log_dir/test${test_id}/merge_um_test${test_id}_sh.sh"
        if [ -f $sh_file ]; then rm $sh_file; fi

        un="$input_dir/test${test_id}/test${test_id}.merged.unique.bam"
        mm="$input_dir/test${test_id}/test${test_id}.merged.mm.bam"

        output_file="$output_dir/test${test_id}.merged.si.bam"
        
        echo "#!/bin/bash
        module load samtools
        
        #set / create tmp dir
        tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"
        
        #merge uniuqe and mm
        samtools merge -f \${tmpdir}/tmp.bam $un $mm
            
        #sort and index
        samtools sort -T \${tmpdir} \${tmpdir}/tmp.bam -o $output_file
        samtools index $output_file
        " > $sh_file
    done
fi

##############################################################################################################################
#Run dedup
input_dir="${project_dir}/merged"
output_dir="${project_dir}/dedup"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi

if [ "$option" == "dedup" ]; then \
    for test_id in {1..9}; do \
        sh_file="$log_dir/test${test_id}/dedup_test${test_id}_sh.sh"

        f1="$input_dir/test${test_id}.merged.si.bam"
        base="test${test_id}.dedup.si"
        output_file="$output_dir/test${test_id}.dedup.si.bam"

        echo "#!/bin/bash
        module load samtools umitools
        
        #set / create tmp dir
        tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"

        #dedup
        umi_tools dedup \
        -I $f1 \
        --method unique \
        --multimapping-detection-method=NH \
        --umi-separator="rbc:" \
        -S \${tmpdir}/$base.bam\
        --log2stderr;

        #sort and index
        samtools sort -T  \${tmpdir} \${tmpdir}/$base.bam -o $output_file
        samtools index $output_file
        " > $sh_file
    done
fi

# ##############################################################################################################################
#Run create_beds_safs
input_dir="${project_dir}/dedup"
output_dir="${project_dir}/beds"
bed_dir="$output_dir/01_bed"
saf_dir="$output_dir/02_SAF"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
if [[ ! -d "$bed_dir" ]]; then mkdir $bed_dir; fi
if [[ ! -d "$saf_dir" ]]; then mkdir $saf_dir; fi

if [ "$option" == "beds" ]; then
    for test_id in {1..9}; do
        sh_file="$log_dir/test${test_id}/beds_test${test_id}_sh.sh"
        input_file="$input_dir/test${test_id}.dedup.si.bam"
        output_bed="$bed_dir/test${test_id}"
        output_saf="$saf_dir/test${test_id}"
        base="test${test_id}.dedup"

        echo "#!/bin/bash
        module load samtools bedtools

        #set / create tmp dir
        tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"

        #create header for unique
        samtools view -H $input_file > \${tmpdir}/$base.header.txt;
        
        #create unique bam
        samtools view $input_file | grep -w 'NH:i:1' | cat \${tmpdir}/$base.header.txt - |  samtools sort -T \${tmpdir} | samtools view -Sb > \${tmpdir}/$base.unique.bam;
        
        #index
        cp \${tmpdir}/$base.unique.bam \${tmpdir}/$base.unique.i.bam; 
        samtools index \${tmpdir}/$base.unique.i.bam;
        
        #create SAFS
        bedtools bamtobed \
            -split -i ${input_file} | bedtools sort -i - > "${output_bed}_all.bed"; 
        bedtools bamtobed \
            -split -i \${tmpdir}/$base.unique.i.bam | bedtools sort -i - > "${output_bed}_unique.bed";
        bedtools merge \
            -c 6 -o count,distinct -bed -s -d 50 \
            -i "${output_bed}_all.bed"| \
            awk '{OFS=\"\t\"; print \$1\":\"\$2\"-\"\$3\"_\"\$5,\$1,\$2,\$3,\$5}'| \
            awk 'BEGIN{{print \"ID\",\"Chr\",\"Start\",\"End\",\"Strand\"}}1' > "${output_saf}_all.SAF"
        bedtools merge \
            -c 6 -o count,distinct -bed -s -d 50 \
            -i "${output_bed}_unique.bed" | \
            awk '{OFS=\"\t\"; print \$1\":\"\$2\"-\"\$3\"_\"\$5,\$1,\$2,\$3,\$5}'| \
            awk 'BEGIN{print \"ID\",\"Chr\",\"Start\",\"End\",\"Strand\"}1' > "${output_saf}_unique.SAF"
        " > $sh_file
    done
fi

# ##############################################################################################################################
#Run counts
input_dir_dedup="${project_dir}/dedup"
input_dir_saf="${project_dir}/beds/02_SAF"
output_dir="${project_dir}/counts"
all_dir="$output_dir/03_allreadpeaks"
unique_dir="$output_dir/03_uniquereadpeaks"

if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
if [[ ! -d "$output_dir/03_uniquereadpeaks" ]]; then mkdir "$output_dir/03_uniquereadpeaks"; fi
if [[ ! -d "$output_dir/03_allreadpeaks" ]]; then mkdir "$output_dir/03_allreadpeaks"; fi


if [ "$option" == "counts" ]; then
    for test_id in {1..9}; do

        sh_file="$log_dir/test${test_id}/counts_test${test_id}_sh.sh"
        if [ -f $sh_file ]; then rm $sh_file; fi

        input_file_dedup="$input_dir_dedup/test${test_id}.dedup.si.bam"
        input_file_all="$input_dir_saf/test${test_id}_all.SAF"
        input_file_unique="$input_dir_saf/test${test_id}_unique.SAF"

        echo "#!/bin/bash
        module load subread

        #run for allreadpeaks
        featureCounts -F SAF \
            -a $input_file_all \
            -O \
            -J \
            --fraction \
            --minOverlap 1 \
            -s 1 \
            -T 32 \
            -o "$all_dir/test${test_id}_uniqueCounts.txt" \
            $input_file_dedup;
        featureCounts -F SAF \
            -a $input_file_all \
            -M \
            -O \
            -J \
            --fraction \
            --minOverlap 1 \
            -s 1 \
            -T 32 \
            -o "$all_dir/test${test_id}_allFracMMCounts.txt" \
            $input_file_dedup
        #run for uniquereadpeaks
        featureCounts -F SAF \
            -a $input_file_unique \
            -O \
            -J \
            --fraction \
            --minOverlap 1 \
            -s 1 \
            -T 32 \
            -o "$unique_dir/test${test_id}_uniqueCounts.txt" \
            $input_file_dedup;
        featureCounts -F SAF \
            -a $input_file_unique \
            -M \
            -O \
            -J \
            --fraction \
            --minOverlap 1 \
            -s 1 \
            -T 32 \
            -o "$unique_dir/test${test_id}_allFracMMCounts.txt" \
            $input_file_dedup
        " > $sh_file
    done
fi

# ##############################################################################################################################
#Run project_annotations
script="/home/sevillas2/git/iCLIP/workflow/scripts/08_annotation.R"
config="/home/sevillas2/git/iCLIP/config/annotation_config.txt"
output_dir="${project_dir}/annotations"
if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi
if [[ ! -d "$output_dir/01_project" ]]; then mkdir "$output_dir/01_project"; fi

if [ "$option" == "project_annotations" ]; then
        module load R

        Rscript /home/sevillas2/git/iCLIP/workflow/scripts/08_annotation.R \
        --ref_species mm10 \
        --refseq_rRNA TRUE \
        --alias_path $alias_path \
        --gencode_path $gencode_path \
        --refseq_path $refseq_path \
        --canonical_path $canonical_path \
        --intron_path $intron_path \
        --rmsk_path $rmsk_path \
        --custom_path $custom_path \
        --out_dir "$output_dir/01_project/" \
        --reftable_path $config
fi

# ##############################################################################################################################
#Run peak_annotations
input_dir="${project_dir}/counts"
all_dir="$input_dir/03_allreadpeaks"
unique_dir="$input_dir/03_uniquereadpeaks/"

output_dir="${project_dir}/annotations"
if [[ ! -d "$output_dir/02_peaks" ]]; then mkdir "$output_dir/02_peaks"; fi


if [ "$option" == "peak_annotations" ]; then
    for test_id in {1..9}; do

        sh_file="$log_dir/test${test_id}/peak_annotations_test${test_id}_sh.sh"
        if [ -f $sh_file ]; then rm $sh_file; fi
        
        input_unique=""
        input_all=""

        echo "#!/bin/bash
        module load R bedtools

        #set / create tmp dir
        tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"

        Rscript /home/sevillas2/git/iCLIP/workflow/scripts/09_peak_annotation.R \
        --rscript /home/sevillas2/git/iCLIP/workflow/scripts/09_peak_annotation_functions.R \
        --peak_type ALL \
        --peak_unique "$all_dir/test${test_id}_uniqueCounts.txt" \
        --peak_all "$all_dir/test${test_id}_allFracMMCounts.txt" \
        --join_junction TRUE \
        --condense_exon TRUE \
        --read_depth 3 \
        --demethod NONE \
        --sample_id "test${test_id}" \
        --ref_species mm10 \
        --anno_dir "$output_dir/01_project/" \
        --reftable_path /data/RBL_NCI/Wolin/Sam/mapq_recalc/config/annotation_config.txt \
        --gencode_path $gencode_path \
        --intron_path $intron_path \
        --rmsk_path $rmsk_path \
        --tmp_dir \${tmpdir} \
        --out_dir "$output_dir/02_peaks/" \
        --out_dir_DEP "$output_dir/02_peaks/" \
        --output_file_error "$output_dir/02_peaks/read_depth_error.txt"
        " > $sh_file
    done
fi

##############################################################################################################################
#Run annotation_report
input_dir="${project_dir}/annotations/02_peaks"
output_dir="${project_dir}/annotations"
if [[ ! -d "$output_dir/03_reports" ]]; then mkdir "$output_dir/03_reports"; fi

if [ "$option" == "annotation_report" ]; then
    for test_id in {1..9}; do

        sh_file="$log_dir/test${test_id}/annotation_report_test${test_id}_sh.sh"
        if [ -f $sh_file ]; then rm $sh_file; fi

        echo "#!/bin/bash
        module load R
        
        Rscript -e 'library(rmarkdown); \
        rmarkdown::render(\"/home/sevillas2/git/iCLIP/workflow/scripts/10_annotation.Rmd\",
            output_file = \"$output_dir/03_reports/test${test_id}_annotation_final_report.html\",
            params= list(samplename = \"test${test_id}\",
            peak_in = \"$input_dir/test${test_id}_peakannotation_complete.txt\",
            output_table = \"$output_dir/test${test_id}_annotation_final_table.txt\",
            readdepth = \"0\",
            PeakIdnt = \"ALL\"))'
        " > $sh_file
    done
fi