
#!/bin/sh

#set dir
log_dir=$1
project_dir=$2
doc="/data/CCBR_Pipeliner/iCLIP/container/USeq_8.9.6/Apps/SamTranscriptomeParser"

#set args from command line
rulename=$3
sample_name=$4
maxn=$5

#mkdirs if necessary
if [[ ! -d "$log_dir" ]]; then mkdir $log_dir; fi

#set alignment types
declare -a a_types=("masked" "unmasked" "unaware")

##############################################################################################################################
#Run cleanup
input_dir="${project_dir}/01_preprocess/02_alignment"
output_dir_m="${project_dir}/01_preprocess/03_genomic"
output_dir_u="${project_dir}/01_preprocess/04_unmapped"

if [[ ! -d "$output_dir_m" ]]; then mkdir $output_dir_m; fi
if [[ ! -d "$output_dir_u" ]]; then mkdir $output_dir_u; fi

#Create SH
if [ "$rulename" == "cleanup" ]; then \
    
    #for each sample
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            #if the alignment file exists, BUT the cleanup file doesnt exist then create an SH to create it
            if [[ -f $input_dir/${sample_name}.${a_id}.split.$i.sam.gz ]] && [[ ! -f $output_dir_m/${sample_name}.${a_id}.split.$i.sam.gz ]] ; then
                
                #SH file name
                sh_file="$log_dir/cleanup_${sample_name}.${a_id}.split.${i}_sh.sh"
                if [ -f $sh_file ]; then rm $sh_file; fi

                #set file names
                input_file="$input_dir/${sample_name}.${a_id}.split.${i}.sam.gz"
                base="${sample_name}.${a_id}.split.${i}"
                output_file_m="$output_dir_m/${sample_name}.${a_id}.split.${i}.sam"
                output_file_u="$output_dir_u/${sample_name}.${a_id}.split.${i}.final.si.bam"

                echo "#!/bin/sh
                module load samtools java;
                
                #set / create tmp dir
                tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"            

                if [ ${a_id} == "unaware" ]; then
                    cp ${input_file} ${output_file_m};
                    gzip ${output_file_m};
                    touch ${output_file_u};
                else
                    #run cleanup
                    zcat ${input_file} | samtools view | awk '{{ if ((\$4 == 1 && \$6 ~/(^[0-9]+)H|S([0-9]+I)/ && \$1~/:/ ) || (\$4 == 1 && \$6 ~/^[0-9]I/ && \$1~/:/ )) {{}} else  print }}' > \${tmpdir}/$base.tmp.sam
                    zcat ${input_file} | samtools view -H | cat - \${tmpdir}/$base.tmp.sam > \${tmpdir}/$base.tmp.mapped.sam
                    zcat ${input_file} | samtools view -f 4 > \${tmpdir}/$base.tmp.unmapped.sam

                    #run genomic conversion
                    java -Djava.io.tmpdir=\${tmpdir} -jar -Xmx100G $doc \
                    -f \${tmpdir}/$base.tmp.mapped.sam \
                    -a 50000 -n 25 -u -s $output_file_m
                    
                    #converts sam to bam
                    gunzip -c $output_file_m.gz | samtools view -S -b - > \${tmpdir}/$base.tmp.unmasked.bam
                    
                    #merge unmapped sam and unmasked bam into final output
                    FILESIZE=\$(stat -c%s \${tmpdir}/$base.tmp.unmapped.sam)
                    
                    #if there are unmapped reads
                    if [ \"\$FILESIZE\" -gt 100 ]; then
                        #merge unmasked reads with unmapped reads
                        samtools merge -f \${tmpdir}/$base.merged.bam \${tmpdir}/$base.tmp.unmapped.sam \${tmpdir}/$base.tmp.unmasked.bam
                        
                        #pull header from merged and add to final output
                        samtools view -h -o \${tmpdir}/$base.final.bam \${tmpdir}/$base.merged.bam;
                        
                        #sort and index final outupt
                        samtools sort -T \${tmpdir} \${tmpdir}/$base.final.bam -o $output_file_u;
                        samtools index $output_file_u

                    #if there are no unmapped reads then final file is unmasked bam
                    else
                        samtools sort -T \${tmpdir} \${tmpdir}/$base.tmp.unmasked.bam -o $output_file_u
                        samtools index $output_file_u
                    fi
                fi
                " > $sh_file
            fi
        done
    done
fi

##############################################################################################################################
#Run unmapped merge
input_dir="${project_dir}/01_preprocess/04_unmapped"
output_dir="${project_dir}/02_bam/01_unmapped"

if [[ ! -d "$output_dir" ]]; then mkdir -p $output_dir; fi

#Create SH
if [ "$rulename" == "unmapped" ]; then \
    #set alignment
    a_id="unmasked"

    #set output
    output_file="$output_dir/${sample_name}.${a_id}.complete.bam"                
        
    #if the unmapped sam file exists, BUT the unmapped bam file doesnt exist then create an SH to create it
    if [[ -f $input_dir/${sample_name}.${a_id}.split.$i.final.si.bam ]] && [[ ! -f $output_file ]] ; then

        #SH file name
        sh_file="$log_dir/unmapped_${sample_name}.${a_id}_sh.sh"
        if [ -f $sh_file ]; then rm $sh_file; fi
        
        echo "#!/bin/sh
        module load samtools samtools;
                    
        #set / create tmp dir
        tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"    

        samtools merge -f $output_file \\" > $sh_file

        for ((n = 1; n <= $maxn; n++ )); do
                
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi
            
            #echo list of files
            if [ "${n}" -eq $maxn ]; then
                echo "$input_dir/${sample_name}.${a_id}.split.$i.final.si.bam " >> $sh_file
            else
                echo "$input_dir/${sample_name}.${a_id}.split.$i.final.si.bam \\" >> $sh_file
            fi
        done
    fi
fi

##############################################################################################################################
#Run create Unique and MM
input_dir="${project_dir}/01_preprocess/03_genomic"
output_dir_u="${project_dir}/01_preprocess/05_unique"
output_dir_m="${project_dir}/01_preprocess/05_mm"

if [[ ! -d "$output_dir_u" ]]; then mkdir $output_dir_u; fi
if [[ ! -d "$output_dir_m" ]]; then mkdir $output_dir_m; fi

#run unique
if [ "$rulename" == "unique_mm" ]; then \
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            #if the unmapped sam file exists, BUT the unmapped bam file doesnt exist then create an SH to create it
            input_file=$input_dir/${sample_name}.${a_id}.split.$i.sam.gz
            output_file="$output_dir_u/${sample_name}.${a_id}.split.${i}.unique.si.bam"
            if [[ -f $input_file ]] && [[ ! -f $output_file ]] ; then
            
                #create sh file
                sh_file="$log_dir/unique_${sample_name}.${a_id}.split.${i}_sh.sh"

                #set files, base, name
                f1="$input_dir/${sample_name}.${a_id}.split.${i}.sam.gz"
                p_base="${sample_name}.${a_id}.split.${i}.tmp"

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
                samtools sort -T \${tmpdir_u} \"\${tmpdir_u}/${p_base}.unique.bam\" -o $output_file;
                samtools index $output_file;

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

#run mm
if [ "$rulename" == "unique_mm" ]; then \
    for a_id in ${a_types[@]}; do
        for ((n = 1; n <= $maxn; n++ )); do
            
            #pad numbers 1-9 with a "0" to match naming schema of snakemake
            if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

            #if the unmapped sam file exists, BUT the unmapped bam file doesnt exist then create an SH to create it
            input_file=$input_dir/${sample_name}.${a_id}.split.$i.sam.gz
            output_file="$output_dir_m/${sample_name}.${a_id}.split.${i}.mm.si.bam"
            if [[ -f $input_file ]] && [[ ! -f $output_file ]] ; then
            
                #create sh file
                sh_file="$log_dir/mm_${sample_name}.${a_id}.split.${i}_sh.sh"

                #set file, base, names
                f1="$input_dir/${sample_name}.${a_id}.split.${i}.sam.gz"
                p_base="${sample_name}.${a_id}.split.${i}.tmp"

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
                samtools sort -T \${tmpdir_m} \"\${tmpdir_m}/${p_base}.mm.bam\" -o $output_file;
                samtools index $output_file;

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
input_dir_u="${project_dir}/01_preprocess/05_unique"
input_dir_m="${project_dir}/01_preprocess/05_mm"
output_dir="${project_dir}/02_bam/02_merged"

if [[ ! -d "$output_dir" ]]; then mkdir $output_dir; fi

if [ "$rulename" == "merge_splits" ]; then \
    for a_id in ${a_types[@]}; do
        
        #run unique
        #output file doesnt exist then create an SH to create it
        output_file="$output_dir_u/${sample_name}.${a_id}.merged.unique.si.bam"
        if [[ ! -f $output_file ]] ; then
            
            sh_file="$log_dir/merge_splits_u_${sample_name}_${a_id}_sh.sh"
            if [ -f $sh_file ]; then rm $sh_file; fi
            touch $sh_file
        
            echo "#!/bin/bash
            module load samtools
            samtools merge -f $output_file \\" > $sh_file    
        
            for ((n = 01; n <= $maxn; n++ )); do

                #pad numbers 1-9 with a "0" to match naming schema of snakemake
                if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi

                input_file="$input_dir_u/${sample_name}.${a_id}.split.${i}.unique.si.bam"

                #echo list of files
                if [ "${n}" -eq $maxn ]; then
                    echo "$input_file;" >> $sh_file
                else
                    echo "$input_file \\" >> $sh_file
                fi
            done
        fi

        #run mm
        #output file doesnt exist then create an SH to create it
        output_file="$output_dir_m/${sample_name}.${a_id}.merged.mm.si.bam"
        if [[ ! -f $output_file ]] ; then
            
            sh_file="$log_dir/merge_splits_m_${sample_name}_${a_id}_sh.sh"
            if [ -f $sh_file ]; then rm $sh_file; fi
            touch $sh_file
        
            echo "#!/bin/bash
            module load samtools
            samtools merge -f $output_file \\" > $sh_file    
        
            for ((n = 1; n <= $maxn; n++ )); do

                #pad numbers 1-9 with a "0" to match naming schema of snakemake
                if [ "${n}" -lt 10 ]; then i="0$n"; else i=$n; fi
                
                input_file="$input_dir_m/${sample_name}.${a_id}.split.${i}.mm.si.bam"
                
                #echo list of files
                if [ "${n}" -eq $maxn ]; then
                    echo "$input_file;" >> $sh_file
                else
                    echo "$input_file \\" >> $sh_file
                fi
            done
        fi
    done
fi


##############################################################################################################################
#Run merged unique and mm
input_dir="${project_dir}/02_bam/02_merged"
output_dir="${project_dir}/02_bam/02_merged"
output_dir_qc="${project_dir}/qc/01_qc_post"

if [[ ! -d "$output_dir" ]]; then mkdir -p $output_dir; fi
if [[ ! -d "$output_dir_qc" ]]; then mkdir -p $output_dir_qc; fi

if [ "$rulename" == "merge_um" ]; then \
    for a_id in ${a_types[@]}; do

        #output file doesnt exist then create an SH to create it
        output_file="$output_dir/${sample_name}.${a_id}.merged.si.bam"
        output_file_st="$output_dir_qc/${sample_name}.${a_id}_samstats.txt"
        if [[ ! -f $output_file ]] || [[ ! -f $output_file_st ]]; then
        
            #set sh files
            sh_file="$log_dir/merge_um_${sample_name}_${a_id}_sh.sh"
            if [ -f $sh_file ]; then rm $sh_file; fi
            
            echo "#!/bin/bash
            module load samtools
            
            #set / create tmp dir
            tmpdir=\"/lscratch/\${SLURM_JOB_ID}\"
            
            #merge uniuqe and mm
            samtools merge -f \${tmpdir}/tmp.bam $input_dir/${sample_name}.${a_id}.merged.unique.bam $input_dir/${sample_name}.${aid}.merged.mm.bam 
                
            #sort and index
            samtools sort -T \${tmpdir} \${tmpdir}/tmp.bam -o $output_file
            samtools index $output_file

            #run samstats
            samtools view -h $output_file | samtools stats - > $output_file_st
            " > $sh_file
        fi
    done
fi