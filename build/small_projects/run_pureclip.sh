################################
# setup
################################
# set dirs
pureclip_dir="/data/CCBR_Pipeliner/bin/PureCLIP/PureCLIP/build/pureclip"
output_dir="/data/RBL_NCI/Wolin/Sam/pureclip"
if [[ ! -d $output_dir ]]; then mkdir $output_dir; fi

# set samples, dirs
sample_list=("Ro_Clip" "Control_Clip")
v_dir_list=("/data/sevillas2/test_2/v1/02_bam/03_dedup/" "/data/sevillas2/test_2/v4/02_bam/02_dedup")

# controls
run_pureclip="N"
run_purecomparison="Y"

################################
# functions
################################
function pureclip_command()
{
    local pc_cmd="$pureclip_dir \
        -i $merged_bam \
        -bai ${merged_bam}.bai \
        -g $clean_reference \
        -ld -nt 8 \
        -o ${output_dir}/${sample_id}_${v_id}_sites.bed \
        -or ${output_dir}/${sample_id}_${v_id}_regions.bed"
    echo "$pc_cmd"
}

function pureclip_comparison_command()
{
    v1_file=${output_dir}/${sample_id}_v1_${file_type}.bed
    v2_file=${output_dir}/${sample_id}_v2_${file_type}.bed
    v1_sorted=${output_dir}/comparisons/${sample_id}_v1_${file_type}_sorted.bed
    v2_sorted=${output_dir}/comparisons/${sample_id}_v2_${file_type}_sorted.bed

    # sort files
    if [[ ! -f $v1_sorted ]]; then bedtools sort -i $v1_file > $v1_sorted; fi
    if [[ ! -f $v2_sorted ]]; then bedtools sort -i $v2_file > $v2_sorted; fi

    # https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
    # in V1 but not in V2
    bedtools intersect -v -a $v1_sorted -b $v2_sorted > ${output_dir}/comparisons/${sample_id}_notinV2_${file_type}.txt
    
    # in V2 but not in V2
    bedtools intersect -v -a $v2_sorted -b $v1_sorted > ${output_dir}/comparisons/${sample_id}_notinV1_${file_type}.txt

    # intersection of both
    bedtools jaccard -a $v1_sorted -b $v2_sorted > ${output_dir}/comparisons/${sample_id}_inboth_${file_type}.txt
}
################################
# commands
################################

# clean reference file (Y's are in the chrY)
reference_file="/data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/mm10/GRCm38.p6.genome.fa"
clean_reference="/data/RBL_NCI/Wolin/Sam/pureclip/GRCm38.p6.genome.fa"
if [[ ! -f $clean_reference ]]; then
    echo "** Creating REF **"
    cat $reference_file | sed '/^\>/! s/[RYKMSWBVHD]/N/g' > $clean_reference
fi

# run pureclip on RBL samples
if [[ $run_pureclip == "Y" ]]; then
    echo "** Running PURECLIP **"
    for v_dir in "${v_dir_list[@]}"; do
        v_id=`echo $v_dir | cut -f5 -d"/" | sed "s/4/2/g"`
        echo "--$v_id"
        for sample_id in "${sample_list[@]}"; do
            echo "----$sample_id"
            merged_bam="${v_dir}/${sample_id}.dedup.si.bam"
            cmd=`pureclip_command`
            $cmd
        done
    done
fi

# run comparison of pureclip output on RBL samples
if [[ $run_purecomparison == "Y" ]]; then
    module load bedtools
    echo "** RUNNING COMPARISON **"
    for sample_id in "${sample_list[@]}"; do
        echo "--$sample_id"
        
        ##################
        ### SITES
        file_type="sites"
        pureclip_comparison_command

        ##################
        ### REGIONS
        file_type="regions"
        pureclip_comparison_command
    done

    echo "** RESULTS **"
    for sample_id in "${sample_list[@]}"; do
        echo "--$sample_id"
        regions=${output_dir}/comparisons/${sample_id}_inboth_regions.txt
        sites=${output_dir}/comparisons/${sample_id}_inboth_sites.txt
        
        echo "----regions"
        cat $regions

        echo "----sites"
        cat $sites
        echo ""
    done
fi

