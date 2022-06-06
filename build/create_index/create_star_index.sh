############################################################################################
# Controls
############################################################################################
flag_download="Y"
flag_index="Y"

############################################################################################
# USER INPUTS
############################################################################################
parent_dir="/data/CCBR_Pipeliner/iCLIP/index/active/2022_0505"
perm_bk_path="/data/CCBR_Pipeliner/iCLIP/index/active/2021_0607/mm10/01_source/BK000964.3_TPA_exp.fa"

# select mm10 or hg38
ref_species="hg38"

if [[ $ref_species == "mm10" ]]; then
    output_dir="$parent_dir/mm10"
    ref_fa_link="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz"
    ref_gtf_link="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz"
    ref_genome_shorthand="GRCm38.p6"
elif [[ $ref_species == "hg38" ]]; then
    output_dir="$parent_dir/hg38"
    ref_fa_link="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz"
    ref_gtf_link="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"
    ref_genome_shorthand="GRCh38.p13"
fi

############################################################################################
# SETUP
############################################################################################
# set working dirs
log_dir=$output_dir/log
ref_dir=$output_dir/ref
index_dir=$output_dir/index

# make dirs
dir_list=($ref_dir $log_dir $index_dir)
for new_dir in "${dir_list[@]}"; do if [[ ! -d $new_dir ]]; then mkdir -p $new_dir; fi; done

# create shorthand output file names
# example: $ref_dir/GRCm38.p6.fa
ref_bk=$ref_dir/BK000964.3_TPA_exp.fa

ref_gen=`echo $ref_fa_link | cut -f9 -d"/" | sed "s/.gz//g"`
ref_gen="${ref_dir}/${ref_gen}"
    
ref_gtf=`echo $ref_gtf_link | cut -f9 -d"/"| sed "s/.gz//g"`
ref_gtf="${ref_dir}/${ref_gtf}"
  
############################################################################################
# FUNCTIONS
############################################################################################
function star_index_mm10()
{
    local star_command="STAR \
        --runThreadN 32 \
        --runMode genomeGenerate \
        --genomeDir $index_dir \
        --genomeFastaFiles $ref_gen $ref_bk"
    echo "$star_command"
}

function star_index_hg38()
{
    local star_command="module load STAR; STAR \
        --runThreadN 32 \
        --runMode genomeGenerate \
        --genomeDir $index_dir \
        --genomeFastaFiles $ref_gen"
    echo "$star_command"
}

############################################################################################
# RUN
############################################################################################
#analysis
if [[ "$flag_download" == "Y" ]]; then
    echo "*** Running download ***"

    # get FA file
    # if .fa file is missing, but .fa.gz is present, gunzip it
    # if .fa.gz is missing, download it
    if [[ ! -f $ref_gen ]] && [[ -f ${ref_gen}.gz ]]; then
        echo "unzipping $ref_gen"
        gunzip ${ref_gen}.gz
    elif [[ ! -f $ref_gen ]]; then
        wget -P $ref_dir $ref_fa_link
    fi

    # additional fa not included in mm10
    if [[ $ref_species == "mm10" ]] && [[ ! -f $ref_bk ]]; then
        cp $perm_bk_path $ref_bk
    fi

    ## get GTF file
    # if .fa file is missing, but .fa.gz is present, gunzip it
    # if .fa.gz is missing, download it
    if [[ ! -f $ref_gtf ]] && [[ -f "${ref_gtf}.gz" ]]; then
        echo "unzipping $ref_gtf"
        gunzip ${ref_gtf}.gz
    elif [[ ! -f $ref_gtf ]]; then
        wget -P $ref_dir $ref_gtf_link
    fi
fi 

if [[ "$flag_index" == "Y" ]]; then
    echo "*** Indexing ***"
    module load STAR

    #create index
    if [[ $ref_species == "mm10" ]]; then
        cmd=`star_index_mm10`
    elif [[ $ref_species == "hg38" ]]; then
        cmd=`star_index_hg38`
    fi

    echo $cmd > $log_dir/index.sh
    swarm -f $log_dir/index.sh --job-name index_files --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00
fi