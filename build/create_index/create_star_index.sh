############################################################################################
# USER INPUTS
############################################################################################
output_dir="/data/CCBR_Pipeliner/iCLIP/index/active/2022_0505/mm10"

ref_species="mm10"
ref_fa_link="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz"
ref_gtf_link="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz"
ref_genome_shorthand="GRCm38.p6"

############################################################################################
# SETUP
############################################################################################
# create shorthand output file names
ref_bk=$ref_dir/BK000964.3_TPA_exp.fa

# example: $ref_dir/GRCm38.p6.fa
ref_gen=`echo $ref_gtf_link | cut -f8 -d"/" | sed "s/.genome.fa.gz//"`
ref_gen="${ref_dir}/${ref_gen}.fa"
    
ref_gtf=`echo $ref_gtf_link | cut -f8 -d"/" | sed "s/.annotation.gtf.gz//"`
ref_gtf="${ref_gtf}/${ref_gen}.gtf"

# set dirs
perm_bk_path="/data/CCBR_Pipeliner/iCLIP/index/active/2021_0607/mm10/01_source/BK000964.3_TPA_exp.fa"
log_dir=$output_dir/log
ref_dir=$output_dir/ref
index_dir=$output_dir/index

# make dirs
dir_list=($ref_dir $log_dir $index_dir)
for new_dir in "${dir_list[@]}"; do if [[ ! -d $new_dir ]]; then mkdir -p $new_dir; fi; done
  
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

function star_index_human()
{
    local star_command="STAR \
        --runThreadN 32 \
        --runMode genomeGenerate \
        --genomeDir $index_dir \
        --genomeFastaFiles $ref_gen"
    echo "$star_command"
}

############################################################################################
# Controls
############################################################################################
flag_download="Y"
flag_index="N"

############################################################################################
# RUN
############################################################################################
#analysis
if [[ "$flag_download" == "Y" ]]; then
    echo "*** Running download ***"

    # get FA file
    if [[ ! -f $ref_gen ]]; then
        wget -P $ref_dir $ref_fa_link
        gunzip ${ref_gen}.gz
    fi

    # additional fa not included in mm10
    if [[ $species == "mm10" ]] & [[ ! -f $ref_bk ]]; then
        cp $perm_bk_path $ref_bk
    fi

    ## get GTF file
    if [[ ! -f $ref_gtf ]]; then
        wget -P $ref_dir $ref_gtf_link
        gunzip ${ref_gtf}.gz
    fi
fi 

if [[ "$flag_index" == "Y" ]]; then
    echo "*** Indexing ***"
    module load STAR

    #create index
    if [[ $species == "mm10" ]]; then
        cmd=`star_index_mm10`
    elif [[ $species == "human" ]]; then
        cmd=`star_index_human`
    fi

    echo $cmd > $log_dir/index.sh
    swarm -f $log_dir/index.sh --job-name index_files --merge-output --logdir $log_dir --partition=norm -g 40 -t 4 --time 10:00:00

fi