# set args
qc_dir="$1"
single_qc_threshold="$2"
project_qc_threshold="$3"

#qc_dir="/data/RBL_NCI/Wolin/8-09-21-HaCaT_fCLIP_v2.0_2/log/STAR"
qc_failure=$qc_dir/qc_failure.txt
qc_pass=$qc_dir/qc_pass.txt

#single_qc_threshold=95
#project_qc_threshold=50

# set counters
total_unmapped=0
counter=0

# create logs
if [[ -f $qc_failure ]]; then rm $qc_failure; fi; touch $qc_failure
if [[ -f $qc_pass ]]; then rm $qc_pass;  fi; touch $qc_pass

# for each of the STAR log files, determine unmapped read percentages
for f in $qc_dir/*; do
    mismatch=`cat $f | grep "% of reads unmapped: too many mismatches" | cut -f2 -d"|" | cut -f1 -d"." | sed "s/%//g"`
    tooshort=`cat $f | grep "% of reads unmapped: too short" | cut -f2 -d"|" | cut -f1 -d"." | sed "s/%//g"`
    other=`cat $f | grep "% of reads unmapped: other" | cut -f2 -d"|" | cut -f1 -d"." | sed "s/%//g"`

    #check unmapped in any one samples is under treshold
    single_unmapped=$((mismatch + tooshort + other))
    if [[ $single_unmapped -gt $single_qc_threshold ]]; then 
        echo "Sample $f does not meet the requirements ($single_qc_threshold%), with a percent unampped reads of $single_unmapped%. Please review details and if acceptable \
            increase the single_qc_threshold accordingly. If not acceptable, remove sample from sample_manifest and multiplex_manifest to continue" >> $qc_failure
    fi

    # save single to total unmapped
    total_unmapped=$((total_unmapped + single_unmapped))
    counter=$((counter+1))
done

# check total unmapped across all samples is under treshold
average_unmapped=$((total_unmapped / counter))
if [[ $average_unmapped -gt $project_qc_threshold ]]; then 
    echo "Project unmapped percentages do not meet the requirements ($project_qc_threshold%), with an averaged percent unampped reads of $average_unmapped%. Please review details and if acceptable \
    increase the project_qc_threshold accordingly. If not acceptable, remove sample{s} from sample_manifest and multiplex_manifest to continue". >> $qc_failure
fi

# # check if project has any failures, if so create failure log and fail project
length_fails=`cat $qc_failure | wc -l `
if [[ $length_fails -eq 0 ]]; then
    rm $qc_failure
else
    rm $qc_pass
    echo "fail qc"
fi