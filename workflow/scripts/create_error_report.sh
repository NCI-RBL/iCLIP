#Create error report
echo
echo -e "****Creating error log****\n"

#create log file, intermed error files
echo `date +"%Y%m%d_%H%M"` > error_log.txt
touch disc.txt
touch miss.txt
touch other.txt

#cehck if errros are in snakmemake config or in the rules
#do not include errors in rules: error in, error executing, error message
#do not include error parameters: printableError
error_snakefile="$(cat *iCLIP_* | grep -v "output_file_error" | grep -v "Error in" | \
grep -v "Error executing" | grep -v "error message" | grep -v "printableError" | grep -i "error\|Error" | wc -l)"

#do not include output_file_error
#do not include error parameters: printableError
error_rules="$(cat *iCLIP_* | grep -v "output_file_error" | grep -v "printableError" | grep -i "error in" | wc -l)"

#if there are no snakemake errors or rul errors
if [ $error_snakefile == "0" ] && [ $error_rules == "0" ]; then
  echo "-Pipeline has successfully completed" >> error_log.txt
  echo -e "-Pipeline has successfully completed\n"
  
#if there were only snakemake errors
elif [ $error_snakefile != "0" ] && [ $error_rules == "0" ]; then
    echo "- Errors were found in snakemake"
    echo "The following error(s) were found in the snakemake file." >> error_log.txt
    echo "*********************************************" >> error_log.txt
    cat *iCLIP_* | grep -v "output_file_error" | grep -v "Error in" | grep -v "Error executing" | grep -v "error message" | \
    grep -v "printableError" | grep -i "error\|Error" | sort --unique >> error_log.txt
    echo >> error_log.txt
    echo "view error log to determine which rules failed by running: cat error_log.txt"
    echo
#if there are rule errors
else
    #cat the name of the rules with problems
    echo "- Errors were found in 1 or more snakemake rules"
    echo "The following error(s) were found in rules:" >> error_log.txt
    echo "*********************************************" >> error_log.txt
    grep -i "error in" *iCLIP_* | sort --unique >> error_log.txt
    echo >> error_log.txt

    echo "--- Processing individual rule errors..."
    max=($(grep -oP "Error in rule \K\w+" *iCLIP_* | sort --unique | wc -l))
    count=1

    #for each of the rule error_rules, determine if it's a memory error, file missing error, or other error
    grep -oP "Error in rule \K\w+" *iCLIP_* | sort --unique | \
    while read -r result; do
        #if it's a memory issue
        if grep -q "Disk quota exceeded" *$result*; then
            grep " Disk quota exceeded" *$result* | while read -r result2; do echo $result2 >> disc.txt; done
        #if it's a permissions issue
        elif grep -q "Failed to open file" *$result*; then
            grep "Failed to open file" *$result* | while read -r result2; do echo $result2 >> miss.txt; done
        #if there is something else
        else
            grep "Error" *$result* | while read -r result2; do echo $result2 >> other.txt; done
            grep "error" *$result* | while read -r result2; do echo $result2 >> other.txt; done
        fi

        count=$(($count + 1))

        #once complete with all rules, cat the output to the error log and remove intermed files
        if (( count > max )); then
            echo "The following samples are affected by memory and must be deleted:" >> error_log.txt
            echo "*********************************************" >> error_log.txt
            cat disc.txt >> error_log.txt
            echo >> error_log.txt
            rm disc.txt

            echo "The following samples are affected by missing input files/output dir and should be reviewed:" >> error_log.txt
            echo "*********************************************" >> error_log.txt
            cat miss.txt >> error_log.txt
            echo >> error_log.txt
            rm miss.txt

            echo "The following samples are affected by rule failures and should be reviewed:" >> error_log.txt
            echo "--listed by snakemake rule, sample ID, alignment type, split file, any related slurm ids" >> error_log.txt
            echo "*********************************************" >> error_log.txt
            rule_list=(`cat other.txt | awk -F"." '{print$1}' | sort --unique`); \
            for r in ${rule_list[@]}; do \
                echo "*******$r*****" >> error_log.txt
                #for each sample
                sample_list=(`cat other.txt | grep -o "$r[.][0-9]*[.]al=[a-z]*[,]n=[0-9]*[,]sp=.*err:" | awk -F"," '{print$3}' | grep -o -P '(?<=sp[=]).*(?=.err)' | sort --unique`); \
                for s in ${sample_list[@]}; do \
                    echo $s >> error_log.txt
                        
                    #for each alignment file type
                    al_list=(`cat other.txt | grep -o "[0-9]*[.]al=[a-z]*[,]n=[0-9]*[,]sp=$s.err:" | awk -F"." '{print$2}' | grep -o -P '(?<=al[=]).*?(?=[,])' | sort --unique`); \
                    for al in ${al_list[@]}; do \
                        echo "--$al" >> error_log.txt
                            
                        #for each split sample
                        n_list=(`cat other.txt | grep -o "[0-9]*[.]al=$al[,]n=[0-9]*[,]sp=$s.err:" | awk -F"." '{print$2}' | grep -o -P '(?<=n[=]).*?(?=[,])' | sort --unique`); \
                        for n in ${n_list[@]}; do \
                                
                            #for each slurm id
                            s_list=(`cat other.txt | grep -o "[0-9]*[.]al=$al[,]n=$n*[,]sp=$s.err:" | awk -F"." '{print$1}' | sort --unique`); \
                            slurm_ex1="cat *$r*${s_list[0]}*$al*$n*$s*.err"
                            slurm_ex2="cat *$r*${s_list[1]}*$al*$n*$s*.err"
                            echo "----$n" >> error_log.txt
                            echo "------$slurm_ex1" >> error_log.txt
                            echo "------$slurm_ex2" >> error_log.txt
                        done
                    done
                done
            done
        fi
    done
    echo "- Processing complete. View error log to determine which rules failed by running: cat error_log.txt"
    echo
fi
