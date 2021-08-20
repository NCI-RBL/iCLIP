#Create error report
echo
echo "Creating error log"

#create log file, intermed error files
echo `date +"%Y%m%d_%H%M"` > error_log.txt
touch disc.txt
touch miss.txt
touch other.txt

#cehck if errros are in snakmemake config or in the rules
error_snakefile="$(grep -i "error" *_iCLIP.out | wc -l)"
error_rules="$(grep -i "error in" *_iCLIP.out | wc -l)"

#if there are no snakemake errors or rul errors
if [ $error_snakefile == "0" ] && [ $error_rules == "0" ]; then
  echo "Pipeline has successfully completed" >> error_log.txt
  echo "Pipeline has successfully completed"
#if there were only snakemake errors
elif [ $error_snakefile != "0" ] && [ $error_rules == "0" ]; then
    echo "- Errors were found in snakemake"
    echo "The following error(s) were found in the snakemake file." >> error_log.txt
    echo "*********************************************" >> error_log.txt
    grep -i "error" *_iCLIP.out | sort --unique >> error_log.txt
    echo >> error_log.txt
    echo "view error log to determine which rules failed by running: cat error_log.txt"
    echo
#if there are rule errors
else
    #cat the name of the rules with problems
    echo "- Errors were found in 1 or more snakemake rules"
    echo "The following error(s) were found in rules:" >> error_log.txt
    echo "*********************************************" >> error_log.txt
    grep -i "error in" *_iCLIP.out | sort --unique >> error_log.txt
    echo >> error_log.txt

    echo "--- Processing individual rule errors..."
    max=($(grep -oP "Error in rule \K\w+" *_iCLIP.out | sort --unique | wc -l))
    count=1

    #for each of the rule error_rules, determine if it's a memory error, file missing error, or other error
    grep -oP "Error in rule \K\w+" *_iCLIP.out | sort --unique | \
    while read -r result; do
        if grep -q "Disk quota exceeded" *$result*; then
            grep " Disk quota exceeded" *$result* | while read -r result2; do echo $result2 >> disc.txt; done
        elif grep -q "Failed to open file" *$result*; then
            grep "Failed to open file" *$result* | while read -r result2; do echo $result2 >> miss.txt; done
        else
            grep "Error" *$result* | while read -r result2; do echo $result2 >> other.txt; done
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

            echo "The following samples are affected by other error_rules and should be reviewed:" >> error_log.txt
            echo "*********************************************" >> error_log.txt
            cat other.txt >> error_log.txt
            echo >> error_log.txt
            rm other.txt
        fi

    done

    echo "- Processing complete. View error log to determine which rules failed by running: cat error_log.txt"
    echo
fi