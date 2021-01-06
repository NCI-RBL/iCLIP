tsvfile=$1
n_sequences=$2

#determine N lines - each seq is 4 lines
n_lines=$(($n_sequences*4))

#all args after file name will be processed
for file in "${@:3}"; do
  #determine the length of the file, and number of split files to be created
  #based on N of lines
  file_len=$(sed -n \$= "$file")

  if [[ "$(($file_len/$n_lines))" -lt 1 ]]; then
    raw_split=1
  else
    raw_split=$((($file_len/$n_lines) + 1))
  fi

  filenum=$(echo $raw_split |
  awk '{{print (($0)-int($0)>0)?int($0)+1:int($0)}}')

  #determine size of each split file, determine if that value is div by 4
  #since seq chunks are in sets of 4
  chunk=$(($file_len/$filenum))
  chunki=$(echo $chunk | awk '{{$0=int($0/4+1)*4}}1')

  #print file name, n of files, chunk size of each file to output
  echo -e $file "\t" $filenum "\t" $chunki >> $tsvfile;
done
