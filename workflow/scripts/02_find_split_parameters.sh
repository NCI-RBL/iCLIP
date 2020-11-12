tsvfile=$1

#all args after file name will be processed
for file in "${@:2}"; do
  #determine the length of the file, and number of split files to be created
  #based on N of lines
  filelen=$(sed -n \$= "$file")

  filenum=$(echo $filelen |
  awk '{{print (($0/12000)-int($0/12000)>0)?int($0/12000)+1:int($0/12000)}}')

  #determine size of each split file, determine if that value is div by 4
  #since seq chunks are in sets of 4
  chunk=$(($filelen/$filenum))
  chunki=$(echo $chunk | awk '{{$0=int($0/4+1)*4}}1')

  #print file name, n of files, chunk size of each file to output
  echo -e $file "\t" $filenum "\t" $chunki >> $tsvfile;
done
