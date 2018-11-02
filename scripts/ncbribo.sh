#!/bin/bash

OUTPUT_DIR=$1

#Download reads for a specific organism from NCBI based on results of ggrrndbsra.R
#INPUT_FILE containing accessions obtained from ggrrndbsra.R result for chosen organism.

if [ -z "$OUTPUT_DIR" ]
  then
echo "Directory exists"
  else
echo "Creating directory"
  mkdir $OUTPUT_DIR
fi


grep "Pediococcus acidilactici" /Users/alexandranolan/Desktop/16scandidate/srafind/srapure >$OUTPUT_DIR/srafind 

while read sra gs 
  do echo $sra
    cat /Users/alexandranolan/Desktop/16scandidate/srafind/srapure | grep -A 1 ">$sra$"
    fastq-dump --split-files $sra  > $OUTPUT_DIR/sra.fastq 
done < $OUTPUT_DIR/srafind 


exit

#Use pyani to find the best reference genome for each genome

for sra in $OUTPUT_DIR/sra.fastq
  do 
    run_ANI.sh -e ORGANISMNAME -f FORWARD -r REVERSE  -n 10 -o $OUTPUT_DIR
   echo $sra
 done

#Use riboseed to assemble the reads 

#for file in riboseed_files 
 #do
 #ribo run -r  -F  -R  --cores 16 --threads 1 -v 1 --serialize -o 
  #echo $file
#done