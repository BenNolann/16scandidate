#!/bin/bash
set -e 
OUTPUT_DIR=$1
ORGANISM=$2

#Download reads for a specific organism from NCBI based on results of .R
#INPUT_FILE containing accessions obtained from ggrrndbsra.R result for chosen organism

if [ -z "$OUTPUT_DIR" ]
  then
    echo "Error: Directory exists!"
    exit 1
  else
    echo "Creating directory"
  mkdir $OUTPUT_DIR
fi


grep "${ORGANISM}" ./srapure > ${OUTPUT_DIR}/srafind 

# I have no idea why this wont work in a while loop but does runniong a list of cmds

while read sra gs 
  do 
    echo "fastq-dump --split-files $sra -O ${OUTPUT_DIR}" >> ${OUTPUT_DIR}/runcmds
    #fastq-dump --split-files "$sra" -O ${OUTPUT_DIR}
done < ${OUTPUT_DIR}/srafind
bash ${OUTPUT_DIR}/runcmds

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