#!/bin/bash

INPUT_FILE=$1
OUTPUT_DIR=$2

#Download reads for a specific organism from NCBI based on results of ggrrndbsra.R
#INPUT_FILE containing accessions obtained from ggrrndbsra.R result for chosen organism.
while read acc
  do
    fastq-dump --split-files $acc
  echo $acc
  done < $INPUT_FILE 



#Use pyani to find the best reference genome for each genome

for file in genome_files
 do 
   run_ANI.sh -e EXPERIMENTNAME -f FORWARD -r REVERSE  -n 10 -o $OUTPUT_DIR
  echo $file
done

#Use riboseed to assemble the reads 

for file in riboseed_files 
 do
 ribo run -r  -F  -R  --cores 16 --threads 1 -v 1 --serialize -o 
  echo $file
done