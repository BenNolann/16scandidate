#!/bin/bash 
# Given an organism genus and species name, take out the 16S data from greengenes.
# Requires $1 to be a a file only containing the id and accessions of the organism of interest
 


OUTPUT_DIR=$1

if [ -d "$OUTPUT_DIR" ] 
then
    echo "Directory exists." 
else
    echo "Creating directory"
    mkdir $OUTPUT_DIR
fi



# Input organism name after grep
grep "Methanofollis liminatans" /Users/alexandranolan/Desktop/16scandidate/greengenes/greenidacc > $OUTPUT_DIR/green_acid


while read id acc
  do echo $acc
 cat /Users/alexandranolan/Desktop/FYP/datam/greengenes/gg_13_5.fasta | grep -A 1 ">$id$" >> $OUTPUT_DIR/green_acid.fasta 
  done < $OUTPUT_DIR/green_acid