#!/bin/bash
#Combination of sort16greengenes.sh and muscletree.sh

OUTPUT_DIR=$1


#sort16greengenes.sh
if [ -d "$OUTPUT_DIR" ]
 then
   echo "Directory exists"
 else
   echo "Creating directory"
   mkdir $OUTPUT_DIR
   mkdir $OUTPUT_DIR/greengenes/
   mkdir $OUTPUT_DIR/muscle/
   mkdir $OUTPUT_DIR/iqtree/
fi

grep "Methanofollis liminatans" /Users/alexandranolan/Desktop/16scandidate/greengenes/greenidacc > $OUTPUT_DIR/greengenes/green_acid

while read id acc
  do echo $acc
   cat /Users/alexandranolan/Desktop/FYP/datam/greengenes/gg_13_5.fasta | grep -A 1 ">$id$" >> $OUTPUT_DIR/greengenes/green_acid.fasta
   done < $OUTPUT_DIR/greengenes/green_acid

#muscletree.sh

muscle -in $OUTPUT_DIR/greengenes/green_acid.fasta -out $OUTPUT_DIR/muscle/muscle.fasta

iqtree -s $OUTPUT_DIR/muscle/muscle.fasta -nt AUTO  
 
open -a "FigTree v1.4.3.app" $OUTPUT_DIR/muscle/muscle.fasta.treefile