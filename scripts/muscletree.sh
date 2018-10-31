#!/bin/bash

INPUT_FILE=$1
OUTPUT_DIR=$2

#Usage: bash muscletree.sh *.fasta .fasta

#muscle requires a series of fasta file to align, output file must contain .fasta
muscle -in $INPUT_FILE -out $OUTPUT_DIR

#iqtree requires a multiple sequence aligned fasta file
iqtree -s $OUTPUT_DIR -nt AUTO 

#figtree requires a tree
open -a "FigTree v1.4.3.app" $OUTPUT_DIR.treefile 