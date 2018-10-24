#!/bin/bash

INPUT_FILE=$1
OUTPUT_DIR=$2


if [ -z "$OUTPUT_DIR" ]
then
echo "USAGE: extract_16s_from_contigs.sh input_contigs.fasta output_directory/"
exit 1
fi

# make output dir
mkdir $OUTPUT_DIR

# identify 16Ss
barrnap $INPUT_FILE | grep "16S"  > $OUTPUT_DIR/ribo16.gff

cut -f 1,4,5,7  $OUTPUT_DIR/ribo16.gff  >  $OUTPUT_DIR/ribo16.coords

while read chrom start end ori
do 
    # set RC or not
    if [ $ori == "-" ]
    then 
	SUFFIX="-RC_"  
    else 
	SUFFIX=""
    fi
    #echo \'"$chrom$SUFFIX :$start:$end"\'
    ~/open_utils/extractRegion/extractRegion.py \'"$chrom$SUFFIX :$start:$end"\' $INPUT_FILE -v 1 >> $OUTPUT_DIR/ribo16.fasta
done < $OUTPUT_DIR/ribo16.coords