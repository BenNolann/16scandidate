# Given an organism genus and species name, take out the 16S data from greengenes.

INPUT= $1
OUTPUT= $2


 while read id acc; do echo $acc;
 cat ./Desktop/FYP/datam/greengenes/gg_13_5.fasta | grep -A 1 ">$id" >> $2; 
  done < $1