#!/bin/bash
for fname in `ls ~/MAGIC16/annotation/genome*/genes.gff`; 
do 
	echo $fname
	awk -F'\t' '{if($3=="mRNA")print}' $fname | head 
done
