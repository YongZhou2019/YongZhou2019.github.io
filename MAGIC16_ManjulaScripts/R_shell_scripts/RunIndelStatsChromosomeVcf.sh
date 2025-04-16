#!/bin/bash
echo "genome chromosome insertions deletions both total"
#for fname in `ls ./genome*_*/*_annotated_col1_9.vcf`
for fname in `ls ./genome33_*/*_annotated_col1_9.vcf`
do
        echo $fname
	f=`basename ${fname}`
	IFS='.' read -ra ADDR <<< "$f"
	for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12
	do
		INS=$(awk -F'\t' -v chr="$chr" '$3~chr"_" && $8~/ins/ && $8!~/del/ {print $3,$8}' "$fname" | wc -l)
		DELS=$(awk -F'\t' -v chr="$chr" '$3~chr"_" && $8!~/ins/ && $8~/del/ {print $3,$8}' "$fname" | wc -l)
		BOTH=$(awk -F'\t' -v chr="$chr" '$3~chr"_" && $8~/ins/ && $8~/del/ {print $3,$8}' "$fname" | wc -l)

		#INS=`awk -F'\t' '{if($3~/"$chr_"/ && $8~/ins/ && $8!~/del/)print $3,$8}' $fname | wc -l`
        	#DELS=`awk -F'\t' '{if($3~/${chr}_/ && $8!~/ins/ && $8~/del/)print $3,$8}' $fname | wc -l`
        	#BOTH=`awk -F'\t' '{if($3~/${chr}_/ && $8~/ins/ && $8~/del/)print $3,$8}' $fname | wc -l`
		TOTAL=$((INS+DELS+BOTH))
		echo ${ADDR[0]} $chr $INS $DELS $BOTH $TOTAL
	done
done

