#!/bin/bash
echo "genome insertions deletions both total"
#for fname in `ls ./genome*_*/*_annotated_col1_9.vcf`
for fname in `ls ./*INDELs.withID.PlinkFormat.vcf`
do
        #echo $fname
	f=`basename ${fname}`
	IFS='.' read -ra ADDR <<< "$f"
	TOTAL=`awk -F'\t' '{if($3~/Chr/)print}' $fname| cut -f1-9 | wc -l`
	#INS=`awk -F'\t' '{if($8~/ins/&&$8!~/del/)print $3,$8}' $fname | wc -l`
        #DELS=`awk -F'\t' '{if($8!~/ins/&&$8~/del/)print $3,$8}' $fname | wc -l`
        #BOTH=`awk -F'\t' '{if($8~/ins/&&$8~/del/)print $3,$8}' $fname | wc -l`
	#TOTAL=$((INS+DELS+BOTH))
	#echo ${ADDR[0]} $INS $DELS $BOTH $TOTAL
	echo ${ADDR[0]} $TOTAL
done

