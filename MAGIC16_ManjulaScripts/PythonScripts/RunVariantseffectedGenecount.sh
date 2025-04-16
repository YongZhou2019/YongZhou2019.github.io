#!/bin/bash
#echo "genome insertions deletions both total"
#for fname in `ls ./genome*_*/*_annotated_col1_9.vcf`
for fname in `ls /home/thimmamp/MAGIC16/tmp/vcf/genome28_processed/genome*variants_effect*.txt`
do
        #echo $fname
	#f=`basename ${fname}`
	#IFS='_' read -ra ADDR <<< "$fname"
	#numgenes=`cut -f2 $fname| cut -d'_' -f2| sort | uniq | wc -l`
	#numgenes=`cut -f2 $fname| cut -d'_' -f2 | sort | uniq |wc -l` # genome1
	#numgenes=`cut -f2 $fname| cut -d'_' -f2 | sort | uniq | wc -l` #genome4
	#numgenes=`cut -f2 $fname | cut -d'_' -f2,3 | cut -d'.' -f1 | sort | uniq | wc -l` #genome5 - 16
	#numgenes=`cut -f2 $fname| cut -d'.' -f1| sort | uniq | wc -l`
	numgenes=`cut -f2 $fname| sort | uniq | wc -l` #genome27 and 28
	#INS=`awk -F'\t' '{if($8~/ins/&&$8!~/del/)print $3,$8}' $fname | wc -l`
        #DELS=`awk -F'\t' '{if($8!~/ins/&&$8~/del/)print $3,$8}' $fname | wc -l`
        #BOTH=`awk -F'\t' '{if($8~/ins/&&$8~/del/)print $3,$8}' $fname | wc -l`
	#TOTAL=$((INS+DELS+BOTH))
	#echo ${ADDR[0]} $INS $DELS $BOTH $TOTAL
	echo $fname $numgenes
done

