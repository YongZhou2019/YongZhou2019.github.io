#!/bin/bash
echo "genome Numofgenes"
for g in genome1 genome15 genome5 genome14 genome7 genome28 genome33 genome13 genome16 genome8 genome11 genome10 genome9 genome27 genome4 genome12
do
        #echo $fname
	#f=`basename ${fname}`
	#IFS='.' read -ra ADDR <<< "$f"
	TOTAL=`awk '{if($3=="gene")print}' /home/thimmamp/MAGIC16/annotation/${g}/genes.gff| cut -f9 | cut -d';' -f1 | cut -d':' -f2 | grep -v Un | sort | uniq | wc -l`
	echo ${g} $TOTAL
done


### awk '{if($3=="gene")print}' ../annotation/genome1/genes.gff| cut -f9 | cut -d';' -f1 | cut -d':' -f2 | grep -v Un | sort | uniq | wc -l
#37007
