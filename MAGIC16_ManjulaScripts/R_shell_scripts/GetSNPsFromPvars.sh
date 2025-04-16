#!/bin/bash
for fname in `ls /home/zhouy0e/10K-RGP/venn_plots/*.pvar`
#for fname in `ls /home/zhouy0e/10K-RGP/venn_plots/*prune2.bim`
outprefix=${f%%.gfa}
do
	echo $fname
	f=`basename ${fname}`
        IFS='.' read -ra ADDR <<< "$f"
        #echo ${ADDR[0]}
        #echo ${ADDR[1]}
        #echo ${ADDR[2]}
        #echo ${ADDR[3]}
        prefix=${ADDR[0]}.${ADDR[1]}.${ADDR[2]}.${ADDR[3]}.${ADDR[4]}.prune2
	echo $prefix
	## for pvar files
	grep -v "^#" $fname | cut -f3| cut -d'_' -f1,2 > $prefix
	## for bim files
	#cut -f2 $fname| cut -d'_' -f1,2 > $prefix 
done	
