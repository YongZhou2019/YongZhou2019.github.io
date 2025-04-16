#!/bin/bash
for fname in `ls /home/thimmamp/MAGIC16/indels_results/*transcripts*.bed`
do
	f=`basename ${fname%%.bed}`
        #temp=${f%%.bed}
        IFS='_' read -ra g <<< "$f"
	annot="/home/thimmamp/MAGIC16/annotation/${g[4]}/sequences.fa"
	#echo $fname $annot
	echo $fname ${g[4]} $annot
	#/home/thimmamp/bioinformatics_tools/bedtools2-2.31.1/bin/bedtools getfasta -fi $annot -bed $fname -nameOnly -fo /home/thimmamp/MAGIC16/indels_results/${f}_transcripts_extractedSeqs.fa
	#echo $cmd
	#`${cmd}`
	## for pvar files
	#grep -v "^#" $fname | cut -f3| cut -d'_' -f1,2 > $prefix
	## for bim files
	#cut -f2 $fname| cut -d'_' -f1,2 > $prefix 
done	
