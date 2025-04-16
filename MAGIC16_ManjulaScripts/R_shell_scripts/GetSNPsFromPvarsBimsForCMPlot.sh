#!/bin/bash
##usage
## thimmamp@kw61043:~/MAGIC16/scripts$ sh GetSNPsFromPvarsBimsForCMPlot.sh /home/zhouy0e/10K-RGP/venn_plots/

#for fname in `ls $1/*.pvar`
for fname in `ls $1/*.bim`
do
	echo $fname
	f=`basename ${fname}`
	#outprefix=${f%%.genomewide.SNPs.withID.PlinkFormat.pvar}
	outprefix=${f%%.bim}
	echo ${outprefix}_forCMPlot
	## for pvar files
	#grep -v "^#" $fname | cut -f3 > temp1
	#grep -v "^#" $fname | cut -f3| cut -d'_' -f1 > temp2
	#grep -v "^#" $fname | cut -f3| cut -d'_' -f2 > temp3
	#paste temp1 temp2 temp3 >${1}/${outprefix}_forCMplot.txt 
	
	## for bim files
	cut -f2 $fname > temp1
	cut -f4 $fname > temp3
	cut -f2 $fname | cut -d'_' -f1 > temp2
	paste temp1 temp2 temp3 >${1}/${outprefix}_forCMplot.txt 
	rm temp1 temp2 temp3
done	
