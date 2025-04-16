#!/bin/bash
##usage
#sh GetSNPVariants2BedSort.sh /home/zhouy0e/10K-RGP/venn_plots/
## thimmamp@kw61043:~/MAGIC16/scripts$ sh GetSNPsFromPvarsBimsForCMPlot.sh /home/zhouy0e/10K-RGP/venn_plots/
# /home/zhouy0e/10K-RGP/venn_plots$ awk -F'\t' '{print $2"\t"$3-1"\t"$3}' genome1.10K.biallelic.base_forCMplot.txt | sort -k1,1 -k2,2n > outfile
#for fname in `ls $1/*.pvar`
outfolder=/home/thimmamp/10K/data
for fname in `ls $1/genome33.*_forCMplot.txt`
do
	echo $fname
	f=`basename ${fname}`
	#outprefix=${f%%.genomewide.SNPs.withID.PlinkFormat.pvar}
	outprefix=${f%%_forCMplot.txt}
	echo ${outprefix}_variants2sorted.bed
	#awk -F'\t' '{print $2"\t"$3-1"\t"$3}' $fname | grep -v ChrUN | sort -k1,1 -k2,2n > ${outfolder}/${outprefix}_variants2sorted.bed
	bedtools intersect -a ~/MAGIC16/annotation/genome33_chromosomes_1Mb_Windows.bed -b ${outfolder}/${outprefix}_variants2sorted.bed -c > ${outfolder}/${outprefix}"_SNPs_count_1Mbwindow.txt"
done	
