#!/bin/bash
## MOREF=`sudo run command against "$VAR1" | grep name | cut -c7-`
##echo ${1}
#file="/home/thimmamp/MAGIC16/tmp/vcf/genome6_sample_names_from_header.txt"
total=`wc -l $1`
echo "${total}"
num=$((total*0.1))
echo "${num}"

#echo $total
for i in {1..2} ## Number of random runs planned for 100
do
	echo "rand count ${i}"
	for percent in 0.1 0.2  ## percentage of samples from 3024
	do
		echo "${total} ${percent}"
		num=$((total * percent)) ## get number of percentage of samples
		echo "${num}"
#		#`shuf -n "$num" "$1"` > subsample ## get randomly shuffled percentage of samples
#		#`cat "$subsample"
		numdels=bcftools view -S $subsample genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type  | grep -oP "TYPE=\w+" | sort | uniq -c | grep del | cut -d' ' -f1
#1558785
		numins=bcftools view -S $subsample genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type  | grep -oP "TYPE=\w+" | sort | uniq -c | grep ins | cut -d' ' -f1
		echo "${percent}\t${i}\t${num}\t${numins}\t${numdels}\n"
	done
done	
