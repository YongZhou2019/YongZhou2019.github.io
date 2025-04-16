#!/bin/bash

## Usage : ./Randomisation_90.sh ../dataset/SNPs/20KSamples.txt
percentage=90
fname="Randomisation_${percentage}.txt"
#for percentage in {10..100..10}; 
#do 
	#for i in {1..100}; do
for i in {1..100}; do
	sample_size=$(($(wc -l < $1) * percentage / 100))
        shuf -n "$sample_size" $1 > "sample_${percentage}_${i}.txt"
	plink2 --pfile /home/projects/10KRGP/dataset/SNPs/genome1.genomewide.SNPs.withID.PlinkFormat --keep "sample_"${percentage}"_"${i}".txt" --make-pgen pvar-cols=xheader,qual,filter,info --sort-vars --allow-extra-chr --min-ac 1 --out "rs_"${percentage}"_"${i}"_plink"
	snpcount=$(grep -v "^#" "rs_${percentage}_${i}_plink.pvar" | grep -v "ChrUN" | cut -f3 | wc -l)
	echo $percentage $i $snpcount >> $fname
	rm "rs_${percentage}_${i}_plink.pvar"
	rm "rs_${percentage}_${i}_plink.pgen"
	rm "rs_${percentage}_${i}_plink.psam"
	rm "rs_${percentage}_${i}_plink.log"
	rm "sample_${percentage}_${i}.txt"
#done
done

