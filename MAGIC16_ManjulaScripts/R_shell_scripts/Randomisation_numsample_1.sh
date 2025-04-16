#!/bin/bash

## Usage : ./Randomisation_60.sh ../dataset/SNPs/20KSamples.txt <percentage/numofsample>
numofsample=${2}
fname="Randomisation_numsample_${numofsample}_biallelic.txt"

#for percentage in {10..100..10}; 
#do 
	#for i in {1..100}; do
#for i in {1..10}; do
for i in {1..2}; do
	#sample_size=$(($(wc -l < $1) * percentage / 100))
	#sample_size=${numofsample}
        #shuf -n "$sample_size" $1 > "sample_${numofsample}_${i}.txt"
        shuf -n "$2" $1 > "sample_${2}_${i}.txt"
	#plink2 --pfile /home/projects/10KRGP/dataset/SNPs/genome1.genomewide.SNPs.withID.PlinkFormat --keep "sample_"${2}"_"${i}".txt" --make-pgen pvar-cols=xheader,qual,filter,info --sort-vars --allow-extra-chr --min-ac 1 --out "rs_"${2}"_"${i}"_plink"
	plink2 --pfile /home/projects/10KRGP/dataset/SNPs/genome1.biallelic.base.missing0.2.maf0.01.genomewide.SNPs.withID.PlinkFormat --keep "sample_"${2}"_"${i}".txt" --make-pgen pvar-cols=xheader,qual,filter,info --sort-vars --allow-extra-chr --min-ac 1 --out "rs_"${2}"_"${i}"_plink"
	newsnpfname="Randomisation_numsample_${2}_${i}_biallelic_snps.txt"
	grep -v "^#" "rs_${2}_${i}_plink.pvar" | grep -v "ChrUN" | cut -f3 > $newsnpfname 
	snpcount=$(grep -v "^#" "rs_${2}_${i}_plink.pvar" | grep -v "ChrUN" | cut -f3 | wc -l)
	#echo $2 $i $snpcount >> $fname
	rm "rs_${2}_${i}_plink.pvar"
	rm "rs_${2}_${i}_plink.pgen"
	rm "rs_${2}_${i}_plink.psam"
	rm "rs_${2}_${i}_plink.log"
	rm "sample_${2}_${i}.txt"
	echo $2 $i $snpcount >> $fname
#done
done

