#!/bin/bash


for g in genome1 genome4 genome5 genome7 genome8 genome9 genome10 genome11 genome12 genome13 genome14 genome15 genome16 genome33 genome27 genome28
do
	python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/annotation/${g}/genes.gff /home/thimmamp/MAGIC16/tmp/vcf/${g}_processed/${g}.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/gramene_genome33_Only_YongUpdatedFromZhiChao/specific/${g}.*.gset.specific > ~/MAGIC16/indels_results/${g}_specific_indelscount_test.txt
	#python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/annotation/${g}/genes.gff /home/thimmamp/MAGIC16/tmp/vcf/${g}_processed/${g}.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/gramene_genome33_Only_YongUpdatedFromZhiChao/specific/${g}.*.gset.specific
	#python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/annotation/genome6_genes.gff /home/thimmamp/MAGIC16/tmp/vcf/genome6_processed/genome6.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/magic16_homologs_20230125/specific/genome6.IR64.gset.specific
	#echo $cmd
done
