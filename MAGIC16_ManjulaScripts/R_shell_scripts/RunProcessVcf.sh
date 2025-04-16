#!/bin/bash
#awk '{if($5!~/\*/)print}' genome12.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome12.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome7.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome7.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome15.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome15.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome5.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome5.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome14.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome14.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome28.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome28.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome13.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome13.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome16.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome16.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome8.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome8.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome9.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome9.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome10.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome10.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome11.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome11.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome27.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome27.nostar.bimultiallele.50bpindel_annotated.vcf
#awk '{if($5!~/\*/)print}' genome33.filterCount5.genomewide.INDELs.withID.PlinkFormat_conversion.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type | grep -v "^##" > genome33_processed/genome33.nostar.bimultiallele.50bpindel_annotated.vcf

## Get first nine columns of vcf to input snpEff
#cut -f1-9 genome10.nostar.bimultiallele.50bpindel_annotated.vcf > genome10.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome11.nostar.bimultiallele.50bpindel_annotated.vcf > genome11.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome12.nostar.bimultiallele.50bpindel_annotated.vcf > genome12.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome13.nostar.bimultiallele.50bpindel_annotated.vcf > genome13.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome14.nostar.bimultiallele.50bpindel_annotated.vcf > genome14.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome15.nostar.bimultiallele.50bpindel_annotated.vcf > genome15.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome16.nostar.bimultiallele.50bpindel_annotated.vcf > genome16.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome27.nostar.bimultiallele.50bpindel_annotated.vcf > genome27.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome28.nostar.bimultiallele.50bpindel_annotated.vcf > genome28.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome4.nostar.bimultiallele.50bpindel_annotated.vcf > genome4.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome5.nostar.bimultiallele.50bpindel_annotated.vcf > genome5.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome7.nostar.bimultiallele.50bpindel_annotated.vcf > genome7.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome8.nostar.bimultiallele.50bpindel_annotated.vcf > genome8.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome9.nostar.bimultiallele.50bpindel_annotated.vcf > genome9.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf
#cut -f1-9 genome33.nostar.bimultiallele.50bpindel_annotated.vcf > genome33.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf

###

#for fname in `ls ./genome*/*_annotated_col1_9.vcf`
for fname in `ls ./genome33_processed/*_annotated_col1_9.vcf`
do
        echo $fname
        #java -jar  /home/zhouy0e/tools/SnpEff/snpEff/snpEff.jar ann genome6 genome6.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_col1_9_all.vcf2 -s genome6.nostar.bimultiallele.50bpindel_annotated_col1_9_all_snpeff_summary_all.html2 -csvStats genome6.nostar.bimultiallele.50bpindel_annotated_col1_9_all.csv2
        IFS='.' read -ra ADDR <<< "$fname"
        #echo ${ADDR[0]}
        #echo ${ADDR[1]}
        #echo ${ADDR[2]}
        #echo ${ADDR[3]}
        prefix=${ADDR[0]}.${ADDR[1]}.${ADDR[2]}.${ADDR[3]}
	awk '{if($8~/ins/&&$8!~/del/)print}' $fname > ${prefix}_insertions.vcf #only ins
	awk '{if($8!~/ins/&&$8~/del/)print}' $fname > ${prefix}_deletions.vcf  #only dels
        awk '{if($8~/ins/&&$8~/del/)print}' $fname > ${prefix}_insertions_deletions.vcf	 #both ins and dels
        #java -jar /home/zhouy0e/tools/SnpEff/snpEff/snpEff.jar ann ${ADDR[0]} $fname > ${prefix}_all.vcf2 -s ${prefix}_all_snpeff_summary.html2 -csvStats ${prefix}_all.csv2
done

