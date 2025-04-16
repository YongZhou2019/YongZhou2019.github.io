#!/bin/bash
## Splitting by chromosome
#awk -F'\t' '{if($3~/Chr01_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr01.vcf
#awk -F'\t' '{if($3~/Chr02_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr02.vcf
#awk -F'\t' '{if($3~/Chr03_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr03.vcf
#awk -F'\t' '{if($3~/Chr04_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr04.vcf
#awk -F'\t' '{if($3~/Chr05_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr05.vcf
#awk -F'\t' '{if($3~/Chr06_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr06.vcf
#awk -F'\t' '{if($3~/Chr07_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr07.vcf
#awk -F'\t' '{if($3~/Chr08_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr08.vcf
#awk -F'\t' '{if($3~/Chr09_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr09.vcf
#awk -F'\t' '{if($3~/Chr10_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr10.vcf
#awk -F'\t' '{if($3~/Chr11_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr11.vcf
#awk -F'\t' '{if($3~/Chr12_/)print}' genome6.nostar.bimultiallele.50bpindel_annotated.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_Chr12.vcf

## adding header for parser
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr01.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr01.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr02.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr02.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr03.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr03.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr04.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr04.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr05.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr05.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr06.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr06.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr07.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr07.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr08.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr08.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr09.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr09.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr10.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr10.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr11.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr11.vcf
cat genome6_vcf_header.txt genome6.nostar.bimultiallele.50bpindel_annotated_Chr12.vcf > genome6.nostar.bimultiallele.50bpindel_annotated_withheader_Chr12.vcf


