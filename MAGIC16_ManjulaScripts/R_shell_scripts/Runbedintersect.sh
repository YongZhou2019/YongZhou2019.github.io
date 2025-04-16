#!/bin/bash
#bedtools intersect -a genome1_processed/genome1_exons.bed -b genome1_processed/genome1.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome1_processed/genome1_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome4_processed/genome4_exons.bed -b genome4_processed/genome4.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome4_processed/genome4_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome5_processed/genome5_exons.bed -b genome5_processed/genome5.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome5_processed/genome5_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome6_processed/genome6_exons.bed -b genome6_processed/genome6.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome6_processed/genome6_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome7_processed/genome7_exons.bed -b genome7_processed/genome7.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome7_processed/genome7_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome8_processed/genome8_exons.bed -b genome8_processed/genome8.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome8_processed/genome8_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome9_processed/genome9_exons.bed -b genome9_processed/genome9.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome9_processed/genome9_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome10_processed/genome10_exons.bed -b genome10_processed/genome10.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome10_processed/genome10_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome11_processed/genome11_exons.bed -b genome11_processed/genome11.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome11_processed/genome11_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome12_processed/genome12_exons.bed -b genome12_processed/genome12.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome12_processed/genome12_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome13_processed/genome13_exons.bed -b genome13_processed/genome13.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome13_processed/genome13_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome14_processed/genome14_exons.bed -b genome14_processed/genome14.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome14_processed/genome14_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome15_processed/genome15_exons.bed -b genome15_processed/genome15.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome15_processed/genome15_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome16_processed/genome16_exons.bed -b genome16_processed/genome16.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome16_processed/genome16_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome27_processed/genome27_exons.bed -b genome27_processed/genome27.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome27_processed/genome27_exons_gffVsvariant_intersect_output.txt
#bedtools intersect -a genome28_processed/genome28_exons.bed -b genome28_processed/genome28.nostar.bimultiallele.50bpindel_annotated_col1_9.bed -wo > genome28_processed/genome28_exons_gffVsvariant_intersect_output.txt

### Run exons.bed intersect with variants, and get uniq genes
#cut -f4 genome1_processed/genome1_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'_' -f1| sort | uniq -c > genome1_processed/genome1_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
#cut -f4 genome4_processed/genome4_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'_' -f1| sort | uniq -c > genome4_processed/genome4_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome5_processed/genome5_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome5_processed/genome5_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome6_processed/genome6_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome6_processed/genome6_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome7_processed/genome7_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome7_processed/genome7_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome8_processed/genome8_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome8_processed/genome8_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome9_processed/genome9_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome9_processed/genome9_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome10_processed/genome10_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome10_processed/genome10_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome11_processed/genome11_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome11_processed/genome11_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome12_processed/genome12_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome12_processed/genome12_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome13_processed/genome13_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome13_processed/genome13_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome14_processed/genome14_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome14_processed/genome14_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome15_processed/genome15_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome15_processed/genome15_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome16_processed/genome16_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome16_processed/genome16_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
#cut -f4 genome27_processed/genome27_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome27_processed/genome27_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
#cut -f4 genome28_processed/genome28_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1| sort | uniq -c > genome28_processed/genome28_exons_gffVsvariant_UniqGeneswithNumofVariants.txt

### change for genome27 and genome28
cut -f4 genome27_processed/genome27_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1,2 | cut -d'_' -f1 | uniq -c > genome27_processed/genome27_exons_gffVsvariant_UniqGeneswithNumofVariants.txt
cut -f4 genome28_processed/genome28_exons_gffVsvariant_intersect_output.txt| cut -d'=' -f2 | cut -d'.' -f1,2 | cut -d'_' -f1 | uniq -c > genome28_processed/genome28_exons_gffVsvariant_UniqGeneswithNumofVariants.txt

