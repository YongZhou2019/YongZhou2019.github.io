library(CMplot)
#data("cattle50K")
#data("pig60K")
rm(list=ls())

# head(pig60K)
# CMplot(pig60K,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
#        main="illumilla_60K",verbose=TRUE,width=9,height=6)
# CMplot(pig60K,type="p",plot.type="c",chr.labels=paste("Chr",c(1:18,"X","Y"),sep=""),r=0.4,
#        cir.legend=TRUE, outward=FALSE, cir.legend.col="black", cir.chr.h=1.3,chr.den.col="black",
#        file="jpg", memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)
## inhouse data
#genome1 <- read.table("/home/thimmamp/MAGIC16/tmp/vcf/old/genome1VsMagic16genomes_deletionsonly_forcmplot.txt", sep="\t", header = T)
#genome <- read.table("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/genome5_GJ_trop1_Azucenamodified.filtered.deletions.only.forcmplot.txt", sep="\t")
#g1ins <- read.table("/home/thimmamp/MAGIC16/tmp/vcf/genome33_processed/IR64RS2_deletions_forcmplot.txt", 
#                    sep="\t")

snps <- read.table("/home/thimmamp/10K/data/genome1_SNPS_forcmplot.txt", sep="\t")
setwd("/home/thimmamp/10K/Visualisation/")
filename="genome1Vs10KRGP_DensityPlot"
title="genome1Vs10KRGP_genomewide_SNPs_DensityPlot"
CMplot(snps,type="p",plot.type="d",bin.size=1e6,
       #chr.den.col=c("darkgreen", "yellow", "red"),
       #chr.den.col=c("orange", "blue", "black"),
       chr.den.col=c("lightgrey", "grey", "black"),
       file="pdf", memo=filename, dpi=300,
       main=title,
       verbose=TRUE,width=9,height=6)

