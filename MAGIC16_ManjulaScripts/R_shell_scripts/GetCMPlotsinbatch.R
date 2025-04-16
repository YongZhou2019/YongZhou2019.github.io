library(CMplot)
rm(list=ls())
setwd("/home/thimmamp/10K/Visualization/")
#args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]

## inhouse data
fns  <- list.files("/home/zhouy0e/10K-RGP/venn_plots/", pattern="*_forCMplot.txt", full.names = TRUE )
#infile <- "/home/thimmamp/10K/data/genome1_filtered_forCMplot.txt"

  
getcmplot <- function(infile) {
  fn <- basename(infile)
  temp <- sub('\\_forCMplot\\.txt$', '', fn)
  outprefix <- sub("\\.genomewide\\.SNPs\\.withID\\.PlinkFormat", "",  temp)
  print(outprefix)
  snps <- read.table(infile, sep="\t")
  title=paste0(outprefix,"Vs10KRGP_SNPs_DensityPlot")
  CMplot(snps,type="p",plot.type="d",bin.size=1e6,
         #chr.den.col=c("darkgreen", "yellow", "red"),
         #chr.den.col=c("orange", "blue", "black"),
         chr.den.col=c("lightgrey", "grey", "black"),
         file="pdf", memo=outprefix, dpi=300,
         main=title,
         verbose=TRUE,width=12,height=8)
}
#getcmplot(filename)
out <- lapply(fns, getcmplot)
