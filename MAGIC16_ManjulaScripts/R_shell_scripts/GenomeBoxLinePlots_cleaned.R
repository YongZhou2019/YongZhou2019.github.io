## Visualise genome1 vs 3k genomes INDELS randomisation
rm(list= ls())
pacman::p_load(tidyverse,ggplot2, cowplot)

g1_absentgenes <- data.table::fread("/home/thimmamp/MAGIC16/indels_results/GenesAbsentinNipponbare_presentin_othergenomes_indeldet_blasthits_irgsp_kwevi_withUniqueGoMolecularFunction.txt",
                                       sep="\t", header=TRUE)
mydata <- g1_absentgenes %>% 
              mutate(genome = sapply(strsplit(as.character(Geneid), "[._]"), `[`, 1)) %>% 
              filter(GONames!="NA") %>% 
              select(INS, DELS, genome, GONames)

mydata_long <- mydata %>% 
               pivot_longer(
                 cols =c("INS", "DELS"),
                 names_to = "INDELS_CATEGORY",
                 values_to = "INDELS"
               ) %>% 
              filter(INDELS > 0)


#genome=="Osmh63",INDELS_CATEGORY=="DELS",
for (gen in c("OsAzu", "OsCMeo", "OsKYG", "OsLima", "OsGoSa", "OsKeNa", 
              "OsNaBo", "OsARC", "OsPr106", "OsLiXu", "OsLaMu","Oszs97", "Osmh63") ) { 
  gen <- "Osmh63"
png(file=paste0("/home/thimmamp/MAGIC16/Visualisation/GenesAbsentinIRGSPpresentin_",gen, "_IndelswithMolecularFunctions.png"))
 p <- mydata_long %>% 
   filter(genome == gen) %>% 
  ggplot(., aes(x=GONames, y=INDELS, fill=INDELS_CATEGORY))+
   geom_bar(stat = "identity",position = "stack") +
   geom_text(aes(label=INDELS), vjust = 1)+
   labs(title = paste0("Bar Plot for ", gen, " with INDELS for MF"),
       x = "MF",
       y = "INDELS") +
  theme_minimal()+
   theme(axis.text.x = element_text(angle = 90))
 print(p)
 dev.off()
  #+ facet_wrap(~genome)
}
  
# permutated_snps <- read.table("/home/thimmamp/10K/Randomization/Randomisation_filtered_percentage_numsnps_biallelic.txt",
#                                                             sep=" ", header = FALSE)
#permutated_snps_newsnps <- read.table("/home/thimmamp/10K/data/genome1_basemissing_Randomisation_newsnps_numsnps.txt",
#                                      sep="\t", header = TRUE)

#tragus <- read.table("/home/thimmamp/SVillosa/references/Salsola_tragus_chr_size.txt", sep="\t", header = TRUE)

#twentyk_newsnps <- read.table("/home/thimmamp/10K/Randomization/NumofNewSnps.txt", sep="\t",
#                              header=FALSE)
# chrsnps_1mbwindow <- read.table("/home/thimmamp/10K/data/genome1_filtered_SNPs_count_1Mbwindow.txt", sep="\t")
# 
# permuted_newsnps  <- read.table("/home/thimmamp/10K/data/genome1_filtered_permutation_newsnps_all_counts_new.txt", sep="\t")


#colnames(permutated_snps) <- c("num_samples", "randomisation_num", "numsnps", "file")
  #c("numsamples", "randomisation_num", "numsnps") 
#colnames(chrsnps) <- c("chrnum", "chr", "numsnps")
# colnames(chrsnps_1mbwindow)  <- c("chr", "start", "end", "numsnps")
# colnames(permuted_newsnps)  <- c("group", "randomisation_num", "numsnps")

pdf(file=paste0("/home/thimmamp/10K/Visualization/genome1basemissingVs10K_NumSNPs.pdf"))
permutated_snps_newsnps %>% 
  #gather("snps", "numbers", -c(percentage, Randomstep)) %>% 
  ggplot(., aes(x=as.factor(percentage), y=numsnps)) + geom_boxplot()+  
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="orange", aes(group=1),
               geom="line", lwd=1, lty=1) +
  #theme_classic() +
  labs(title="Number of SNPs in genome1 basemissing Vs 20K genomes by randomization", 
       x="percentage of rice accessions sampled", y = "Number of SNPs")
dev.off()

pdf(file=paste0("/home/thimmamp/10K/Visualization/genome1basemissingVs10K_NumandnewNPsfacetwrap.pdf"))
permutated_snps_newsnps %>% 
  gather("snps", "numbers", -c(percentage, Randomstep)) %>% 
  ggplot(., aes(x=as.factor(percentage), y=numbers)) + geom_boxplot(aes(fill=snps))+  
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_point(aes(y=numbers, group=snps), position = position_dodge(width=0.75))+
  facet_wrap( ~ percentage, scales="free")+
 xlab("Percentage of samples permutated") + ylab("Number of SNPS") + ggtitle("Trend in the number of snps and new snps observed during \ndifferent percentage of samples permutated")+
 guides(fill=guide_legend(title="Legend_Title"))
dev.off()

################
permutated_snps_newsnps %>% 
  gather("snps", "numbers", -c(percentage, Randomstep)) %>%
  ggplot(., aes(x=as.factor(percentage), y=numbers)) + 
  geom_boxplot(aes(fill=snps)) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="blue", aes(group=1),
               geom="line", lwd=1, lty=1) 

### two lines in one plot with two different y-axis
percent <- permutated_snps_newsnps %>%
  select(percentage, numsnps, newsnps) %>% 
  group_by(percentage) %>% 
  summarise_all("mean") %>% 
  select(percentage)
  
avgnew  <- permutated_snps_newsnps %>%
  select(percentage, numsnps, newsnps) %>% 
  group_by(percentage) %>% 
  summarise_all("mean") %>% 
  select(newsnps)

avgnum  <- permutated_snps_newsnps %>%
  select(percentage, numsnps, newsnps) %>% 
  group_by(percentage) %>% 
  summarise_all("mean") %>% 
  select(numsnps)

data <-data.frame(percent, avgnew, avgnum)

## set up some fake test data
# time <- seq(0,72,12)
# betagal.abs <- c(0.05,0.18,0.25,0.31,0.32,0.34,0.35)
# cell.density <- c(0,1000,2000,3000,4000,5000,6000)

## add extra space to right margin of plot within frame
# par(mar=c(5, 4, 4, 6) + 0.1)
# 
# ## Plot first set of data and draw its axis
# plot(data$percent, data$avgnum, pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", 
#      type="b",col="black", main="genome1 basemissing Number and New SNPs")
# axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
# mtext("Number of new snps",side=2,line=2.5)
# box()
# 
# ## Allow a second plot on the same graph
# par(new=TRUE)
# 
# ## Plot the second plot and put axis scale on right
# plot(data$percent, data$avgnum, pch=15,  xlab="", ylab="", ylim=c(0,7000), 
#      axes=FALSE, type="b", col="red")
# ## a little farther out (line=4) to make room for labels
# mtext("Cell Density",side=4,col="red",line=4) 
# axis(4, ylim=c(0,7000), col="red",col.axis="red",las=1)
# 
# ## Draw the time axis
# axis(1,pretty(range(time),10))
# mtext("Time (Hours)",side=1,col="black",line=2.5)  
# 
# ## Add Legend
# legend("topleft",legend=c("Beta Gal","Cell Density"),
#        text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

# Draw first plot using axis y1 
par(mar = c(7, 3, 5, 4) + 0.3)               
plot(data$percent,  data$numsnps, type="l",pch=21, col = 2)   

# set parameter new=True for a new axis 
par(new = TRUE)          

# Draw second plot using axis y2 
plot(data$percent, data$newsnps, type="l",pch=21, col = 3, axes = FALSE, xlab = "", ylab = "") 

axis(side = 4, at = pretty(range(data$newsnps)))       
mtext("Number of new snps", side = 4, line = 3)

################

## Line plot
pdf(file=paste0("/home/thimmamp/10K/Visualization/genome1filteredVs10K_NumSNPs_1MbWindow_facet.pdf"))
p <- ggplot(data = chrsnps_1mbwindow, aes(x = end, y = numsnps)) + geom_line(linetype = "dashed", color="red") +
  geom_point() +
  labs(title="Chromosomewise Number of SNPs for genome1 Vs 20K genomes", 
       x="Genomic loci 1Mb window", y = "Number of SNPs")
p + facet_wrap(~chr, nrow=12)
dev.off()

## boxplot
pdf(file="/home/thimmamp/10K/Randomization/genome1filteredVs20K_10random_newsnps_boxplot.pdf")
ggplot(permutated_snps, aes(x=as.factor(num_samples), y=numsnps)) +
  geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="orange", aes(group=1),
               geom="line", lwd=1, lty=1) +
  #theme_classic() +
  labs(title="Number of SNPs in genome1 filtered Vs 20K genomes by randomization", 
       x="number of rice accessions sampled", y = "Number of SNPs")
dev.off()

###################
## 1Mbwindow line plots in batch

fns  <- list.files("/home/thimmamp/10K/data/", pattern="genome33.*_1Mbwindow.txt", full.names = TRUE )

getlineplot <- function(infile) {
  fn <- basename(infile)
  temp <- sub('\\_1Mbwindow\\.txt$', '', fn)
  outprefix <- sub("\\_SNPs\\_count", "",  temp)
  print(outprefix)
  chrsnps_1mbwindow <- read.table(infile, sep="\t")
  colnames(chrsnps_1mbwindow)  <- c("chr", "start", "end", "numsnps")
  pdf(file=paste0("/home/thimmamp/10K/Visualization/",outprefix,"NumSNPs_1MbWindow_facet.pdf"))
  p <- ggplot(data = chrsnps_1mbwindow, aes(x = end, y = numsnps)) + geom_line(linetype = "dashed", color="red") +
    geom_point() +
    labs(title=paste0("Chromosomewise Number of SNPs for",outprefix),
         x="Genomic loci 1Mb window", y = "Number of SNPs") + 
    facet_wrap(~chr, nrow=12)
  print(p)
  dev.off()
  rm(fn,temp,outprefix,chrsnps_1mbwindow,p)
}
## call function for list of files
out <- lapply(fns, getlineplot)
###################


pdf(file="/home/thimmamp/10K/20KRGP_Lineplotofnewsnps.pdf")
ggplot(twentyk_newsnps,aes(x=numsamples, y=numnewsnps, group=1)) +
  geom_line(linetype = "dashed", color="red") +
  geom_point() +
  labs(title="Number of new SNPs in 20K genomes", 
       x="Number of samples increased", y = "Number of SNPs")
dev.off()

ggplot(dat,aes(x=end, y=numsnps, group=1)) +
  geom_line(linetype = "dashed", color="red") +
  geom_point() +
  labs(title="Number of new SNPs in 20K genomes", 
       x="Number of samples increased", y = "Number of SNPs")

ggplot(chrsnps,aes(x=chr, y=numsnps))+
  geom_bar()

ggplot(twentyk_newsnps,aes(x=numsamples, y=numnewsnps, group=1)) +
  geom_line(linetype = "dashed", color="red") +
  geom_point() +
  labs(title="Number of new SNPs in 20K genomes", 
       x="Number of samples increased", y = "Number of SNPs")

mydata <- rand_indels %>% 
  #group_by(randomisation_num, sample_percent, Numofsamples) %>% 
  #summarise(INS = sum(numins), DELS = sum(numdels)) %>% 
  #mutate(percent = sample_percent*10) %>% 
  select(randomisation_num, sample_percent, Numofsamples, INS, DELS) 

write.table(mydata, file="/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_100random_indels.csv",sep = ",", row.names = FALSE)

pdf(file="/home/thimmamp/10K/Visualization/genome1_filtered_Vs20KRGP_10random_SNPs_randomisation_moreintervals_boxplot.pdf")
ggplot(permutated_snps, aes(x=as.factor(numsamples), y=numsnps, color=numsamples)) +
  geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2))+
   stat_summary(fun=mean, colour="orange", aes(group=1),
               geom="line", lwd=1, lty=1) +
  #theme_classic() +
  labs(title="Number of SNPs in genome1 filtered Vs 10K genomes by randomization", 
       x="Number of 10K data sampled", y = "Number of SNPs")
dev.off()



## boxplot
pdf(file="/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_100random_indels_randomisation_Insertion_boxplot.pdf")
ggplot(permutated_snps, aes(x=as.factor(sample_percent), y=INS, color=sample_percent)) +
  geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="orange",
               geom="point",position=position_dodge(width=0.75)) +
  stat_summary(fun=mean, colour="orange", aes(group=1),
               geom="line", lwd=1, lty=1) +
  theme_classic() +
  labs(title="Number of insertions in genome6 (IR64) by randomisation", 
       x="Percentage of 3K data sampled", y = "Number of Insertions")
#geom_boxplot(aes(fill=sample_percent))
dev.off()

pdf(file="/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_100random_indels_randomisation_Deletions_boxplot.pdf")
ggplot(mydata, aes(x=as.factor(sample_percent), y=DELS, color=sample_percent)) +
  geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="green",
               geom="point",position=position_dodge(width=0.75)) +
  stat_summary(fun=mean, colour="green", aes(group=1),
               geom="line", lwd=1, lty=1) +
  theme_classic() +
  labs(title="Number of deletions in genome6(IR64) by randomisation", 
       x="Percentage of 3K data sampled", y = "Number of Deletions")
#geom_boxplot(aes(fill=sample_percent))
dev.off()  


ggplot(indels, aes(x=PercentageofSamples, y=NumberofDELS))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(title="Number of deletions by randomisation", 
       x="Percentage of 3K data sampled", y = "Number of Deletions")+
  facet_grid(~Repetition)


pdf(file="/home/thimmamp/MAGIC16/indels_results/genome1Vs3K_indels_randomisation.pdf")
ggplot(indels, aes(x=PercentageofSamples, y=Number, fill=INDELS))+
  geom_bar(stat="identity", position = position_dodge())+
  theme_classic()+
  labs(title="Number of indels by randomisation", 
       x="Percentage of 3K data sampled", y = "Number of indels")+
  facet_grid(~Repetition)
dev.off()

ggplot(tragus, aes(x=name, y=length,fill=chromosome)) + 
  geom_bar(stat="identity")+
  theme_classic()+
  labs(title="Number of indels by randomisation", 
       x="Percentage of 3K data sampled", y = "Number of indels")+
  facet_grid(haplotype~chromosome)
