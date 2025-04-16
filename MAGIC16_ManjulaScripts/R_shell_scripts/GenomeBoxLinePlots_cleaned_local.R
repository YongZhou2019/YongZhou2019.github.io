## Visualise genome1 vs 3k genomes INDELS randomisation
rm(list = ls())
pacman::p_load(tidyverse, patchwork)

df <- read.delim("/home/thimmamp/10K/data/raw_Randomisation_percentage_newsnps_numsnps_no3_5.txt")
          
#/home/thimmamp/10K/data/basemissing_Randomisation_newsnps_numsnps_no3_5.txt
#/home/thimmamp/10K/data/base_Randomisation_newsnps_numsnps_no3_5.txt
#/home/thimmamp/10K/data/genome1_basemaf_Randomisation_newsnps_numsnps_no3_5.txt
#/home/thimmamp/10K/data/genome1_biallele_Randomisation_newsnps_numsnps_no3_5.txt
#/home/thimmamp/10K/data/raw_Randomisation_percentage_newsnps_numsnps_no3_5.txt
#----------

#genome1_biallele_Randomisation_newsnps_numsnps.txt")
#genome1_base_Randomisation_newsnps_numsnps.txt")
#genome1_raw_Randomisation_percentage_newsnps_numsnps.txt")
#genome1_basemissing_Randomisation_newsnps_numsnps.txt")

# df %>% 
#   group_by(percentage) %>% 
#   summarize(mean(newsnps), mean(numsnps))

# df <- read.table("/home/thimmamp/10K/data/genome1_basemissing_Randomisation_newsnps_numsnps.txt",
#                   sep="\t", header = TRUE)
options(scipen=999)
options(repr.plot.width =5, repr.plot.height =9) 
library(scales)

pdf(file="/home/thimmamp/10K/Randomization/genome1_raw_NewSNPs_NumberofSNPS.pdf", width = 3, height=4)
g1 <- ggplot(df, aes(x=as.factor(percentage), y=newsnps)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="blue", aes(group=1),
               geom="line", lwd=1, lty=1) +
  theme_classic() +
  scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6))+
  #scale_y_continuous(labels = scales::comma)+
  labs(title="Permutation Test", 
       x="percentage of samples permuted", y = "Number of new SNPs")

g2 <- ggplot(df, aes(x=as.factor(percentage), y=numsnps)) + geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=mean, colour="orange", aes(group=1),
               geom="line", lwd=1, lty=1) +
  theme_classic() +
  scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6))+
  labs(x="percentage of samples permuted", y = "Number of SNPs")
  #labs(title="Permutation Test", 
  #     x="percentage of samples permuted", y = "Number of SNPs")
g1 / g2
dev.off()
#geom_dotplot(binaxis = "y", stackdir = "center")


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



pdf(file=paste0("/home/thimmamp/10K/Visualization/genome1basemissingVs10K_NumSNPs.pdf"))
