## Visualise genome1 vs 3k genomes INDELS randomisation
rm(list= ls())
pacman::p_load(tidyverse,ggplot2)
twentyk_snps <- read.table("/home/thimmamp/10K/Randomization/Randomisation_numsample_all_intervals.txt", sep=" ",
                         header=FALSE)
twentyk_newsnps <- read.table("/home/thimmamp/10K/Randomization/NumofNewSnps.txt", sep="\t",
                           header=FALSE)
#tenk_snps <- read.table("/data_b/IRGSP_SNP_InDels_Venn/data/10KRGP_SNPs_Randomisation.txt", sep = " ", header = FALSE)
#rand_indels <- read.table("/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_indels_100rand_allchr_output.txt", sep="\t", header = FALSE)
#colnames(rand_indels) <- c("sample_percent", "randomisation_num", "Numofsamples", "chromosome", "numins", "numdels")
colnames(twentyk_snps) <- c("numsamples", "randomisation_num", "numsnps") 
colnames(twentyk_newsnps) <- c("numsamples", "numnewsnps") 

pdf(file="/home/thimmamp/10K/20KRGP_Lineplotofnewsnps.pdf")
ggplot(twentyk_newsnps,aes(x=numsamples, y=numnewsnps, group=1)) +
  geom_line(linetype = "dashed", color="red") +
  geom_point() +
  labs(title="Number of new SNPs in 20K genomes", 
                    x="Number of samples increased", y = "Number of SNPs")
dev.off()

mydata <- rand_indels %>% 
  #group_by(randomisation_num, sample_percent, Numofsamples) %>% 
  #summarise(INS = sum(numins), DELS = sum(numdels)) %>% 
  #mutate(percent = sample_percent*10) %>% 
  select(randomisation_num, sample_percent, Numofsamples, INS, DELS) 

write.table(mydata, file="/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_100random_indels.csv",sep = ",", row.names = FALSE)
  
pdf(file="/data_b/IRGSP_SNP_InDels_Venn/20KRGP_10random_SNPs_randomisation_moreintervals_boxplot.pdf")
ggplot(twentyk_snps, aes(x=as.factor(numsamples), y=numsnps, color=numsamples)) +
  geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000), 
                     labels = expression(10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8),
                     expand = c(0, 0))+
  stat_summary(fun=mean, colour="orange", aes(group=1),
               geom="line", lwd=1, lty=1) +
  #theme_classic() +
  labs(title="Number of SNPs in 10K genomes by randomization", 
       x="Number of 10K data sampled", y = "Number of SNPs")
dev.off()



## boxplot
pdf(file="/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_100random_indels_randomisation_Insertion_boxplot.pdf")
ggplot(mydata, aes(x=as.factor(sample_percent), y=INS, color=sample_percent)) +
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
  
#indels <- read.table("genome1_INDELS_randomisation.csv", header = TRUE, sep = ",") 
#%>% 
#       pivot_longer(!c(PercentageofSamples, Repetition,NumberofSamples), names_to = "INDELS",values_to = "Number")


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
