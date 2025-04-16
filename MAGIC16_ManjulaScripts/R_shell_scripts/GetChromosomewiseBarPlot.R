## Create bar plot for chromosome-wise ins and dels count
## reference http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization

pacman::p_load(tidyverse, MASS, reshape, reshape2, ggplot2, janitor, data.table)
rm(list = ls())
#x <- count.fields("/home/thimmamp/MAGIC16/indels_results/core_out_test1.txt", sep = "\t")
#which(x!=8)
#dat <- read.csv(file="indels_results/genome1_startlost_stopgain_stoplost_withclone_datasetsvenn.txt",sep = '\t',
#                  header = TRUE)
dat   <- fread("/home/thimmamp/MAGIC16/indels_results/core_out_test1.txt",fill = TRUE, sep = "\t", header=FALSE)
#dat <- read.table("~/MAGIC16/tmp/vcf/genome33_processed/genome33_chromosomewise_indelstats.txt",sep =' ', header = TRUE)
dim(dat)
colnames(dat) <- c("geneid", "variants", "ins", "dels", "RAPid", "cloneid", "gff_desc", "mapman_desc")
dat %>% 
  filter(ins==0, dels==0) %>% 
  count()
## 1739/24747 do not have INDELS
temp <- dat %>% 
  filter(ins!=0, dels!=0, gff_desc!="") %>% 
  dplyr::select(geneid, gff_desc, mapman_desc) %>% 
  group_by(gff_desc) %>% 
  summarise(Numdesc=n())%>% 
  filter(Numdesc > 5)


dat %>% 
  filter(Function!="no annotation") %>% 
  nrow()
#184 no annotation
# 504 clones with MF annotation


temp <- dat %>% 
  filter(Function!="no annotation") %>% 
  #roup_by(Function, CloneId,INDELType) %>%
  group_by(Function, INDELType, Events) %>% 
  summarise(NumEvents=n()) %>% 
  filter(NumEvents > 1)

temp <- dat %>% 
    filter(Function!="no annotation") %>% 
    group_by(Function, Events, INDELType) %>% 
    summarise(NumEvents=n()) %>% 
    filter(NumEvents > 1)

ggplot(temp, aes(x=Function, y=NumEvents, fill=INDELType)) +
  #geom_bar(stat="identity", position=position_dodge())+
  geom_bar(stat="identity")+
  labs(title="Functions of clones with INDELS", x="Function", y="Number of Clones")+
  #scale_fill_brewer(palette="Blues")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  
mdata <- melt(dat, id=c("genome", "chromosome")) 
colnames(mdata)[3] = "INDELS"
colnames(mdata)[4] = "Number"
#mdata <- mdata %>% rename("INDELS" = "variable" , "Number" = "value")
mdata %>% head()
d <- mdata %>% 
        filter(genome=="genome33", INDELS %in% c("insertions", "deletions")) 
pdf(file="~/MAGIC16/indels_results/genome33_Chromosomewise_INDELS_Numbers.pdf")       
ggplot(data=d, aes(x=chromosome, y=Number, fill=INDELS)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(title="Chromosomewise INDELs Number for IR64RS2", x="Chromosome", y="Number of INDELS")+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
dev.off()
