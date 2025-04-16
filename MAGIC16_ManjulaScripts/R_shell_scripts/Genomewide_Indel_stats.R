# R script to get stats on INDELs variant
pacman::p_load( tidyverse, data.table, janitor)
# stop_gained <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_gained_details.txt", 
#                           sep="\t", header = T, stringsAsFactors = F)

stop_gained <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_gained_details_MAPmanandRAPannotation.txt",
                          sep="\t", header = T, stringsAsFactors = F)

stop_lost <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_lost_details_MAPmanandRAPannotation.txt",
                          sep="\t", header = T, stringsAsFactors = F)

start_lost <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_start_lost_details_MAPmanandRAPannotation.txt",
                         sep="\t", header = T, stringsAsFactors = F)

# ann <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_INS_Annotation_Summary.txt",
#                   header=T, stringsAsFactors = F)

# colnames(ann) <- c("Variant", "INDEL", "Location", "annotation")
# ins_summary <- 
#   ann %>% 
#   #filter(INDEL=="insertion") %>% 
#   group_by(INDEL,Location ) %>% 
#   count() #%>% 
#   ggplot( ins_summary, aes(x=Location,y=Total) ) + geom_bar(stat="identity", position = "dodge")
#stop_gained_gene_vars <- 
pdf("/home/thimmamp/MAGIC16/indels_results/stopgained_insertions_top10Funcs.pdf")
stop_gained %>%
  filter(RAPGeneSymbol!='NA', Function!=" no annotation", INDEL_classification=="insertion") %>% 
  select(RAPGeneSymbol, Function) %>% 
  group_by(Function) %>% 
  summarise(COUNT=n()) %>% 
  arrange(desc(COUNT)) %>% 
  head(10) %>% 
  ggplot(.,aes(x=Function, y=COUNT,fill=Function))+
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("Most abundant functions for insertions \n by stop codon gained") +
   ylab("Number of times found")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("/home/thimmamp/MAGIC16/indels_results/stopgained_deletions_top10Funcs.pdf")
stop_gained %>%
  filter(RAPGeneSymbol!='NA', Function!=" no annotation", INDEL_classification=="deletion") %>% 
  select(RAPGeneSymbol, Function) %>% 
  group_by(Function) %>% 
  summarise(COUNT=n()) %>% 
  arrange(desc(COUNT)) %>% 
  head(10) %>% 
  ggplot(.,aes(x=Function, y=COUNT,fill=Function))+
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("Most abundant functions for deletions \n by stop codon gained") +
  ylab("Number of times found")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

######
stop_lost %>%
  filter(RAPGeneSymbol!='NA', Function!=" no annotation", INDEL_classification=="insertion") %>% 
  select(RAPGeneSymbol, Chromosome) %>%
  group_by(RAPGeneSymbol, Chromosome) %>% 
  summarise(COUNT=n()) %>% 
  arrange(desc(COUNT)) %>% 
  head(10) %>% 
  write.table(., file = "/home/thimmamp/MAGIC16/indels_results/stoplost_top10Genes_withinsertions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######

start_lost %>%
  filter(GeneDescription!="NA", RAPGeneSymbol!='NA', Function!=" no annotation", INDEL_classification=="deletion") %>% 
  select(Chromosome, RAPGeneSymbol, Function, GeneDescription) %>% 
  group_by(Chromosome, RAPGeneSymbol, Function, GeneDescription) %>% 
  summarise(COUNT=n()) %>% 
  filter(COUNT>=5) %>% 
  arrange(Chromosome, desc(COUNT)) %>% 
  ggplot(.,aes(x=Function, y=COUNT,fill=Function))+
  geom_bar(stat="identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

stop_gained %>% 
  group_by(gene, INDEL_classification) %>% 
  select(featureid, gene, variant, Chromosome, INDEL_classification, Function, Cellular_Component) %>% 
  arrange(gene,INDEL_classification) %>% 
  tail()

pdf("/home/thimmamp/MAGIC16/indels_results/GenefuncsForGenesinChrwithmorethan15Variants.pdf")
stop_gained %>% 
  filter(Function!=" no annotation") %>% 
  group_by(Chromosome, gene, INDEL_classification, Function) %>% 
  count() %>% 
  rename(NumVariants = n) %>% 
  filter(NumVariants > 15)%>% 
  ggplot(.,aes(x=Chromosome, y=NumVariants,fill=Function))+
  geom_bar(stat="identity", position = "dodge") 
dev.off()
#%>% 
#  rename(NumVariants = n) 

stop_gained_gene_vars %>% 
  arrange(desc(NumVariants))

stop_gained_feature_classfn_summary <- stop_gained %>% 
  #filter(INDEL_classification %in% c("insertion", "deletion")) %>% 
  group_by(featureid, INDEL_classification) %>% 
  count()

stop_gained_feature_classfn_summary <-  stop_gained_feature_classfn_summary %>% 
                                          rename(Number = n) 

stop_gained_feature_classfn_summary %>% 
  group_by(featureid) %>% 
  arrange(desc(Number))
#####  

 stop_gain_func_summary <- stop_gained %>% 
  group_by(INDEL_classification, Function) %>% 
  count() 

   #filter(INDEL_classification %in% c("insertion", "deletion")) %>% 
   group_by(INDEL_classification,Function) %>% 
   count()
 
  stop_gain_func_summary <-  stop_gain_func_summary %>% 
   rename(Total = n) 
  
  stop_gain_func_summary <- stop_gain_func_summary %>% 
    filter(Total > 100)
  #######
  stop_lost_summary <- stop_lost %>% 
    filter(INDEL_classification %in% c("insertion", "deletion")) %>% 
    group_by(Chromosome, INDEL_classification) %>% 
    count()
  
  stop_lost_summary <-  stop_lost_summary %>% 
    rename(Total = n)
  
  #######
  start_lost_summary <- start_lost %>% 
    filter(INDEL_classification %in% c("insertion", "deletion")) %>% 
    group_by(Chromosome, INDEL_classification) %>% 
    count()
  
  start_lost_summary <-  start_lost_summary %>% 
    rename(Total = n)
  
  # tofile <- ins_summary %>% 
  #  arrange(desc(Sum))
  # 
  pdf(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_gained_barplot.pdf")
  ggplot(stop_gain_summary, aes(x=Chromosome, y=Total,fill=INDEL_classification)) + 
    geom_bar(stat="identity", position = "dodge") +
    ggtitle("Chromosome wise Number of exons with INDELS\n for stop_codon_gained") +
    xlab("Chromosomes") + ylab("Number")
  ggplot(stop_gain_func_summary, aes(x=INDEL_classification, y=Total,fill=Function)) + 
    geom_bar(stat="identity", position = "dodge") +
    ggtitle("Chromosome wise Number of exons with INDELS\n for stop_codon_gained") +
    xlab("Chromosomes") + ylab("Number")
  dev.off()
  pdf(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_lost_barplot.pdf")
  ggplot(stop_lost_summary, aes(x=Chromosome, y=Total,fill=INDEL_classification)) + 
    geom_bar(stat="identity", position = "dodge") +
    ggtitle("Chromosome wise Number of exons with INDELS\n for stop_codon_lost") +
    xlab("Chromosomes") + ylab("Number")
  dev.off()
  pdf(file="/home/thimmamp/MAGIC16/indels_results/genome1_start_lost_barplot.pdf")
  ggplot(start_lost_summary, aes(x=Chromosome, y=Total,fill=INDEL_classification)) + 
    geom_bar(stat="identity", position = "dodge") +
    ggtitle("Chromosome wise Number of exons with INDELS\n for start_codon_lost") +
    xlab("Chromosomes") + ylab("Number")
  dev.off()

  write.table(tofile,file="/home/thimmamp/MAGIC16/indels_results/genome1_INS_Annotation_Summary_afterupdate.txt", sep = "\t",
            quote = F, row.names = F)
