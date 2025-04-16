# R script to get stats on INDELs variant
pacman::p_load( tidyverse, data.table, janitor)
# # stop_gained <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_gained_details.txt", 
# #                           sep="\t", header = T, stringsAsFactors = F)
# 
# stop_gained <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_gained_details_MAPmanandRAPannotation.txt",
#                           sep="\t", header = T, stringsAsFactors = F)
# 
# stop_lost <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_stop_lost_details_MAPmanandRAPannotation.txt",
#                           sep="\t", header = T, stringsAsFactors = F)
# 
# start_lost <- read.table(file="/home/thimmamp/MAGIC16/indels_results/genome1_start_lost_details_MAPmanandRAPannotation.txt",
#                          sep="\t", header = T, stringsAsFactors = F)
data <- read.csv(file="/home/thimmamp/MAGIC16/indels_results/MergedAgritraitgenes_in_genome1_Dispensablewithduplicategenes.txt",
                   sep="\t", header = TRUE, stringsAsFactors = FALSE)
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

data %>% 
          select(Description, DELS) %>% 
          filter(DELS >= 10) %>% 
          group_by(Description) %>% 
          summarise(COUNT=n()) %>% 
          arrange(desc(COUNT)) %>% 
          head(20) %>% 
          ggplot(.,aes(x=Description, y=COUNT)) +
          geom_bar(stat="identity", position = "dodge") +
          ggtitle("Most abundant functions for deletions in Agri trait genes") +
           ylab("Number of times abundant functions with dels > 20 ")+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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

