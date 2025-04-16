## Boxplot for genomewide insertions and deletions
pacman::p_load(tidyverse, stringr, janitor, ggrepel, GGally, broom, knitr,
               ggpubr, corrplot, ggstatsplot, reshape, rstatix)
setwd("/home/thimmamp/MAGIC16/DEGforInversions/scripts/")
rm(list = ls())

data <- read.table("Indelsshorterthan50pb_output.txt", sep="\t", header = TRUE) 
myd <- melt(data, id=c("SampleName"), variable_name="Indel") %>% 
  mutate(valueinM=round(value/1000000, digits = 2))
myd %>% head()
myd %>% 
  filter(Indel %in% c("Insertions", "Deletions")) %>% 
  ggplot(., aes(x=SampleName, y=valueinM, fill=Indel, label=valueinM)) +
  geom_bar(stat="identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  labs(x="genome", y = "Number of Indels", title="Distribution of Indels shorter than50bp per genome MAGIC16") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
ggsave(filename="GenomewiseIndelDistributionMagic16.png", height=10, width=7)
    
#     
#     + scale_color_manual(values = c("#E7B800", "#FC4E07", "#00AFBB")) + 
#       geom_jitter(shape=16, position=position_jitter(0.2)) + 
#       labs(title= paste0("Percentage of DEGs for ",my.accession," from ", my.tissue)) +
#       xlab("Location") + ylab("% of DEGs") +
#       theme(plot.title = element_text(hjust = 0.5)) +
#       stat_compare_means(method = "t.test", label="p.signif", tip.length = 0.01,
#                          symnum.args=list(cutpoints=c(0,0.05,1), symbols=c("*","ns")),
#                          comparisons = list(c("Genomewide",   "FlankingOnly"),c("Genomewide", "InversionOnly")) )
#     ggsave(filename=paste0("../results/Percentage_of_DEGs_for_", my.accession, "_", my.tissue, "_100Mbwindow.png"),
#            height=10, width=7)
#   }
# }
# 
# ## for each location compared with Genomewide
# 
# for(my.tissue in unique(data$tissue)){
#   
#   for(my.accession in unique(data$accession)){
#     #     
#     #     for(my.location in unique(data$location)[2:4]) {
#     #       
#     #       mycomparison <- list( c("Genomewide", my.location) )
#     #       
#     #       df <- data %>%
#     #         filter(accession==my.accession, tissue==my.tissue) %>%
#     #         filter(location %in% c("Genomewide", my.location) )
#     #       
#     fit <-  data %>%
#       filter(accession==my.accession, tissue==my.tissue) %>%
#       #filter(location %in% c("Genomewide", my.location) ) %>%
#       t.test(value ~ location, data=.)
#     #       
#     #       p <- ggplot( df, aes(x=location, y=value, fill=location)) +
#     #         geom_boxplot(alpha=0.4, outlier.shape = 1) +
#     #         facet_grid( tissue ~ accession, scales="free_y") +
#     #         #scale_y_continuous(expand = expansion(add = c(0,1)), limits = c(0,NA)) +
#     #         theme_bw(base_size=18) +
#     #         theme(legend.position="none") +
#     #         labs(x=NULL, y = "Percentage of DEGs", caption="**** means p-value < 2.2e-16") +
#     #         stat_compare_means(label = "p.signif", comparisons = mycomparison, method="t.test") #+
#     #       #scale_fill_manual(values = c(Genomewide="green", mycomparison[[1]][2]="red") )
#     #       
#     #       ggsave(filename=paste0("Percentage_of_DEGs_for_",my.accession,"_",my.tissue,"_", my.location,  ".png"),
#     #              height=10, width=7)
#     #       #
#     cat(my.tissue, my.accession, fit$statistic, fit$p.value, "\n")
#     rm(fit)
#     #       #print(p)
#     #       
#     #     }
#   }
# }
# 
# 
# #########
# 
# 
# mycomparison <- list( c("Genomewide", "Inversion") )
# # mycomparison <- list( c("Genomewide", "All") )
# #
# data$location <- factor(data$location, levels=c("Genomewide", "FlankingOnly", "InversionOnly"))
# 
# ggplot(data, aes(x=location, y=value, fill=location)) +
#   geom_boxplot(alpha=0.4, outlier.shape = 1) +
#   facet_grid( tissue ~ accession, scales="free_y") +
#   scale_y_continuous(expand = expansion(add = c(0,1)), limits = c(0,NA)) +
#   theme_bw(base_size=18) +
#   theme(legend.position="none") +
#   labs(x=NULL, y = "Percentage of DEGs", caption="ns: p > 0.05\n**: p <= 0.01\n ***: p <= 0.001\n ****: p-value < 2.2e-16") +
#   #stat_compare_means(label = "p.signif", comparisons = mycomparison, method="t.test")+
#   stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Genomewide",   "FlankingOnly"),c("Genomewide", "InversionOnly")),
#   )
# scale_fill_manual(values = c(Genomewide="green", FlankingOnly="red", InversionOnly="blue"))
# 
# ggsave(filename="../results/Percentage_of_DEGs_Both_100Mb.png", height=10, width=7)