## Boxplot for random sampled inversions
pacman::p_load(tidyverse, stringr, janitor, ggrepel, GGally, broom, knitr,
               ggpubr, corrplot, ggstatsplot, reshape, rstatix)
setwd("/Users/thimmamp/Magic16Vis/scripts/")
rm(list = ls())

data <- read.table("../data/Alltissue1000sampling_new.txt", sep="\t", header = TRUE) %>%
  pivot_longer(cols=everything(), names_to = "tissueloc") %>%
  mutate(tissueloc = str_replace(tissueloc, "Genomewide", "-Genomewide-"),
         #tissueloc = str_replace(tissueloc, "Inversion", "-Inversion-")) %>%
         #separate(tissueloc, c("accession", "location", "tissue"), remove = FALSE)
         tissueloc = str_replace(tissueloc, "All", "-All-"),
         tissueloc = str_replace(tissueloc, "FlankingOnly", "-FlankingOnly-"),
         tissueloc = str_replace(tissueloc, "InversionOnly", "-InversionOnly-") ) %>%
  separate(tissueloc, c("accession", "location", "tissue"), remove = FALSE)
#%>%
#filter(accession=="MH", tissue=="Root")

# compare_means(value ~ location,
#               data = data,
#               group.by = c("accession", "tissue"),
#               method = "t.test")

# for(my.tissue in unique(data$tissue)){
#  
#   for(my.accession in unique(data$accession)){
#    
#     #for(my.location in unique(data$location)) {
#    
#     fit <-  data %>%
#       filter(accession==my.accession, tissue==my.tissue) %>%
#       #filter(my.location in c("Genomewide", in.location))
#       t.test(value ~ location, data=.)
#    
#     cat(my.tissue, my.accession, fit$statistic, fit$p.value, "\n")
#     rm(fit)
#    
#   }
# }

# ML MH -20.42082 4.790713e-84
# ML ZS -49.0099 1.403146e-321
# Panicle MH -98.39163 0
# Panicle ZS -55.63897 0
# Root MH -27.87929 2.143517e-144
# Root ZS 20.00735 5.81588e-78
# YL MH -66.55441 0
# YL ZS -55.70455 0

#########

for(my.tissue in unique(data$tissue)){
  
  for(my.accession in unique(data$accession)){
    
    for(my.location in unique(data$location)[2:4]) {
      
      mycomparison <- list( c("Genomewide", my.location) )
      
      df <- data %>%
        filter(accession==my.accession, tissue==my.tissue) %>%
        filter(location %in% c("Genomewide", my.location) )
      
      fit <-  data %>%
        filter(accession==my.accession, tissue==my.tissue) %>%
        filter(location %in% c("Genomewide", my.location) ) %>%
        t.test(value ~ location, data=.)
      
      p <- ggplot( df, aes(x=location, y=value, fill=location)) +
        geom_boxplot(alpha=0.4, outlier.shape = 1) +
        facet_grid( tissue ~ accession, scales="free_y") +
        #scale_y_continuous(expand = expansion(add = c(0,1)), limits = c(0,NA)) +
        theme_bw(base_size=18) +
        theme(legend.position="none") +
        labs(x=NULL, y = "Percentage of DEGs", caption="**** means p-value < 2.2e-16") +
        stat_compare_means(label = "p.signif", comparisons = mycomparison, method="t.test") #+
      #scale_fill_manual(values = c(Genomewide="green", mycomparison[[1]][2]="red") )
      
      ggsave(filename=paste0("Percentage_of_DEGs_for_",my.accession,"_",my.tissue,"_", my.location,  ".png"),
             height=10, width=7)
      #
      #cat(my.tissue, my.accession, fit$statistic, fit$p.value, "\n")
      #rm(fit)
      #print(p)
      
    }
  }
}


#########


#mycomparison <- list( c("Genomewide", "Inversion") )
# mycomparison <- list( c("Genomewide", "All") )
#
# ggplot(data, aes(x=location, y=value, fill=location)) +
#   geom_boxplot(alpha=0.4, outlier.shape = 1) +
#   #facet_grid( tissue ~ accession, scales="free_y") +
#   scale_y_continuous(expand = expansion(add = c(0,1)), limits = c(0,NA)) +
#   theme_bw(base_size=18) +
#   theme(legend.position="none") +
#   labs(x=NULL, y = "Percentage of DEGs", caption="**** means p-value < 2.2e-16") +
#   stat_compare_means(label = "p.signif", comparisons = mycomparison, method="t.test")
# #+
#   #scale_fill_manual(values = c(Genomewide="green", mycomparison[2]="red"))
#
# ggsave(filename="Percentage_of_DEGs.png", height=10, width=7)