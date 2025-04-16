## Generate heatmap
rm(list = ls())
pacman::p_load(tidyverse, patchwork, ggplot2)

df <- read.delim("/home/thimmamp/MAGIC16/indels_results/genome1_core_dispensable_specific_withRAP_MapmanAnnot_vareffects_forheatmap.txt", 
                 header = TRUE, sep = "\t")
df1 <- df %>% 
      filter(type=="core")
df_long <- pivot_longer(data = df,
                          cols = -c("NiporGrameneId", "Name", "Description", "FunctionalAnnotation", "INS", "DELS", "type"),
                          names_to = "Event",
                          values_to = "EventValues")
pdf(file="/home/thimmamp/MAGIC16/Visualisation/genomicevents_forcoredisspecific.pdf")
mine_heatmap <- ggplot(data = df_long, mapping = aes(x = Event,
                                                       y = EventValues, fill = type)) +
  geom_tile() +
  xlab(label = "GenomicEvent") + 
  facet_grid(~type)
mine_heatmap
dev.off()
