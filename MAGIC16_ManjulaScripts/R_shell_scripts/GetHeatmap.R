pacman::p_load(tidyverse, ggplot2, RColorBrewer,pheatmap)
# dat <- read.table("/home/thimmamp/MAGIC16/annotation/MAGIC16_genomes_orthologs_final_forheatmap.txt", 
#                   header = TRUE, sep = "\t")
dat <- read.delim("/home/thimmamp/MAGIC16/indels_results/genome1_core_dispensable_specific_withRAP_MapmanAnnot_vareffects_forheatmap.txt", 
                 header = TRUE, sep = "\t")
dat <- dat %>% 
       filter(type=="specific", FunctionalAnnotation!=" no annotation")
#%>% column_to_rownames("Type")
df <- dat[,7:10]
d <- as.matrix(df)
row.names(d) <- dat[,1]
my_group <- sort(as.numeric(as.factor(dat[,11])))
colSide <- brewer.pal(3, "Set1")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(16)
pdf(file="/home/thimmamp/MAGIC16/Visualisation/GenesVsGenomicEvent_specific_heatmap_testing.pdf")
p <- pheatmap(d, Rowv = NA, Colv = NA, 
        RowSideColors=colSide, col=colMain )
p
dev.off()
