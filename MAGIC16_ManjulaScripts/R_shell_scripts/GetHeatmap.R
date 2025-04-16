pacman::p_load(tidyverse, ggplot2, RColorBrewer,pheatmap)
dat <- read.table("/home/thimmamp/MAGIC16/annotation/MAGIC16_genomes_orthologs_final_forheatmap.txt", 
                  header = TRUE, sep = "\t")
#%>% column_to_rownames("Type")
df <- dat[,2:ncol(dat)]
d <- as.matrix(df)
row.names(d) <- dat[,1]
my_group <- sort(as.numeric(as.factor(dat[,1])))
colSide <- brewer.pal(3, "Set1")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(16)
heatmap(d, Rowv = NA, Colv = NA, 
        RowSideColors=colSide, col=colMain )
