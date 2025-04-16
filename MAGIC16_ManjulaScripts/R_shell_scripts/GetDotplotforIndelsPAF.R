## Test dotplot for IRGSP genome vs genome5 dot plot using pafr
library(pafr)
library(tidyverse)

setwd("/home/thimmamp/MAGIC16/tmp/pafiles/")
ali <- read_paf("genome1.genome5.best.paf")
ali
ali %>% class()
#prim_alignment <- filter_secondary_alignments(ali)
#ali_chr <- subset(ali, qname=="Chr01", tname="Chr01")
#ali_chr <- subset(ali, qname=="Chr03", tname="Chr03")
nrow(ali)
query="Azucena"
target = "IRGSP"
#dotplot(ali_chr1, label_seqs=TRUE, order_by = "qstart", xlab= query, ylab=target ) +theme_bw()
dotplot(ali, label_seqs = TRUE, xlab= query, ylab=target ) +theme_bw()

plot_synteny(ali,  q_chrom = "Chr01", t_chrom = "Chr03", centre=TRUE) + theme_bw()
ali$qname %>% unique()

plot_synteny(ali = ali, q_chrom = "Chr01", t_chrom = "Chr01", centre = TRUE, rc = TRUE)

plot_coverage(ali)

plot_coverage(ali, target = FALSE)

plot_coverage(ali, fill = "qname")+ scale_fill_brewer(palette = "Set1")
