## Plot Training and Validation loss plots for ML on inversion images classification results
pacman::p_load(tidyverse, stringr, janitor, ggrepel, GGally, broom, knitr,
               ggpubr, corrplot, ggstatsplot, reshape, rstatix, cowplot)
setwd("/home/thimmamp/MAGIC16/DEGforInversions/scripts/")
rm(list = ls())

df <- read.table("/home/thimmamp/MAGIC16/ML_TrainingValidationData_28112022.txt", header = TRUE, sep = "\t")
pdf(file="../results/InversionML_TrainingValidation_Accuracy.pdf")
p1 <- ggplot(df, aes(x=epochnum, y=accuracy, color=type)) +
  geom_line() +
  labs(title= "\nTraining and Validation accuracy ") +
  xlab("epoch numbers") + ylab("accuracy") 
dev.off()

pdf(file="../results/InversionML_TrainingValidation_Loss.pdf")  
p2 <- ggplot(df, aes(x=epochnum, y=loss, color=type)) +
  geom_line() +
  labs(title= "\nTraining and Validation Loss ") +
  xlab("epoch numbers") + ylab("loss") 
dev.off()
plot_grid(p1,p2, align="h") %>% 
  ggsave(
    filename="acculoss.pdf", 
    plot = ., 
    #device = cairo_pdf, 
    width = 210, 
    height = 297, 
    units = "mm")
