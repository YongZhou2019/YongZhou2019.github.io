## Get histograms for genes with variants on their exons
## Get histograms for genes with variants on their exons
pacman::p_load(ggplot2, cowplot)

## load data
g1 <- read.table(file="/Users/thimmamp/Magic16Vis/indels/genome1_numvariants_pertxpt.txt", header = FALSE)
colnames(g1) <- c("NumberofVariants", "Gene")
g4 <- read.table(file="/Users/thimmamp/Magic16Vis/indels/genome4_numvariants_pertxpt.txt", header = FALSE)
colnames(g4) <- c("NumberofVariants", "Gene")
g5 <- read.table(file="/Users/thimmamp/Magic16Vis/indels/genome5_numvariants_pertxpt.txt", header = FALSE)
colnames(g5) <- c("NumberofVariants", "Gene")
g6 <- read.table(file="/Users/thimmamp/Magic16Vis/indels/genome6_numvariants_pertxpt.txt", header = FALSE)
colnames(g6) <- c("NumberofVariants", "Gene")

#mean_variant <- mean(df$NumberofVariants)
# Basic histogram

#geom_histogram(ylim=c(0,100))
p1 <- ggplot(g1, aes(x=NumberofVariants)) +
  geom_histogram( color="darkblue", fill="lightblue", breaks = seq(0, 200, by = 10))
p2 <- ggplot(g4, aes(x=NumberofVariants)) +
  geom_histogram( color="darkblue", fill="lightblue", breaks = seq(0, 200, by = 10))
p3 <- ggplot(g5, aes(x=NumberofVariants)) +
  geom_histogram( color="darkblue", fill="lightblue", breaks = seq(0, 200, by = 10))
p4 <- ggplot(g6, aes(x=NumberofVariants)) +
  geom_histogram( color="darkblue", fill="lightblue", breaks = seq(0, 200, by = 10))

plot_grid(p1,p2,p3,p4,ncol = 2, nrow = 2)

#geom_vline(aes(xintercept = mean_variant, color = "red", linewidth = 0.02))

## with density
#geom_histogram(aes(y=..density..), colour="black", fill="white")+
#  geom_density(alpha=.2, fill="#FF6666")

# Change the width of bins
#ggplot(df, aes(x=NumberofVariants)) +
#  geom_histogram(binwidth=1)
# Change colors
#p<-ggplot(df, aes(x=NumberofVariants)) +
#  geom_histogram(color="black", fill="white")
