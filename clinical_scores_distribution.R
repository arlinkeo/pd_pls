setwd("C:/Users/dkeo/surfdrive/pd_imaging")
library(ggplot2)
library(reshape2)

tab <- read.table("SPSS from excel from Pat IDs with clin and MRI info_for Arlin4_small.txt", 
                  sep = '\t', header = TRUE)
colnames(tab)

hist(tab$SENSPDSC, breaks = 50)
hist(tab$MDS_UPDRS_3, breaks = 50)
hist(tab$MDS_UPDRS_3_converted, breaks = 50)

hist(tab$MMSE, breaks = 50)

t <- tab[, -c(1:8, 61:64)]
df <- melt(t, na.rm = TRUE)
ggplot(df) + geom_histogram(aes(x=value), bins = 20) +
  facet_grid(variable, scales = "free", space = "free")
