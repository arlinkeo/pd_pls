library(ggplot2)
library(reshape2)

########## Mean thickness between conditions ##########

prepare.file <- function(file){
  df <- read.delim(file)
  filter_rows <- grep("Leeftijd|Mean|bankssts", df$Group)
  df <- df[-c(filter_rows, filter_rows+1), ]
  df[seq(2, nrow(df), by=2), "Group"] <- df[seq(1, nrow(df), by=2), "Group"]
  df$hemisphere[grep("lh", df$Group)] <- "L"
  df$hemisphere[grep("rh", df$Group)] <- "R"
  df
}

plot.mean <- function(df){
  mean_ct  <- df[, c("X", "Mean", "hemisphere")]
  df <- melt(mean_ct)
  ggplot(df) + geom_boxplot(aes(x=interaction(X, hemisphere), y=value)) +
    theme_classic()
}

group_stats1 <- prepare.file('../cortical_thickness_ttest_jeroen/group_statistics_cortical_thickness_large.txt')
p1 <- plot.mean(group_stats1)
p1

group_stats2 <- prepare.file('../cortical_thickness_ttest_jeroen/group_statistics_cortical_thickness_small.txt')
p2 <- plot.mean(group_stats2)
p2

########## Equality tests between conditions ##########

eq_tests <- read.delim("../cortical_thickness_ttest_jeroen/equality_tests_cortical_thickness.txt", stringsAsFactors = FALSE)
colnames(eq_tests) <- eq_tests[1, ]
eq_tests <- eq_tests[-c(1,2),]
colnames(eq_tests)[c(1,2,10,11)] <- c("Group", "variance_assumption", "lower95", "upper95")
filter_rows <- grep("Leeftijd|Mean", eq_tests$Group)
eq_tests <- eq_tests[-c(filter_rows, filter_rows+1), ]
eq_tests[seq(2, nrow(eq_tests), by=2), "Group"] <- eq_tests[seq(1, nrow(eq_tests), by=2), "Group"]

eq2 <- eq_tests[eq_tests$variance_assumption == "Equal variances not assumed", -c(3,4)]
eq2$BH <- p.adjust(eq2$`Sig. (2-tailed)`, method = "BH")
signif_rows <- eq2$BH < 0.05

tab <- data.frame(eq2[signif_rows, c("Group", "Mean Difference", "Std. Error Difference", "BH")])
colnames(tab) <- c("Region", "Mean Difference", "Std. Error Difference", "BH-adjusted P")
tab$Region <- gsub("_", " ", tab$Region)
tab <- tab[order(tab$`BH-adjusted P`),]
tab[, 2] <- round(as.numeric(tab[, 2]), digits = 3)
tab[, 3] <- round(as.numeric(tab[,3]), digits = 3)
tab[, 4] <- format(tab[, 4], scientific = TRUE, digits = 4)
write.table(tab, file = "output/cortical_thickness_conditions.txt", quote = FALSE, row.names = FALSE, sep = "\t")
