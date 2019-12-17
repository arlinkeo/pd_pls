# Cortical thickness was compared between condition control and PD
# A t-test was performed in SPSS 23, here we read in the results and output a table

########## Mean thickness between conditions ##########

group_stats <- read.delim('../cortical_thickness_ttest_jeroen/group_statistics_cortical_thickness_large.txt')
filter_rows <- grep("Leeftijd|Mean", group_stats$Group)
group_stats <- group_stats[-c(filter_rows, filter_rows+1), ]
group_stats[seq(2, nrow(group_stats), by=2), "Group"] <- group_stats[seq(1, nrow(group_stats), by=2), "Group"]
group_stats$Group <- gsub("_thickness", "", group_stats$Group)

########## Equality tests between conditions (SPSS results) ##########

eq_tests <- read.delim("../cortical_thickness_ttest_jeroen/equality_tests_cortical_thickness.txt", stringsAsFactors = FALSE)
colnames(eq_tests) <- eq_tests[1, ]
eq_tests <- eq_tests[-c(1,2),]
colnames(eq_tests)[c(1,2,10,11)] <- c("Group", "variance_assumption", "lower95", "upper95")
filter_rows <- grep("Leeftijd|Mean", eq_tests$Group)
eq_tests <- eq_tests[-c(filter_rows, filter_rows+1), ]
eq_tests[seq(2, nrow(eq_tests), by=2), "Group"] <- eq_tests[seq(1, nrow(eq_tests), by=2), "Group"]
eq_tests$Group <- gsub("_thickness", "", eq_tests$Group)

ct_test <- eq_tests[eq_tests$variance_assumption == "Equal variances not assumed", -c(3,4)]
ct_test$mean_control <- group_stats[group_stats$X == "control", "Mean"]
ct_test$mean_pd <- group_stats[group_stats$X == "PD", "Mean"]
ct_test$BH <- p.adjust(ct_test$`Sig. (2-tailed)`, method = "BH")
signif_rows <- ct_test$BH < 0.05

tab <- data.frame(ct_test[signif_rows, c("Group", "Mean Difference", "mean_control", "mean_pd", "BH")])
colnames(tab) <- c("Region", "Mean Difference", "Mean control", "Mean PD", "BH-adjusted P")
tab$Region <- gsub("_", " ", tab$Region)
tab <- tab[order(tab$`BH-adjusted P`),]
tab$`Mean Difference` <- as.numeric(tab$`Mean Difference`)
tab[, c(2:4)] <- round(tab[, c(2:4)], digits = 3)
tab$`BH-adjusted P` <- format(tab$`BH-adjusted P`, scientific = TRUE, digits = 3)
tab[, "Region"] <- gsub("lh", "Left hemisphere", tab[, "Region"])
tab[, "Region"] <- gsub("rh", "Right hemisphere", tab[, "Region"])
write.table(tab, file = "output/cortical_thickness_conditions.txt", quote = FALSE, row.names = FALSE, sep = "\t")
