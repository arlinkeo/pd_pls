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

# Write table to create figure of delta CT
tab <- data.frame(ct_test[, c("Group", "t", "Mean Difference", "mean_control", "mean_pd", "BH")])
colnames(tab) <- c("Region name", "T-score", "Delta CT", "Mean CT control", "Mean CT PD", "BH-corrected P")
tab <- cbind('Region ID' = parcel_id[match(tab$Region, parcel_id$name), "id"], tab)
tab <- tab[order(tab$`Region ID`), ]
write.table(tab, file = "output/cortical_thickness_conditions1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# For supplementary 
tab <- tab[order(tab$`BH-corrected P`), ]
tab[, c(3)] <- as.numeric(tab[, c(3)])
tab[, c(4)] <- as.numeric(tab[, c(4)])
tab[, c(3:6)] <- round(tab[, c(3:6)], digits = 3)
signif_rows <- tab$`BH-corrected P` < 0.05
tab$`BH-corrected P` <- format(tab$`BH-corrected P`, scientific = TRUE, digits = 3)
write.table(tab, file = "output/cortical_thickness_conditions2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# Write table with only significant findings for in manuscript
tab <- tab[signif_rows, ]
tab <- tab[, -1]
write.table(tab, file = "output/cortical_thickness_conditions3.txt", quote = FALSE, row.names = FALSE, sep = "\t")
