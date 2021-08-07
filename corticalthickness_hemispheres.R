# Comparing the cortical thickness between hemispheres for each of the 34 cortical brain regions in PD

# Load cortical thickness data
ct <- read.table("../Corticale_dikte_PD_FS_.txt", header = TRUE)
colnames(ct) <- gsub("_thickness", "", colnames(ct))
rownames(ct) <- ct$ID

# CT differences between men and women
men <- which(ct$Sex == 1)
women <- which(ct$Sex == 0)
ttest_gender <- data.frame(t(apply(ct[,5:ncol(ct)], 2, function(x){
  m <- x[men]
  w <- x[women]
  t <- t.test(m, w)
  c(t$statistic,
    'mean difference' = unname(t$estimate[1]) - unname(t$estimate[2]),
    ct.men = unname(t$estimate[1]), 
    ct.women = unname(t$estimate[2]), 
    pvalue = t$p.value)
})))
ttest_gender$BH <- p.adjust(ttest_gender$pvalue)
ttest_gender <- cbind('Region' = rownames(ttest_gender), ttest_gender)
ttest_gender[ttest_gender$BH<0.05, ]

# Regions of interest (ROI's) per hemisphere
rois <- intersect(colnames(ct), parcel_id$name)
rois_lh <- rois[grep("lh_", rois)]
rois_rh <- rois[grep("rh_", rois)]
ct_lh <- ct[, rois_lh]
ct_rh <- ct[, rois_rh]
identical(gsub("lh_", "", rois_lh), gsub("rh_", "", rois_rh))# Order of samples is identical between hemispheres

# Comparing cortical thickness between left and right hemisphere
ttest_h <- data.frame(t(sapply(c(1:length(rois_lh)), function(i) {
  r <- ct_rh[, i]
  l <- ct_lh[, i]
  t <- t.test(r,l)
  c(t$statistic,
    'mean difference' = unname(t$estimate[1]) - unname(t$estimate[2]),
    ct.rh = unname(t$estimate[1]), 
    ct.lh = unname(t$estimate[2]), 
    pvalue = t$p.value)
})))
ttest_h$BH <- p.adjust(ttest_h$pvalue)
ttest_h <- cbind('Region' = rois_lh, ttest_h)

# Write table to create figure of delta CT
tab <- ttest_h[, -6]
colnames(tab) <- c("Region name", "T-score", "Delta CT", "Mean CT left hemisphere", "Mean CT right hemisphere", "BH-corrected P")
tab <- cbind('Region ID' = parcel_id[match(tab$Region, parcel_id$name), "id"], tab)
tab <- tab[order(tab$`Region ID`), ]
write.table(tab, file = "output/cortical_thickness_hemispheres1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# For supplementary 
tab <- tab[order(tab$`BH-corrected P`), ]
tab[, c(3:6)] <- format(round(tab[, c(3:6)], digits = 2), nsmall = 2)
signif_rows <- tab$`BH-corrected P` < 0.05
tab$`BH-corrected P` <- format(tab$`BH-corrected P`, scientific = TRUE, digits = 3)
write.table(tab, file = "output/cortical_thickness_hemispheres2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# Write table with only significant findings for in manuscript
tab <- tab[signif_rows, ]
tab <- tab[, -1]
write.table(tab, file = "output/cortical_thickness_hemispheres3.txt", quote = FALSE, row.names = FALSE, sep = "\t")
