# Comparing the cortical thickness between hemispheres for each of the 34 cortical brain regions in PD

# Load cortical thickness data
ct <- read.table("../Corticale_dikte_PD_FS_.txt", header = TRUE)
colnames(ct) <- gsub("_thickness", "", colnames(ct))
rownames(ct) <- ct$ID

# Regions of interest (ROI's)
rois <- intersect(colnames(ct), parcel_id$name)
rois_lh <- rois[grep("lh_", rois)]
rois_rh <- rois[grep("rh_", rois)]
ct_lh <- ct[, rois_lh]
ct_rh <- ct[, rois_rh]
# rois_lh <- gsub("lh_", "", rois_lh)
identical(gsub("lh_", "", rois_lh), gsub("rh_", "", rois_rh))# Order of samples are identical

# Comparing cortical thickness between left and right hemisphere
ttest_h <- data.frame(t(sapply(c(1:length(rois_lh)), function(i) {
  l <- ct_lh[, i]
  r <- ct_rh[, i]
  t <- t.test(l,r)
  c('mean difference' = unname(t$estimate[1]) - unname(t$estimate[2]),
    ct.lh = unname(t$estimate[1]), 
    ct.rh = unname(t$estimate[2]), 
    pvalue = t$p.value)
})))
ttest_h$BH <- p.adjust(ttest_h$pvalue)
ttest_h$pvalue <- NULL
ttest_h <- cbind(ROI = gsub("lh_", "", rois_lh), ttest_h)
tab <- ttest_h[ttest_h$BH < 0.05, ]
tab <- tab[order(tab$BH),]
tab[,c(2:4)] <- round(tab[,c(2:4)], digits = 2)
tab[,5] <- format(tab[,5], digits = 3, scientific = TRUE)
colnames(tab)[2:5] <- c('Mean Difference', 'Cortical thickness left hemisphere', 'Cortical thickness right hemisphere', 'BH-corrected P')
write.table(tab, file = "output/corticalthickness_hemispheres.txt", sep = "\t", row.names = FALSE, quote = FALSE)
