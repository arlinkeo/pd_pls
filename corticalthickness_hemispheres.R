# Load cortical thickness data
ct <- read.table("../Corticale_dikte_PD_FS_.txt", header = TRUE)
colnames(ct) <- gsub("_thickness", "", colnames(ct))
parcel_id <- read.table("../parcel_info.txt", col.names = c("id", "name"))
parcel_id$name <- gsub("ctx-", "", parcel_id$name)
parcel_id$name <- gsub("-", "_", parcel_id$name)

# Regions of interest (ROI's)
rois <- intersect(colnames(ct), parcel_id$name)
rois_lh <- rois[grep("lh_", rois)]
rois_rh <- rois[grep("rh_", rois)]
rois2 <- gsub("lh_", "", rois_lh)
ct_lh <- ct[, rois_lh]
ct_rh <- ct[, rois_rh]
identical (gsub("lh_", "", rois_lh), gsub("rh_", "", rois_rh))# Order of samples are identical

# Comparing cortical thickness between left and right hemisphere
ttest_h <- data.frame(t(sapply(c(1:length(rois_lh)), function(i) {
  l <- ct_lh[, i]
  r <- ct_rh[, i]
  t <- t.test(l,r)
  c(ct.lh = unname(t$estimate[1]), ct.rh = unname(t$estimate[2]), pvalue = t$p.value)
})))
ttest_h$BH <- p.adjust(ttest_h$pvalue)
ttest_h <- cbind(ROI = rois2, ttest_h)
tab <- ttest_h[order(ttest_h$BH),]
write.table(tab, file = "output/corticalthickness_hemispheres.txt", sep = "\t", row.names = FALSE, quote = FALSE)
