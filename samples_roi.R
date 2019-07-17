
# Sample selection from Arnatkeviciute et al. 2019
coord <- read.delim("../sample_coordinates_roi.txt", col.names = c("sample", "x", "y", "z"))
parcels <- read.table("../parcel_info.txt", col.names = c("ROI", "Name"))

# Match samples by MNI coordinates
samples_roi <- lapply(donorNames, function(d){
  x <- sample_info[[d]]
  s1 <- x[, c("mni_x", "mni_y", "mni_z")]
  s1 <- paste(as.data.frame(t(s1)))
  s2 <- coord[,-1]
  s2 <- paste(as.data.frame(t(s2)))
  rows <- match(s1, s2)
  cbind(fs_roi = coord[rows,1], x)
})
# saveRDS(samples_roi, file = "output/samples_roi.rds")

t <- sapply(samples_roi, function(x) {
  t <- x[!is.na(x$fs_roi), "fs_roi"]
  sapply(setNames(c(1:34), c(1:34)), function(r) sum(t==r), USE.NAMES = TRUE)
})
colnames(t) <- gsub("donor", "Donor ", colnames(t))
t <- cbind(parcels[rownames(t),], t, Total = apply(t, 1, sum))
t$Name <- gsub("ctx-lh-", "", t$Name)
write.table(t, file = "output/number_of_samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)
