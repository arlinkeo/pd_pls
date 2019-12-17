# Use the mapping of Allen samples to the Desikan-Killiany atlas, curated by Arnatkeviciute et al. 2019
# Using our own processed Allen samples, they are matched by MNI coordinates per donor

# Sample selection from Arnatkeviciute et al. 2019 annotated with MNI-coordinates
coord <- read.delim("../sample_coordinates_roi.txt", col.names = c("sample", "x", "y", "z"))

# Match samples by MNI-coordinates for each AHBA donor and add information to sample_info
sample_info <- lapply(donorNames, function(d){
  x <- sample_info[[d]]
  s1 <- x[, c("mni_x", "mni_y", "mni_z")]
  s1 <- paste(as.data.frame(t(s1)))
  s2 <- coord[,-1]
  s2 <- paste(as.data.frame(t(s2)))
  rows <- match(s1, s2) # rows in s2
  cbind(fs_roi = coord[rows,1], x)
})
# saveRDS(samples_info, file = "output/samples_info.rds")

# Sample size per region and donor
t <- sapply(sample_info, function(x) {
  t <- x[!is.na(x$fs_roi), "fs_roi"] # roi IDs
  sapply(setNames(c(1:34), c(1:34)), function(r) sum(t==r))
})
colnames(t) <- gsub("donor", "Donor ", colnames(t))
t <- cbind(parcel_id[rownames(t),], t, Total = apply(t, 1, sum))
t$name <- gsub("lh_", "", t$name)
colnames(t)[1:2] <- c("ROI", "Name")
write.table(t, file = "output/number_of_samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)
