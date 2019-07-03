setwd("C:/Users/dkeo/surfdrive/pd_imaging_pls")
options(stringsAsFactors = FALSE)

# Brain donor names
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# Load gene expression data 
brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")

# Sample info per donor
ontology <- read.csv("../AHBA_Arlin/Ontology.csv")
sampleInfo <- lapply(donorNames, function(d){
  read.csv(paste0("../AHBA_Arlin/sample_info_", d, "_2018-11-18.csv"))
})

# Sample selection from Arnatkeviciute et al. 2019
coord <- read.delim("sample_coordinates_roi.txt", col.names = c("sample", "x", "y", "z"))
parcels <- read.table("parcel_info.txt", col.names = c("ROI", "Name"))

# Match samples by MNI coordinates
sampleInfo <- lapply(donorNames, function(d){
  x <- sampleInfo[[d]]
  s1 <- x[, c("mni_x", "mni_y", "mni_z")]
  s1 <- paste(as.data.frame(t(s1)))
  s2 <- coord[,-1]
  s2 <- paste(as.data.frame(t(s2)))
  rows <- match(s1, s2)
  cbind(fs_roi = coord[rows,1], x)
})
saveRDS(sampleInfo, file = "resources/sampleInfo.rds")

t <- sapply(sampleInfo2, function(x) {
  t <- x[!is.na(x$fs_roi), "fs_roi"]
  sapply(setNames(c(1:34), c(1:34)), function(r) sum(t==r), USE.NAMES = TRUE)
})
colnames(t) <- gsub("donor", "Donor ", colnames(t))
t <- cbind(parcels[rownames(t),], t, Total = apply(t, 1, sum))
t$Name <- gsub("ctx-lh-", "", t$Name)
write.table(t, file = "number_of_samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)