# Load clinical scores

# clin_dir <- "//vf-DataSafe/DataSafe$/div3/neur/SLEUTELBESTANDEN_794/nieuwe_export_pd"
# clin_files  <- list.files(clin_dir)
# clin_files <- clin_files[grep(".SAV", clin_files)]
# data <- lapply(paste0(clin_dir, "/", clin_files), function(x) read.spss(x, to.data.frame = TRUE))
# names(data) <- clin_files
# 
# clin_scores <- data$BA_.SAV
# 
# x <- lapply(data, function(t){
#   c <- colnames(t)
#   grep("SCOPA", c)
#   # apply(t, 2, function(x){
#   #   # length(intersect(x, clinicalScores$ID))
#   #   length(grep("3163", x))
#   # })
# })
# lapply(x, max)

clinicalScores <- read.table("../clinical_data/SPSS from excel from Pat IDs with clin and MRI info_for Arlin4_small.txt", 
                  sep = '\t', header = TRUE)
rownames(clinicalScores) <- clinicalScores$ID

# Overlapping patient IDs between data of cortical thickness and clinical measures
patientIDs <- as.character(intersect(ct$ID, clinicalScores$ID))

# Select clinical features with (non-nominal) numeric values 
clinicalScores <- clinicalScores[patientIDs, c(10:13, 15,17,22, 24:27)]

# Check available scores (non-NA's)
number_cs <- paste0(colnames(clinicalScores), "\n(n=", apply(clinicalScores, 2, function(x) sum(!is.na(x))), ")")

# Histograms of clinical scores
df <- clinicalScores
colnames(df) <- number_cs
df <- melt(df)
pdf("output/histogram_clinical_scores.pdf", 8, 5)
ggplot(df) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = "white") +
  theme_classic() +
  facet_wrap(~ variable, scales = "free")
dev.off()
