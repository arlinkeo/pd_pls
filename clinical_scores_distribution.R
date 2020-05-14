# Load clinical scores

clinicalScores <- read.table("../clinical_data/SPSS from excel from Pat IDs with clin and MRI info_for Arlin4_small.txt", 
                  sep = '\t', header = TRUE)
rownames(clinicalScores) <- clinicalScores$ID

# Overlapping patient IDs between data of cortical thickness and clinical measures
patientIDs <- as.character(intersect(ct$ID, clinicalScores$ID))

# Select clinical features with (non-nominal) numeric values 
clinicalScores <- clinicalScores[patientIDs, c(10:13,17,22, 24:27)]

# Check available scores (non-NA's)
number_cs <- paste0(colnames(clinicalScores), "\n(n=", apply(clinicalScores, 2, function(x) sum(!is.na(x))), ")")

# Histograms of clinical scores
df <- clinicalScores
colnames(df) <- number_cs
df <- melt(df[,-1])
pdf("output/histogram_clinical_scores.pdf", 8, 5)
ggplot(df) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = "white") +
  theme_classic() +
  facet_wrap(~ variable, scales = "free")
dev.off()
