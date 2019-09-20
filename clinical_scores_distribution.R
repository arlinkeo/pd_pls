library(reshape2)

clinicalScores <- read.table("../SPSS from excel from Pat IDs with clin and MRI info_for Arlin4_small.txt", 
                  sep = '\t', header = TRUE)

# Overlapping patient IDs between data of cortical thickness and clinical measures
patientIDs <- as.character(intersect(ct$ID, clinicalScores$ID))
clinicalScores <- clinicalScores[clinicalScores$ID %in% patientIDs, 10:50]

# Check availble scores (non-NA's)
number_cs <- paste0(colnames(clinicalScores), "\n(n=", apply(clinicalScores, 2, function(x) sum(!is.na(x))), ")")

# Histograms of clinical scores
df <- clinicalScores
colnames(df) <- number_cs
df <- melt(df)
pdf("output/histogram_clinical_scores.pdf", 12, 12)
ggplot(df) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = "white") +
  theme_classic() +
  facet_wrap(~ variable, scales = "free", nrow = 7)
dev.off()
