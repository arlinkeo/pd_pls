# PLS
library(mixOmics)

# Read cortical thickness data
ct <- read.table("Corticale_dikte_PD_FS_.txt", header = TRUE)
colnames(ct) <- gsub("_thickness", "", colnames(ct))

# Read clinical data
clinical <- read.table("SPSS from excel from Pat IDs with clin and MRI info_for Arlin4_small.txt", header = TRUE, sep = "\t")

ct_id <- ct$ID
cl_id <- clinical$ID

overlap <- intersect(ct_id, cl_id)
