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

# Match samples by MNI coordinates
samples <- sapply(sampleInfo, function(s1){
  s1 <- s1[, c("mni_x", "mni_y", "mni_z")]
  s1 <- paste(as.data.frame(t(s1)))
  s2 <- coord[,-1]
  s2 <- paste(as.data.frame(t(s2)))
  
  rows <- which(s2 %in% s1)
  m
})


