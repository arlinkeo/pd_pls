setwd("C:/Users/dkeo/surfdrive/pd_imaging_pls/pd_pls")
options(stringsAsFactors = FALSE)

# Useful variables
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# AHBA data directory and data
ahba_dir <-"C:/Users/dkeo/surfdrive/AHBA_Arlin"
probeInfo <- read.csv(paste0(ahba_dir, "/probe_info_2018-11-18.csv"))
brainExpr <- readRDS(paste0(ahba_dir, "/gene_expr.RDS"))
ontology <- read.csv(paste0(ahba_dir, "/Ontology.csv"))
sample_info <- lapply(donorNames, function(d){ # Sample info per donor
  read.csv(paste0(ahba_dir, "/sample_info_", d, "_2018-11-18.csv"))
})

# Source directory with functions
fun_dir <- paste0(getwd(), "/functions/")
R.utils::sourceDirectory(fun_dir, modifiedOnly = FALSE)

# Make output folder
dir.create("output")

# Run scripts
source("corticalthickness_conditions.R")
source("corticalthickness_hemispheres.R")
source("samples_roi.R")
source("pls.R")


source("clinical_scores_distribution.R")