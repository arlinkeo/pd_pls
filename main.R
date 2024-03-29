setwd("C:/Users/Arlin/surfdrive/pd_imaging_pls/pd_pls")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(reshape2)
library(pls)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(DescTools)
library(ggpubr)
library(abind)
library('ReactomePA')
library(foreign)
library(ggrepel)
library(gridExtra)
library(gplots)

# Useful variables
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# AHBA data directory and data
ahba_dir <-"C:/Users/Arlin/surfdrive/AHBA_Arlin"
probeInfo <- read.csv(paste0(ahba_dir, "/probe_info_2018-11-18.csv"))
brainExpr <- readRDS(paste0(ahba_dir, "/gene_expr.RDS"))
ontology <- read.csv(paste0(ahba_dir, "/Ontology.csv"))
sample_info <- lapply(donorNames, function(d){ # Sample info per donor
  read.csv(paste0(ahba_dir, "/sample_info_", d, "_2018-11-18.csv"))
})

# Load parcel info (Arnatkevic̆iūtė et al., 2019)
parcel_id <- read.table("../parcel_info.txt", col.names = c("id", "name"))
parcel_id$name <- gsub("ctx-", "", parcel_id$name)
parcel_id$name <- gsub("-", "_", parcel_id$name)

# Source directory with functions
fun_dir <- paste0(getwd(), "/functions/")
R.utils::sourceDirectory(fun_dir, modifiedOnly = FALSE)

# Make output folder
dir.create("output")

# Run scripts
source("corticalthickness_conditions.R")
source("corticalthickness_hemispheres.R")
source("samples_roi.R")
source("pls_model1.R")
source("clinical_scores_distribution.R")
source("ct_clinicalscores.R")
source("pls_model2.R")
