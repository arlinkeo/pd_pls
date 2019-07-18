# Partial Least Squares

library(mixOmics)
library(ggplot2)
library(ggrepel)

# Parcel info
parcel_lh <- parcels[grep("lh-", parcels$Name), ]
parcel_lh$Name <- gsub("ctx-", "", parcel_lh$Name)
parcel_lh$Name <- gsub("-", "_", parcel_lh$Name)

# X: gene expression in cortical regions
roi_expr <- lapply(donorNames, function(d){
  e <- brainExpr[[d]]
  s <- samples_roi[[d]]
  e <- sapply(parcel_lh$ROI, function(r){ # Mean expression within ROI
    cols <- which(s$fs_roi == r)
    if (length(cols) > 1) {
      e <- e[, cols]
      apply(e, 1, mean)
    } else if (length(cols) == 1){
      e[, cols]
    } else {
      rep(NA, nrow(e))
    }
  }) # genes x samples
  colnames(e) <- parcel_lh$ROI
  t(e)
})
roi_expr <- apply(simplify2array(roi_expr), c(1,2), function(x) mean(x, na.rm = TRUE)) # Mean across donors
  
# Y: t-stats of comparing cortical thickness between PD and control,
tstats <- eq2[, c("Group", "t")]
rows <- match(parcel_lh$Name, tstats$Group)
tstats <- tstats[rows, ]
rownames(tstats) <- parcel_lh$ROI
tstats$t <- as.numeric(tstats$t)

# PLS
res <- pls(roi_expr, tstats$t, ncomp = 10, mode = "regression")
# tune <- perf(res, validation = "Mfold", folds = 10)
# plot(tune$Q2.total)
# abline(h = 0.0975)

# Plot
df <- data.frame(
  roi = parcel_lh$ROI, 
  name = gsub("lh_", "", parcel_lh$Name), 
  pls1 = res$variates$X[,1],
  tstats = tstats$t)
p <- ggplot(df, aes(pls1, tstats)) +
  geom_point(color = "blue") + geom_smooth(method = "lm") +
  geom_text_repel(aes(label = name), size = 3) +
  labs(x = "PLS1", y = "t-statistics") +
  ggtitle(paste("r =", round(cor(df$pls1, df$tstats), digits = 2))) +
  theme_classic()
pdf("output/PLS1_tstats.pdf", 8, 6)
p
dev.off()

# Top ranked genes per pls component
loadings <- res$loadings$X[,1]
loadings <- rev(sort(loadings))
topgenes <- names(loadings)[1:100]
df <- data.frame(
  'Gene_name' = entrezId2Name(topgenes), 
  'Gene_ID' = topgenes, 
  'PLS1_loadings' = round(loadings[topgenes], digits = 4))
df <- cbind(df[1:50, ], df[51:100,])
write.table(df, "output/top100genes_pls1.txt", quote = FALSE, row.names = FALSE, sep = "\t")

library("RDAVIDWebService")

# RDavid
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
setTimeOut(david, 200000)
bg <- addList(david, ahba.genes(), idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg

output_dir <- "output/Functional_enrichment/"
dir.create(output_dir)

result <- addList(david, topgenes, idType = "ENTREZ_GENE_ID", listName = "top100_pls1", listType = "Gene")
print(result)
setCurrentBackgroundPosition(david, 1)
getFunctionalAnnotationChartFile(david, paste0(output_dir, "top100_pls1_goterms.txt"), threshold=0.05, count=2L)
getClusterReportFile(david, paste0(output_dir, "top100_pls1_termclusters.txt"), type = c("Term"))

# Benjamini-corrected GO-terms
t <- read.csv(paste0(output_dir, "top100_pls1_goterms.txt"), header = TRUE, sep = "\t", colClasses = "character")
rows <- t$Benjamini < 0.05
t <- t[rows, c("Term", "Benjamini")]
t$Benjamini <- as.numeric(t$Benjamini)
t <- t[order(t$Benjamini), ]
t$Term <- sapply(strsplit(t$Term, "~"), function(x) x[2])
t$Benjamini <- format(t$Benjamini, digits = 2, scientific = TRUE)
write.table(t, file = paste0(output_dir, "top100_pls1_goterms_BHcorrected.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
