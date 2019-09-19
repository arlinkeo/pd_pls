# Partial Least Squares
library(pls)
library(ggrepel)
library(ComplexHeatmap)

# Parcel info
parcel_lh <- parcel_id[grep("lh_", parcel_id$name), ]

# X: gene expression in cortical regions
roi_expr <- lapply(donorNames, function(d){
  e <- brainExpr[[d]]
  s <- sample_info[[d]]
  e <- sapply(parcel_lh$id, function(r){ # Mean expression within ROI
    cols <- which(s$fs_roi == r)
    if (length(cols) > 1) {
      apply(e[, cols], 1, mean)
    } else if (length(cols) == 1){
      e[, cols]
    } else {
      rep(NA, nrow(e))
    }
  }) # genes x samples
  colnames(e) <- parcel_lh$id
  t(e)
})
roi_expr <- apply(simplify2array(roi_expr), c(1,2), function(x) mean(x, na.rm = TRUE)) # Mean across donors
  
# Y: t-stats of comparing cortical thickness between PD and control,
tscores <- ct_test[grep("lh_", ct_test$Group), c("Group", "t")]
rownames(tscores) <- parcel_lh$id
tscores$t <- as.numeric(tscores$t)

# PLS
pls1 <- plsr(tscores$t ~ roi_expr, scale = TRUE)

# PLS
# pls1 <- pls(roi_expr, tscores$t, ncomp = 1, mode = "regression")
# # tune <- perf(pls1, validation = "loo") # Leave-one-out cross validation
# # plot(tune$Q2.total)
# # abline(h = 0.0975)
# pls1_scores <- pls1$variates$X[,1]
# pls1_scores <- rev(sort(pls1_scores))

# Plot
df <- data.frame(
  roi = parcel_lh$id, 
  name = gsub("lh_", "", parcel_lh$name), 
  pls1 = pls1$scores[,1],
  tscores = tscores$t)
p <- ggplot(df, aes(pls1, tscores)) +
  geom_point(color = "blue", size = 2) + geom_smooth(method = "lm") +
  geom_text_repel(aes(label = name), size = 3, force = 2) +
  labs(x = "PLS1 score", y = expression("T-score of "~Delta~"CT")) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  ggtitle(paste("r =", round(cor(pls1$scores[,1], tscores$t), digits = 2))) +
  theme_classic()
pdf("output/PLS1_tscores.pdf", 8, 5)
p
dev.off()

# Top ranked genes for PLS1 of X
pls1_coef <- sort(pls1$coefficients[,1,1], decreasing = TRUE)
topgenes <- names(pls1_coef[1:50])

# df <- data.frame(
#   'Gene_name' = entrezId2Name(topgenes), 
#   'Gene_ID' = topgenes, 
#   'PLS1_loadings' = round(loadings[topgenes], digits = 4))
# df <- cbind(df[1:50, ], df[51:100,])
# write.table(df, "output/top100genes_pls1.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Heatmap top 50 genes
exprTopGenes <- scale(roi_expr[, topgenes])
colnames(exprTopGenes) <- entrezId2Name(colnames(exprTopGenes))
rownames(exprTopGenes) <- gsub("lh_", "", rois_lh)
exprTopGenes <- exprTopGenes[order(-pls1$scores[,1]), ]
hm <- Heatmap(exprTopGenes, name = 'Z-Score\nexpression',
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_side = c("left"),
              row_title_rot = 0,
              column_names_side = c("top"),
              column_names_gp = gpar(fontsize = 10, fontface = "italic"),
              column_names_rot = 45,
              width = unit(ncol(exprTopGenes)*.8, "lines"), 
              height = unit(nrow(exprTopGenes)*.8, "lines")
)
pdf("output/heatmap_pls1_topgenes.pdf", 10.8, 6.4)
hm
dev.off()

# # Functional enrichment with Reactome
# library('ReactomePA')
# 
# reactome <- enrichPathway(topgenes, universe = ahba.genes(), readable = TRUE)
# p <- reactome@result
# p <- p[-c(1,3:5,7,8)]
# p$p.adjust <- format(p$p.adjust, digits = 3, scientific = TRUE)
# colnames(p) <- c("Pathway", "BH", "Gene count")
# write.table(p, file = paste0("output/reactomePA_", n, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
# 
# y <- gsePathway(pls1_coef, nPerm=10000,
#                 pvalueCutoff=0.2,
#                 pAdjustMethod="BH", verbose=FALSE)
# res <- as.data.frame(y)
# head(res)

# library("RDAVIDWebService")
# # RDavid
# david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
#                            url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
# setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
# setTimeOut(david, 200000)
# bg <- addList(david, ahba.genes(), idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
# bg
# 
# output_dir <- "output/Functional_enrichment/"
# dir.create(output_dir)
# 
# result <- addList(david, topgenes, idType = "ENTREZ_GENE_ID", listName = "top100_pls1", listType = "Gene")
# print(result)
# setCurrentBackgroundPosition(david, 1)
# getFunctionalAnnotationChartFile(david, paste0(output_dir, "top100_pls1_goterms.txt"), threshold=0.05, count=2L)
# getClusterReportFile(david, paste0(output_dir, "top100_pls1_termclusters.txt"), type = c("Term"))
# 
# # Benjamini-corrected GO-terms
# t <- read.csv(paste0(output_dir, "top100_pls1_goterms.txt"), header = TRUE, sep = "\t", colClasses = "character")
# rows <- t$Benjamini < 0.05
# t <- t[rows, c("Term", "Benjamini")]
# t$Benjamini <- as.numeric(t$Benjamini)
# t <- t[order(t$Benjamini), ]
# t$Term <- sapply(strsplit(t$Term, "~"), function(x) x[2])
# t$Benjamini <- format(t$Benjamini, digits = 2, scientific = TRUE)
# write.table(t, file = paste0(output_dir, "top100_pls1_goterms_BHcorrected.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
