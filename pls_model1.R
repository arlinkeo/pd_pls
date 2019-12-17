# PLS model 1: single response variable
# Gene expression as predictors and T-score of cortical thickness as response

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
x <- roi_expr
y <- tscores$t
pls1 <- plsr(y ~ x, ncomp = 10, validation = "LOO")
pdf("output/pls1_CV_LOO.pdf", 6, 4)
plot(RMSEP(pls1), legendpos = "topright")
dev.off()
summary(pls1)
explVar <- explvar(pls1)
plot(pls1, "loadings", comps = 1:2, legendpos = "topleft", xlab = "gene")
paste(gsub("Comp ", "PLS", names(explVar)), paste0(round(explVar, digits = 2), "%"), sep = ": ", collapse = "; ")
# ncomp.perm <- selectNcomp(pls1, method = "randomization", plot = TRUE)

pls1_scores <- pls1$scores[,1:3]
pls1_scores <- cbind(ID = rownames(pls1_scores), label = rois_lh, pls1_scores)
write.table(pls1_scores, file = "output/pls1_scores.txt", quote = FALSE, row.names = FALSE)


v <- selectNcomp(pls1, method = "randomization", nperm = 1000, alpha = 0.05, ncomp = 10, plot = TRUE)

# PLS model 1 permutation test of components

system.time(
  perm_stat <- t(sapply(1:100, function(i){
    y_perm <- sample(y)
    pls <- plsr(y_perm ~ x, ncomp = 10, scale = TRUE)
    # slope <- pls$Yloadings
    # as.vector(slope)
    r <- sapply(setNames(c(1:10), paste("Comp", c(1:10))), function(c){
      cor(pls$scores[,c], pls$Yscores[,c])
    })
    r
  })
))
slope <- pls1$Yloadings
perm_p <- sapply(setNames(c(1:10), paste("Comp", c(1:10))), function(c){
  s <- slope[c]
  perm_s <- perm_stat[,c]
  t <- (mean(perm_s)-s)/(sd(perm_s)/sqrt(length(y)))
  p <- 2*pt(t, df = length(y) - 1)
  p
})
perm_p

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


# Genes ranked by coefficients of 1st component
pls1_coef <- sort(pls1$coefficients[,1,1], decreasing = TRUE) # sort genes by pls1 coefficients
# write.table(entrezId2Name(names(pls1_coef)), file = "output/pls1_coef.txt", quote = FALSE, sep = "\t", 
#             col.names = FALSE, row.names = FALSE)

# Heatmap top genes for PLS1 of X
topgenes <- names(pls1_coef[1:50])
exprTopGenes <- scale(roi_expr[, topgenes])
colnames(exprTopGenes) <- entrezId2Name(colnames(exprTopGenes))
rownames(exprTopGenes) <- gsub("lh_", "", rois_lh)
exprTopGenes <- exprTopGenes[order(-pls1$scores[,1]), ] # order regions by PLS1 scores of X
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

# Functional enrichment with Reactome PA and GSEA
library('ReactomePA')
gsea <- gsePathway(pls1_coef, organism = "human", pAdjustMethod = "BH")
df <- as.data.frame(gsea)
df <- df[, c("Description", "p.adjust")]
df$p.adjust <- format(df$p.adjust, digits = 3, scientific = TRUE)
write.table(df, file = "output/GSEA_pls1_coef.txt", quote = FALSE, sep = "\t", row.names = FALSE)
options(stringsAsFactors = TRUE)
pdf("output/GSEA_pls1_coef.pdf", 10, 8)
emapplot(gsea, color = "pvalue")
dev.off()
