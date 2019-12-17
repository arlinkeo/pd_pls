# PLS model 1: single response variable
# Gene expression as predictors and T-score of cortical thickness as response

# Y: t-stats of comparing cortical thickness between PD and control,
tscores <- ct_test[grep("lh_", ct_test$Group), c("Group", "t")]
rownames(tscores) <- parcel_lh$id
tscores$t <- as.numeric(tscores$t)

YexplVar <- function(x){ # explained variance of Y in PLS model
  v <- 100*(drop(R2(x, estimate = "train", intercept = FALSE)$val))
  # v <- c(v[1],diff(v)) # reverse cumsum
  names(v) <- paste("Comp", gsub(" comps", "", names(v)))
  v
}

# PLS
x <- roi_expr
y <- tscores$t
pls_model1 <- plsr(y ~ x, ncomp = 10, scale = TRUE, validation = "LOO")
pdf("output/pls1_CV_LOO.pdf", 6, 4)
plot(RMSEP(pls_model1), legendpos = "topright", main = "")
dev.off()
summary(pls_model1)
explVarX <- explvar(pls_model1)
explVarY <- YexplVar(pls_model1)
# plot(pls_model1, "loadings", comps = 1:2, legendpos = "topleft", xlab = "gene")
# paste(gsub("Comp ", "PLS", names(explVarX)), paste0(round(explVarX, digits = 2), "%"), sep = ": ", collapse = "; ")
pls1_scores <- scores(pls_model1)
pls1_scores <- cbind(ID = rownames(pls1_scores), label = rois_lh, pls1_scores)
write.table(pls1_scores, file = "output/pls1_scores.txt", quote = FALSE, row.names = FALSE)

# v <- selectNcomp(pls_model1, method = "randomization", nperm = 1000, alpha = 0.05, ncomp = 10, plot = TRUE)

# PLS model 1 permutation test of components
system.time(
  perm_stat <- t(sapply(1:1000, function(i){
    y_perm <- sample(y)
    pls <- plsr(y_perm ~ x, ncomp = 10, scale = TRUE)
    # YexplVar(pls)
    r <- sapply(paste("Comp", c(1:10)), function(c){
      cor(scores(pls)[,c], Yscores(pls)[,c])
    })
  })
))
perm_p <- sapply(paste("Comp", c(1:10)), function(c){
  sum(perm_stat[,c] >= explVarY[c])/1000
  # r <- cor(scores(pls_model1)[,c], Yscores(pls_model1)[,c])
  # sum(perm_stat[,c] >= r)/1000
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
