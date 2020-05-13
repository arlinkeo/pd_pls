# PLS model 1: single response variable
# Gene expression as predictors and T-score of cortical thickness as response

# Y: t-stats of comparing cortical thickness between PD and control,
tstats <- ct_test[grep("lh_", ct_test$Group), c("Group", "t")]
rownames(tstats) <- parcel_lh$id
tstats$t <- as.numeric(tstats$t)

# YexplVar <- function(x){ # explained variance of Y in PLS model
#   v <- 100*(drop(R2(x, estimate = "train", intercept = FALSE)$val))
#   v <- c(v[1],diff(v)) # reverse cumsum
#   names(v) <- paste("Comp", gsub(" comps", "", names(v)))
#   v
# }

# PLS
x <- roi_expr
y <- tstats$t
# pls-package
pls_model1 <- plsr(y ~ x, ncomp = 10, scale = TRUE, validation = "LOO") # x is scaled
pdf("output/pls1_CV_LOO.pdf", 6, 4)
plot(RMSEP(pls_model1), legendpos = "topright", main = "")
dev.off()
summary(pls_model1)
# explVarX <- explvar(pls_model1)
# explVarY <- YexplVar(pls_model1)
plsmodel1_scores_x <- scores(pls_model1)
tab <- cbind(ID = rownames(plsmodel1_scores_x), label = rois_lh, plsmodel1_scores_x)
write.table(tab, file = "output/plsmodel1_scores_x.txt", quote = FALSE, row.names = FALSE, sep = "\t")
plsmodel1_scores_y <- pls_model1$Yscores

# Plot PLS regression line
df <- data.frame(
  roi = parcel_lh$id, 
  name = gsub("lh_", "", parcel_lh$name), 
  comp1 = plsmodel1_scores_x[,1],
  tstats = y)
df$name[which(!((df$comp1 %in% range(df$comp1)) | (df$tstats %in% range(df$tstats))))] <- ""
r <- round(cor(df$comp1, y), digits = 2)
p <- ggplot(df, aes(comp1, tstats)) +
  geom_point(color = "blue", size = 2) + geom_smooth(method = "lm") +
  geom_text_repel(aes(label = name), size = 3, force = 2) +
  labs(x = bquote('PLS'~italic('component-1')~'score'), y = bquote(italic('t')*'-statistic of '~Delta~'CT')) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  ggtitle(bquote(italic('r')~'='~.(r))) +
  theme_classic()
pdf("output/scatterplot_plsmodel1_tstats.pdf", 4, 3)
p
dev.off()

# Genes ranked by projection to 1st component of X
gene_weights1 <- sort(pls_model1$projection[,1], decreasing = TRUE) # sort genes by pls1 projection
tab <- data.frame(Gene = entrezId2Name(names(gene_weights1)), Weight = gene_weights1)
write.table(tab, file = "output/plsmodel1_geneweights.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Functional enrichment with Reactome PA and GSEA
# abs.geneweights1 <- sort(abs(gene_weights1), decreasing = TRUE)
gsea1 <- gsePathway(geneweights1, organism = "human", pAdjustMethod = "BH")
df <- as.data.frame(gsea1)
df <- df[, c("Description", "p.adjust")]
df$p.adjust <- format(df$p.adjust, digits = 3, scientific = TRUE)
names(df) <- c("Pathway", "P-value")
write.table(df, file = "output/GSEA_plsmodel1_comp1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
options(stringsAsFactors = TRUE)
pdf("output/GSEA_plsmodel1_comp1.pdf", 9, 8)
emapplot(gsea1, color = "pvalue")
dev.off()
options(stringsAsFactors = FALSE)

# Get genesets of significant pathways and gene weights
pathways <- gsea1@geneSets[gsea1@result$ID]
names(pathways) <- gsea1@result$Description
pathways_weight <- lapply(pathways, function(g){
  gene_weights1[intersect(ahba.genes(), g)]
})

# Average expresssion of genes in each significant pathways
exprPathways <- sapply(pathways, function(g){
  g <- intersect(ahba.genes(), g)
  e <- roi_expr[, g]
  apply(e, 1, mean)
  # c <- gene_weights1[g] # weighted gene expression
  # e %*% c
})
exprPathways <- t(scale(exprPathways))
colnames(exprPathways) <- gsub("lh_", "", rois_lh)
col_order <- order(y)
exprPathways <- exprPathways[, col_order]

# Heatmap of pathways for PLS1 of X
pdf("output/heatmap_plsmodel1_comp1.pdf", 11.7, 11.8) #12.2, 17.5)
pls.heatmap(exprPathways, pathways_weight, y[col_order], 'gene weight', 't-statistic of Delta CT')
dev.off()

# Small version of heatmap
pdf("output/heatmap_plsmodel1_comp1_top30.pdf", 12, 5.8)
hm <- pls.heatmap(exprPathways[1:30, ], pathways_weight[1:30], y[col_order], 'gene weight', 't-statistic of Delta CT')
dev.off()