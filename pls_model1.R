# PLS model 1: single response variable
# Gene expression as predictors and T-score of cortical thickness as response

# Y: t-stats of comparing cortical thickness between PD and control,
tscores <- ct_test[grep("lh_", ct_test$Group), c("Group", "t")]
rownames(tscores) <- parcel_lh$id
tscores$t <- as.numeric(tscores$t)

# YexplVar <- function(x){ # explained variance of Y in PLS model
#   v <- 100*(drop(R2(x, estimate = "train", intercept = FALSE)$val))
#   v <- c(v[1],diff(v)) # reverse cumsum
#   names(v) <- paste("Comp", gsub(" comps", "", names(v)))
#   v
# }

# PLS
x <- roi_expr
y <- tscores$t
# pls-package
pls_model1 <- plsr(y ~ x, ncomp = 10, scale = TRUE, validation = "LOO") # x is scaled
pdf("output/pls1_CV_LOO.pdf", 6, 4)
plot(RMSEP(pls_model1), legendpos = "topright", main = "")
dev.off()
summary(pls_model1)
# explVarX <- explvar(pls_model1)
# explVarY <- YexplVar(pls_model1)
pls1_scores_x <- scores(pls_model1)
tab <- cbind(ID = rownames(pls1_scores_x), label = rois_lh, pls1_scores_x)
write.table(tab, file = "output/pls1_scores.txt", quote = FALSE, row.names = FALSE, sep = "\t")
pls1_scores_y <- pls_model1$Yscores

# system.time(
# # PLS model 1 permutation test of components
# perm_stat <- t(sapply(1:1000, function(i){
#   y_perm <- sample(y)
#   # pls <- plsr(y_perm ~ x, ncomp = 10, scale = TRUE)
#   # YexplVar(pls)
#   pls <- plsreg1(x, y_perm, comps = 10, crosval = TRUE)
#   pls$R2
#   # r <- sapply(paste("Comp", c(1:10)), function(c){
#   #   cor(scores(pls)[,c], Yscores(pls)[,c])
#   # })
# })))
# 
# perm_p <- sapply(paste0("t", c(1:10)), function(c){
#   sum(perm_stat[,c] >= explVarY[c])/1000
#   # r <- cor(scores(pls_model1)[,c], Yscores(pls_model1)[,c])
#   # sum(perm_stat[,c] >= r)/1000
# })
# perm_p

# Plot PLS regression line
df <- data.frame(
  roi = parcel_lh$id, 
  name = gsub("lh_", "", parcel_lh$name), 
  pls1 = pls1_scores_x[,1],
  tscores = y)
df$name[which(!((df$pls1 %in% range(df$pls1)) | (df$tscores %in% range(df$tscores))))] <- ""
r <- round(cor(df$pls1, y), digits = 2)
p <- ggplot(df, aes(pls1, tscores)) +
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

# Genes ranked by projection of 1st component of X
gene_weights1 <- sort(pls_model1$projection[,1], decreasing = TRUE) # sort genes by pls1 projection
tab <- data.frame(Gene = entrezId2Name(names(gene_weights1)), Weight = gene_weights1)
write.table(tab, file = "output/plsmodel1_geneweights.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# pls1_coef <- sort(pls_model1$coefficients[,1,1], decreasing = TRUE) # sort genes by pls1 coefficients
# write.table(entrezId2Name(names(pls1_coef)), file = "output/pls1_coef.txt", quote = FALSE, sep = "\t", 
#             col.names = FALSE, row.names = FALSE)

# # Heatmap top genes for PLS1 of X
# topgenes <- names(pls1_coef[1:50])
# exprTopGenes <- scale(roi_expr[, topgenes])
# colnames(exprTopGenes) <- entrezId2Name(colnames(exprTopGenes))
# rownames(exprTopGenes) <- gsub("lh_", "", rois_lh)
# row_order <- order(-y)
# exprTopGenes <- exprTopGenes[row_order, ] # order regions by T-score of delta CT
# ha <- rowAnnotation('T-score' = anno_barplot(y[row_order]), width = unit(4, "cm"))
# hm <- Heatmap(exprTopGenes, name = 'Z-Score\nexpression',
#               cluster_rows = FALSE,
#               cluster_columns = FALSE,
#               row_names_gp = gpar(fontsize = 10),
#               row_names_side = c("left"),
#               row_title_rot = 0,
#               column_names_side = c("top"),
#               column_names_gp = gpar(fontsize = 10, fontface = "italic"),
#               column_names_rot = 45,
#               width = unit(ncol(exprTopGenes)*.8, "lines"), 
#               height = unit(nrow(exprTopGenes)*.8, "lines")
# )
# pdf("output/heatmap_pls1_topgenes.pdf", 12, 6.7)
# hm + ha
# dev.off()

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

# # Cell-type enrichment
# markerlist <- readRDS("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn/output/markerlist.rds")
# celltype_enrichment <- hypertest.or.table(list(names(gene_weights1[1:200])), markerlist, 20017) # 3d-array: cell-types x measures x deg-type
# apply(celltype_enrichment, 3, function(x){
#   x[x[, "bh"] < 0.05, ]
# })
# 
# # Disease enrichment
# disease_list <- readRDS("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn/output/disease_genes.rds")
# celltype_enrichment <- hypertest.or.table(list(names(gene_weights1[1:200])), disease_list, 20017) # 3d-array: cell-types x measures x deg-type
# apply(celltype_enrichment, 3, function(x){
#   x[x[, "bh"] < 0.05, ]
# })
