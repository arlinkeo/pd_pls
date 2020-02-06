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
y <- tscores$t # scaling makes no difference
# pls-package
pls_model1 <- plsr(y ~ x, ncomp = 10, scale = TRUE, validation = "LOO") # x is scaled
pdf("output/pls1_CV_LOO.pdf", 6, 4)
plot(RMSEP(pls_model1), legendpos = "topright", main = "")
dev.off()
summary(pls_model1)
# explVarX <- explvar(pls_model1)
# explVarY <- YexplVar(pls_model1)
pls1_scores <- scores(pls_model1)
tab <- cbind(ID = rownames(pls1_scores), label = rois_lh, pls1_scores)
write.table(tab, file = "output/pls1_scores.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# # plsdepot package
# pls_model1 <- plsreg1(x, y, comps = 10, crosval = TRUE)
# explVarX <- pls_model1$R2Xy
# explVarX <- explVarX[-nrow(explVarX),]
# explVarX <- apply(explVarX, 2, sum)/nrow(explVarX) # average explained variance of X across variables
# explVarY <- pls_model1$R2
# pls_scores <- pls_model1$x.scores
# # plot(pls_model1b$Q2[,1])

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

# Plot # PLS regression line
df <- data.frame(
  roi = parcel_lh$id, 
  name = gsub("lh_", "", parcel_lh$name), 
  pls1 = pls1_scores[,1],
  tscores = y)
df$name[which(!((df$pls1 %in% range(df$pls1)) | (df$tscores %in% range(df$tscores))))] <- ""
p <- ggplot(df, aes(pls1, tscores)) +
  geom_point(color = "blue", size = 2) + geom_smooth(method = "lm") +
  geom_text_repel(aes(label = name), size = 3, force = 2) +
  labs(x = "PLS component 1 score", y = expression("T-score of "~Delta~"CT")) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  ggtitle(paste("r =", round(cor(df$pls1, y), digits = 2))) +
  theme_classic()
pdf("output/PLS1_tscores.pdf", 4, 3)
p
dev.off()

# Genes ranked by coefficients of 1st component
pls1_coef <- sort(pls_model1$coefficients[,1,1], decreasing = TRUE) # sort genes by pls1 coefficients
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
gsea1.1 <- gsePathway(pls1_coef, organism = "human", pAdjustMethod = "BH")
df <- as.data.frame(gsea1.1)
df <- df[, c("Description", "p.adjust")]
df$p.adjust <- format(df$p.adjust, digits = 3, scientific = TRUE)
write.table(df, file = "output/GSEA_pls1_coef.txt", quote = FALSE, sep = "\t", row.names = FALSE)
options(stringsAsFactors = TRUE)
pdf("output/GSEA_pls1_coef.pdf", 9, 8)
emapplot(gsea1.1, color = "pvalue")
dev.off()

# Heatmap of pathways for PLS1 of X
pathways <- gsea1.1@geneSets[gsea1.1@result$ID]
names(pathways) <- gsea1.1@result$Description
pathways_avgcoef <- sapply(pathways, function(g){
  mean(pls1_coef[intersect(ahba.genes(), g)])
})

exprPathways <- sapply(pathways, function(g){
  g <- intersect(ahba.genes(), g)
  e <- roi_expr[, g]
  apply(e, 1, mean)
})
exprPathways <- t(scale(exprPathways))
colnames(exprPathways) <- gsub("lh_", "", rois_lh)
col_order <- order(-y)
row_order <- order(pathways_avgcoef)
pathways_avgcoef <- pathways_avgcoef[row_order]
exprPathways <- exprPathways[row_order, col_order]

# significant_regions <- ct_test[ct_test$BH < 0.05, c("Group", "t", "Mean Difference", "BH")]
# significant_regions <- significant_regions[grep("lh_", significant_regions$Group),]
# significant_regions$Group <- gsub("lh_|rh_", "", significant_regions$Group)
# cols <- as.numeric(colnames(exprPathways) %in% significant_regions$Group)
# cols <- cols*-5+8
ha_col <- HeatmapAnnotation('T-score of Delta CT' = anno_barplot(y[col_order]), #gp = gpar(fill = cols)),
                            height = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 8))

ha_row <- rowAnnotation('average PLS coefficient of genes' = anno_barplot(pathways_avgcoef), 
                        width = unit(3, "cm"), annotation_name_gp = gpar(fontsize = 8))

hm <- Heatmap(exprPathways, name = 'Z-Score\nexpression',
              # split = pathways_avgcoef > 0, 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 6),
              row_names_side = c("right"),
              row_title_rot = 0,
              column_names_side = c("top"),
              column_names_gp = gpar(fontsize = 6),
              column_names_rot = 45,
              width = unit(ncol(exprPathways)*.5, "lines"), 
              height = unit(nrow(exprPathways)*.5, "lines"),
              top_annotation = ha_col,
              right_annotation = ha_row
)
pdf("output/heatmap_pls1_pathways.pdf", 12, 17.5)
draw(hm, heatmap_legend_side = "left")
dev.off()

# Small version of heatmap
rows <- which(rownames(exprPathways) %in% c("Protein folding", "Apoptosis", "Regulation of RAS by GAPs", 
       "Cellular response to hypoxia", "Regulation of mitotic cell cycle",
       "Mitochondrial protein import", "Mitochondrial translation",
       "p53−Independent DNA Damage Response",
       "Stabilization of p53", "APC/C:Cdc20 mediated degradation of mitotic proteins",
       "Signaling by Interleukins", 
       "Circadian Clock", "Transcriptional Regulation by MECP2", 
       "Chromatin organization", "Nucleotide Excision Repair"))
rows <- c(rows, grep("DNA damage|DNA Damage|SUMO", rownames(exprPathways)))
rows <- unique(rows)
rows <- sort(rows)
exprPathways <- exprPathways[rows, ]
ha_row <- rowAnnotation('average PLS coefficient of genes' = anno_barplot(pathways_avgcoef[rows], annotation_name_side = "left"), 
                        width = unit(3, "cm"), annotation_name_gp = gpar(fontsize = 8))
hm <- Heatmap(exprPathways, name = 'Z-Score\nexpression',
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 6),
              row_names_side = c("right"),
              row_title_rot = 0,
              column_names_side = c("top"),
              column_names_gp = gpar(fontsize = 6),
              column_names_rot = 45,
              width = unit(ncol(exprPathways)*.5, "lines"), 
              height = unit(nrow(exprPathways)*.5, "lines"),
              top_annotation = ha_col,
              right_annotation = ha_row
)
pdf("output/heatmap_pls1_pathways_reduced.pdf", 8, 4.5)
hm
dev.off()


# # Cell-type enrichment
# markerlist <- readRDS("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn/output/markerlist.rds")
# celltype_enrichment <- hypertest.or.table(list(names(pls1_coef[1:200])), markerlist, 20017) # 3d-array: cell-types x measures x deg-type
# apply(celltype_enrichment, 3, function(x){
#   x[x[, "bh"] < 0.05, ]
# })

# # Disease enrichment
# disease_list <- readRDS("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn/output/disease_genes.rds")
# celltype_enrichment <- hypertest.or.table(list(names(pls1_coef[1:50])), disease_list, 20017) # 3d-array: cell-types x measures x deg-type
# apply(celltype_enrichment, 3, function(x){
#   x[x[, "bh"] < 0.05, ]
# })
