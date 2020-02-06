# Correlation coefficient of clinical scores

# Relation between ct and cs for at most 123 PD patients
coef <- lapply(clinicalScores, function(cs){ # For each score across patients
  sapply(ct_lh[patientIDs, ], function(ct){ # For each region with ct across patients
    res <- lm(cs ~ ct)
    res  <- summary(res)
    res$coefficients[2, c(1,3,4)]# coefficient + p-value
  })
})
coef <- simplify2array(coef) # 3D-array: measures (estimate, p-value) x regions x clinical features
dimnames(coef)[[2]] <- gsub("lh_", "", dimnames(coef)[[2]])
dimnames(coef)[[3]] <- gsub("_", " ", dimnames(coef)[[3]])
df <- melt(coef["Pr(>|t|)", , ])
df$value <- p.adjust(df$value, method = "BH") # BH-correct for regions and clinical features
bh <- dcast(df, Var1 ~ Var2)[, -1]
coef <- abind(coef, 'BH' = bh, along = 1) # BH added to coef
range(coef["t value", , ])
  sum(coef["BH", , ] < 0.05)


# Heatmaps for coefficients and BH-corrected P-values
t <- t(coef["t value", , ])
q1 <- max(quantile(abs(t), 0.9))
col_fun <- colorRamp2(c(-q1, 0, q1), c("blue", "#EEEEEE", "red"))# Heat colors centered around 0
hm1 <- Heatmap(t, name = "T-score\nof slope",
              col = col_fun,
              cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_side = c("left"),
              row_title_rot = 0,
              column_names_side = c("top"),
              column_names_gp = gpar(fontsize = 10),
              width = unit(ncol(t)*.8, "lines"), 
              height = unit(nrow(t)*.8, "lines"),
              heatmap_legend_param = list(title_position = "topleft"),
              cell_fun = function(i, j, x, y, width, height, fill) {
                if(t(coef["BH", i, j]) < 0.05) {
                  grid.text("*", x = x, y = y)
                } else {
                  grid.text("", x = x, y = y)
                }
              }
)
pdf("output/heatmap_clinical_scores.pdf", 7.8, 3.8)
hm1
dev.off()

# # PLS per feature
# l <- lapply(colnames(coef), function(n){
#   x <- coef[, n]
#   pls <- plsr(x ~ roi_expr, scale = TRUE)
#   df <- data.frame(
#     roi = parcel_lh$id, 
#     name = gsub("lh_", "", parcel_lh$name), 
#     pls1 = pls$scores[,1],
#     pls2 = x)
#   p <- ggplot(df, aes(pls1, pls2)) +
#     geom_point(color = "blue") + geom_smooth(method = "lm") +
#     geom_text_repel(aes(label = name), size = 3) +
#     labs(x = "PLS1", y = n) +
#     geom_hline(yintercept=0, linetype="dashed", color = "gray") +
#     geom_vline(xintercept=0, linetype="dashed", color = "gray") +
#     ggtitle(paste(n, ": r =", round(cor(df$pls1, x), digits = 2))) +
#     theme_classic()
#   p
# })
# pdf("output/pls_clinical_scores.pdf", 32, 32)
# ggarrange(plotlist = l, ncol = 6, nrow = 7)
# dev.off()

YexplVar <- function(x){ # explained variance of Y in PLS model
  v <- R2(x, estimate = "train", intercept = FALSE)$val[1,,]*100
  v <- apply(v, 2, mean)
  v <- c(v[1],diff(v)) # reverse cumsum
  names(v) <- paste("Comp", gsub(" comps", "", names(v)))
  v
}

# PLS with al features
x <- roi_expr
# write.table(x, file = "output/x.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
y <- coef["t value", , ]
# write.table(y, file = "output/y.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
pls_model2 <- plsr(y ~ x, ncomp = 10, scale = TRUE, validation = "LOO")
summary(pls_model2)
# explVar <- R2(pls_model2)
explVarX <- explvar(pls_model2)
explVarY <- YexplVar(pls_model2)
pls2_scores <- scores(pls_model2)
tab <- cbind(ID = c(1:34), label = rois_lh, pls2_scores)
write.table(tab, file = "output/pls2_scores.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# # plsdepot package
# pls_model2b <- plsreg2(x, y, comps = 10)
# pls_model2b$expvar

# Plot explained variance
df <- data.frame(Comp = gsub("Comp ", "", names(explVarX)), explVarX = explVarX, cumexplVarX = cumsum(explVarX))
df <- melt(df, measure.vars = c("explVarX", "cumexplVarX"))
df$Comp <- factor(df$Comp, levels = unique(df$Comp))
df$variable <- factor(df$variable, levels = unique(df$variable), labels = c("Explained variance", "Cumulative explained variance"))
p <- ggplot(df, aes(Comp, value, group = variable)) +
  geom_point(aes(color = variable)) +
  geom_line(aes(color = variable)) +
  labs(x = "Number of components", y = "Explained variance") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "top")
pdf("output/pls_model2_explained_variance.pdf", 4, 3)
p
dev.off()

# # PLS model 2 permutation test of components
# system.time(
# perm_stat <- t(sapply(1:1000, function(i){
#   order <- sample(1:nrow(y))
#   y_perm <- y[order,]
#   pls <- plsr(y_perm ~ x, ncomp = 10, scale = TRUE)
#   YexplVar(pls)
# }))
# )
# perm_p <- sapply(paste("Comp", c(1:10)), function(c){
#   sum(perm_stat[,c] >= explVarY[c])/1000
# })
# perm_p

# Plot
p <- lapply(c(1:4), function(i){
  df <- data.frame(
    roi = parcel_lh$id, 
    name = gsub("lh_", "", parcel_lh$name), 
    p.pls = pls_model2$scores[,i], # predictor PLS1 scores
    r.pls = pls_model2$Yscores[,i] # response PLS1 scores
  )
  ggplot(df, aes(p.pls, r.pls)) +
    geom_point(color = "blue") + geom_smooth(method = "lm") +
    # geom_text_repel(aes(label = name), size = 3) +
    labs(x = paste0("PLS", i, " Gene expression (", round(explVarX[i]), "%)"), 
         y = paste0("PLS", i, " Beta (", round(explVarY[i]), "%)")) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", color = "gray") +
    ggtitle(paste("r =", round(cor(df$p.pls, df$r.pls), digits = 2))) +
    theme_classic()
})
pdf("output/pls_clinical_scores.pdf", 16, 3)
ggarrange(plotlist = p, nrow = 1)
dev.off()

# Genes ranked by coefficients of components
pls2_coef <- sapply(dimnames(pls_model2$coefficients)[[3]], function(i){
  pls2_coef <- sort(pls_model2$coefficients[,1,i], decreasing = TRUE) # sort genes by pls1 coefficients
}, simplify = FALSE)

# # Heatmap top rank genes
# hm <- apply(pls2_coef, 2, function(c){
#   exprTopGenes <- scale(x[, c[1:50]])
#   colnames(exprTopGenes) <- entrezId2Name(colnames(exprTopGenes))
#   rownames(exprTopGenes) <- gsub("lh_", "", rois_lh)
#   # row_order <- order(-y)
#   # exprTopGenes <- exprTopGenes[row_order, ] # order regions by T-score of delta CT
#   # ha <- rowAnnotation('T-score' = anno_barplot(y[row_order]), width = unit(4, "cm"))
#   hm <- Heatmap(exprTopGenes, name = 'Z-Score\nexpression',
#                 cluster_rows = FALSE,
#                 cluster_columns = FALSE,
#                 row_names_gp = gpar(fontsize = 10),
#                 row_names_side = c("left"),
#                 row_title_rot = 0,
#                 column_names_side = c("top"),
#                 column_names_gp = gpar(fontsize = 10, fontface = "italic"),
#                 column_names_rot = 45,
#                 width = unit(ncol(exprTopGenes)*.8, "lines"), 
#                 height = unit(nrow(exprTopGenes)*.8, "lines")
#   )
# })
# pdf("output/heatmap_pls2_topgenes.pdf", 28, 6.7)
# draw(Reduce('+', hm), gap = unit(2, "cm")) # Reduce('+', hm)
# dev.off()

# Functional enrichment with Reactome PA and GSEA
gsea2 <- sapply(names(pls2_coef)[1:4], function(i){
  gsea <- gsePathway(pls2_coef[[i]], organism = "human", pAdjustMethod = "BH")
  df <- as.data.frame(gsea)
  df <- df[, c("Description", "p.adjust")]
  df$p.adjust <- format(df$p.adjust, digits = 3, scientific = TRUE)
  write.table(df, file = paste0("output/GSEA_pls2_comp", gsub(" comps", "", i), ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  gsea
}, simplify = FALSE)

# Overlap pathways
p <- lapply(gsea2, function(gsea){
  gsea@result$ID
})
lengths(p)
overlapping_pathways <- Reduce(intersect, p)
length(overlapping_pathways)
sapply(p, function(x) length(intersect(x, gsea1.1@result$ID)))
# 
# # Heatmap of pathways for PLS2 of X
# pathways <- gsea2[[3]]@geneSets[overlapping_pathways]
# names(pathways) <- gsea2[[3]]@result[names(pathways), "Description"]
# exprPathways <- sapply(pathways, function(g){
#   g <- intersect(ahba.genes(), g)
#   e <- roi_expr[, g]
#   apply(e, 1, mean)
# })
# exprPathways <- t(scale(exprPathways))
# colnames(exprPathways) <- gsub("lh_", "", rois_lh)
# 
# # pathways_avgcoef <- lapply(pls2_coef, function(c){
# # 
# # })
# # 
# #   sapply(overlapping_pathways, function(g){
# #   mean(pls1_coef[intersect(ahba.genes(), g)])
# # })
# # 
# # 
# # col_order <- order(-y)
# # row_order <- order(pathways_avgcoef)
# # pathways_avgcoef <- pathways_avgcoef[row_order]
# exprPathways <- exprPathways[, col_order]
# 
# # significant_regions <- ct_test[ct_test$BH < 0.05, c("Group", "t", "Mean Difference", "BH")]
# # significant_regions <- significant_regions[grep("lh_", significant_regions$Group),]
# # significant_regions$Group <- gsub("lh_|rh_", "", significant_regions$Group)
# # cols <- as.numeric(colnames(exprPathways) %in% significant_regions$Group)
# # cols <- cols*-5+8
# ha_col <- HeatmapAnnotation('Beta' = anno_barplot(y[col_order]), #gp = gpar(fill = cols)),
#                             height = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 8))
# 
# # ha_row <- rowAnnotation('average PLS coefficient of genes' = anno_barplot(pathways_avgcoef), 
# #                         width = unit(3, "cm"), annotation_name_gp = gpar(fontsize = 8))
# 
# hm <- Heatmap(exprPathways, name = 'Z-Score\nexpression',
#               # split = pathways_avgcoef > 0, 
#               cluster_rows = FALSE,
#               cluster_columns = FALSE,
#               row_names_gp = gpar(fontsize = 6),
#               row_names_side = c("right"),
#               row_title_rot = 0,
#               column_names_side = c("top"),
#               column_names_gp = gpar(fontsize = 6),
#               column_names_rot = 45,
#               width = unit(ncol(exprPathways)*.5, "lines"), 
#               height = unit(nrow(exprPathways)*.5, "lines")
#               # top_annotation = ha_col,
#               # right_annotation = ha_row
# )
# pdf("output/heatmap_pls2_pathways.pdf", 12, 17.5)
# draw(hm, heatmap_legend_side = "left")
# dev.off()
# 
# # Small version of heatmap
# rows <- which(rownames(exprPathways) %in% c("Protein folding", "Apoptosis", "Regulation of RAS by GAPs", 
#                                             "Cellular response to hypoxia", "Regulation of mitotic cell cycle",
#                                             "Mitochondrial protein import", "Mitochondrial translation",
#                                             "p53−Independent DNA Damage Response",
#                                             "Stabilization of p53", "APC/C:Cdc20 mediated degradation of mitotic proteins",
#                                             "Signaling by Interleukins", 
#                                             "Circadian Clock", "Transcriptional Regulation by MECP2", 
#                                             "Chromatin organization", "Nucleotide Excision Repair"))
# rows <- c(rows, grep("DNA damage|DNA Damage|SUMO|mitochondrial", rownames(exprPathways)))
# rows <- unique(rows)
# rows <- sort(rows)
# exprPathways <- exprPathways[rows, ]
# # ha_row <- rowAnnotation('average PLS coefficient of genes' = anno_barplot(pathways_avgcoef[rows], annotation_name_side = "left"), 
#                         # width = unit(3, "cm"), annotation_name_gp = gpar(fontsize = 8))
# hm <- Heatmap(exprPathways, name = 'Z-Score\nexpression',
#               cluster_rows = FALSE,
#               cluster_columns = FALSE,
#               row_names_gp = gpar(fontsize = 6),
#               row_names_side = c("right"),
#               row_title_rot = 0,
#               column_names_side = c("top"),
#               column_names_gp = gpar(fontsize = 6),
#               column_names_rot = 45,
#               width = unit(ncol(exprPathways)*.5, "lines"), 
#               height = unit(nrow(exprPathways)*.5, "lines"),
#               top_annotation = ha_col#,
#               # right_annotation = ha_row
# )
# pdf("output/heatmap_pls2_pathways_reduced.pdf", 8, 4.5)
# hm
# dev.off()
