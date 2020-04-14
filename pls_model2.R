# Correlation coefficient of clinical scores

# Relation between ct and cs for at most 123 PD patients
coef_eq1 <- lapply(clinicalScores, function(cs){ # For each score across patients
  sapply(ct_lh[patientIDs, ], function(ct){ # For each region with ct across patients
    res <- lm(cs ~ ct)
    res  <- summary(res)
    res$coefficients[2, c(1,3,4)]# coefficient + p-value
  })
})
coef_eq1 <- simplify2array(coef_eq1) # 3D-array: measures (estimate, p-value) x regions x clinical features
dimnames(coef_eq1)[[2]] <- gsub("lh_", "", dimnames(coef_eq1)[[2]])
dimnames(coef_eq1)[[3]] <- gsub("_", " ", dimnames(coef_eq1)[[3]])
df <- melt(coef_eq1["Pr(>|t|)", , ])
df$value <- p.adjust(df$value, method = "BH") # BH-correct for regions and clinical features
bh <- dcast(df, Var1 ~ Var2)[, -1]
coef_eq1 <- abind(coef_eq1, 'BH' = bh, along = 1) # BH added to coef
range(coef_eq1["t value", , ])
# sum(coef_eq1["BH", , ] < 0.05)

# Heatmaps for coefficients and BH-corrected P-values
t <- t(coef_eq1["t value", , ])
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
                 if(t(coef_eq1["BH", i, j]) < 0.05) {
                   grid.text("*", x = x, y = y)
                 } else {
                   grid.text("", x = x, y = y)
                 }
               }
)
pdf("output/heatmap_clinical_scores.pdf", 7.8, 3.8)
hm1
dev.off()

YexplVar <- function(x){ # explained variance of Y in PLS model
  v <- R2(x, estimate = "train", intercept = FALSE)$val[1,,]*100
  v <- apply(v, 2, mean)
  v <- c(v[1],diff(v)) # reverse cumsum
  names(v) <- paste("Comp", gsub(" comps", "", names(v)))
  v
}

# PLS with al features
x <- roi_expr
y <- coef_eq1["t value", , ]
pls_model2 <- plsr(y ~ x, ncomp = 10, scale = TRUE, validation = "LOO")
summary(pls_model2)
explVarX <- explvar(pls_model2)
explVarY <- YexplVar(pls_model2)
pls2_scores_x <- pls_model2$scores#scores(pls_model2)
tab <- cbind(ID = c(1:34), label = rois_lh, pls2_scores_x)
write.table(tab, file = "output/plsmodel2_scores_x.txt", quote = FALSE, row.names = FALSE, sep = "\t")
pls2_scores_y <- pls_model2$Yscores
tab <- cbind(ID = c(1:34), label = rois_lh, pls2_scores_y)
write.table(tab, file = "output/plsmodel2_scores_y.txt", quote = FALSE, row.names = FALSE, sep = "\t")

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
p <- lapply(c(1:3), function(i){
  df <- data.frame(
    roi = parcel_lh$id, 
    name = gsub("lh_", "", parcel_lh$name), 
    p.pls = pls_model2$scores[,i], # predictor PLS1 scores
    r.pls = pls_model2$Yscores[,i] # response PLS1 scores
  )
  df$name[which(!((df$p.pls %in% range(df$p.pls)) | (df$r.pls %in% range(df$r.pls))))] <- ""
  ggplot(df, aes(p.pls, r.pls)) +
    geom_point(color = "blue") + geom_smooth(method = "lm") +
    geom_text_repel(aes(label = name), size = 3) +
    labs(x = bquote('PLS'~italic('component-'*.(i))~"of gene expression ("*.(round(explVarX[i]))*"%)"), 
         y = bquote('PLS'~italic('component-'*.(i))~"of"~beta*"'s ("*.(round(explVarY[i]))*'%)')) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", color = "gray") +
    ggtitle(paste("r =", round(cor(df$p.pls, df$r.pls), digits = 2))) +
    theme_classic()
})
pdf("output/pls_clinical_scores.pdf", 12, 3)
ggarrange(plotlist = p, nrow = 1)
dev.off()

# Genes ranked by projection of 1st, 2nd, and 3th component of X
gene_weights2 <- sapply(colnames(pls_model2$projection)[1:3], function(n){
  r <- pls_model2$projection[, n]
  r <- sort(r, decreasing = TRUE)
  tab <- data.frame(Gene = entrezId2Name(names(r)), Weight = r)
  write.table(tab, file = paste0("output/plsmodel2_geneweights", gsub("Comp ", "comp", n), ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  r
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
gsea2 <- sapply(names(gene_weights2)[1:3], function(i){
  gsea <- gsePathway(gene_weights2[[i]], organism = "human", pAdjustMethod = "BH")
  df <- as.data.frame(gsea)
  df <- df[, c("Description", "p.adjust")]
  df$p.adjust <- format(df$p.adjust, digits = 3, scientific = TRUE)
  write.table(df, file = paste0("output/GSEA_plsmodel2_", gsub("Comp ", "comp", i), ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  gsea
}, simplify = FALSE)

# Emapplots
options(stringsAsFactors = TRUE)
pdf("output/GSEA_plsmodel2.pdf", 9, 8)
emapplot(gsea2$`Comp 1`, color = "pvalue")
emapplot(gsea2$`Comp 2`, color = "pvalue")
emapplot(gsea2$`Comp 3`, color = "pvalue")
dev.off()
options(stringsAsFactors = FALSE)

# Overlap pathways
p <- lapply(gsea2, function(gsea){
  gsea@result$ID
})
lengths(p)
overlapping_pathways <- Reduce(intersect, p[[1:2]])
length(overlapping_pathways)
sapply(p, function(x) length(intersect(x, gsea1.1@result$ID))) # overlap with PLS model-1

# Heatmap of pathways for components of PLS model-2
lapply(names(gsea2)[-which(lengths(p)==0)], function(i){
  
  # Get genesets of significant pathways abnd average the PLS coefficients of genes
  pathways <- gsea2[[i]]@geneSets[gsea2[[i]]@result$ID]
  names(pathways) <- gsea2[[i]]@result$Description
  pathways_avgweight <- sapply(pathways, function(g){
    mean(gene_weights2[[i]][intersect(ahba.genes(), g)]) # average PLS coefficient of genes for each significant pathway
  })
  
  # Average expresssion of genes in each significant pathways
  exprPathways <- sapply(pathways, function(g){
    g <- intersect(ahba.genes(), g)
    e <- roi_expr[, g]
    apply(e, 1, mean)
  })
  exprPathways <- t(scale(exprPathways))
  colnames(exprPathways) <- gsub("lh_", "", rois_lh)
  col_order <- order(pls2_scores_y[, i])
  row_order <- order(pathways_avgweight)
  pathways_avgweight <- pathways_avgweight[row_order]
  exprPathways <- exprPathways[row_order, col_order]
  
  # Heatmap of pathways for PLS component i
  hm <- pls_heatmap(exprPathways, pathways_avgweight, pls2_scores_y[col_order, i], 
                    'average gene weight', paste0('PLS component-', i, ' score of Betas'))
  pdf(paste0("output/heatmap_plsmodel2_", gsub("Comp ", "comp", i), ".pdf"), 12.2, length(pathways)/10+2.2)
  draw(hm, heatmap_legend_side = "left")
  dev.off()
  
  # Heatmap of top10 +ve and -ve correlated pathways
  idx <- c(c(1:10), c((length(pathways_avgweight)-9):length(pathways_avgweight)))
  hm <- pls_heatmap(exprPathways[idx, ], pathways_avgweight[idx], pls2_scores_y[col_order, i], 
                    'average gene weight', paste0('PLS component-', i, ' score of Betas'))
  pdf(paste0("output/heatmap_plsmodel2_", gsub("Comp ", "comp", i), "_top10.pdf"), 10, 4.2)
  draw(hm, heatmap_legend_side = "left")
  dev.off()
  
})