# PLS model 2: multi-response variables
# Gene expression as predictors and the relatinship between cortical thickness and clinical scores as response

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
plsmodel2_scores_x <- pls_model2$scores#scores(pls_model2)
tab <- cbind(ID = c(1:34), label = rois_lh, plsmodel2_scores_x)
write.table(tab, file = "output/plsmodel2_scores_x.txt", quote = FALSE, row.names = FALSE, sep = "\t")
plsmodel2_scores_y <- pls_model2$Yscores
tab <- cbind(ID = c(1:34), label = rois_lh, plsmodel2_scores_y)
write.table(tab, file = "output/plsmodel2_scores_y.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Plot explained variance
df <- data.frame(Comp = gsub("Comp ", "", names(explVarX)), explVarX = explVarX, cumexplVarX = cumsum(explVarX))
df <- melt(df, measure.vars = c("explVarX", "cumexplVarX"))
df$Comp <- factor(df$Comp, levels = unique(df$Comp))
df$variable <- factor(df$variable, levels = unique(df$variable), labels = c("Explained variance", "Cumulative explained variance"))
p1 <- ggplot(df, aes(Comp, value, group = variable)) +
  geom_point(aes(color = variable)) +
  geom_line(aes(color = variable)) +
  labs(x = "Number of components", y = "Explained variance") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "top")
df <- data.frame(Comp = gsub("Comp ", "", names(explVarY)), explVarY = explVarY, cumexplVarY = cumsum(explVarY))
df <- melt(df, measure.vars = c("explVarY", "cumexplVarY"))
df$Comp <- factor(df$Comp, levels = unique(df$Comp))
df$variable <- factor(df$variable, levels = unique(df$variable), labels = c("Explained variance", "Cumulative explained variance"))
p2 <- ggplot(df, aes(Comp, value, group = variable)) +
  geom_point(aes(color = variable)) +
  geom_line(aes(color = variable)) +
  labs(x = "Number of components", y = "Explained variance") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "top")
pdf("output/pls_model2_explained_variance.pdf", 4, 3)
p1
p2
dev.off()

# Scatter plot
p <- lapply(c(1:3), function(i){
  df <- data.frame(
    roi = parcel_lh$id, 
    name = gsub("lh_", "", parcel_lh$name), 
    p.pls = pls_model2$scores[,i], # predictor PLS1 scores
    r.pls = pls_model2$Yscores[,i] # response PLS1 scores
  )
  df$name[which(!((df$p.pls %in% range(df$p.pls)) | (df$r.pls %in% range(df$r.pls))))] <- ""
  r = round(cor(df$p.pls, df$r.pls), digits = 2)
  ggplot(df, aes(p.pls, r.pls)) +
    geom_point(color = "blue") + geom_smooth(method = "lm") +
    geom_text_repel(aes(label = name), size = 3) +
    labs(x = bquote('PLS'~italic('component-'*.(i))~'of gene expression ('*.(round(explVarX[i]))*'%)'), 
         y = bquote(atop('PLS'~italic('component-'*.(i))~'of relationship', 'between CT and clinical scores'~'('*.(round(explVarY[i]))*'%)'))) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", color = "gray") +
    ggtitle(bquote(italic('r')~'='~.(r))) +
    theme_classic()
})
pdf("output/scatterplot_plsmodel2.pdf", 12, 3)
ggarrange(plotlist = p, nrow = 1)
dev.off()

# Genes ranked by projection of 1st, 2nd, and 3th component of X
gene_weights2 <- sapply(colnames(pls_model2$projection)[1:3], function(n){
  r <- pls_model2$projection[, n]
  r <- sort(r, decreasing = TRUE)
  tab <- data.frame(Gene = entrezId2Name(names(r)), Weight = r)
  write.table(tab, file = paste0("output/plsmodel2_geneweights-", gsub("Comp ", "comp", n), ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  r
}, simplify = FALSE)

# Functional enrichment with Reactome PA and GSEA
gsea2 <- sapply(names(gene_weights2)[1:3], function(i){
  gsea <- gsePathway(gene_weights2[[i]], organism = "human", pAdjustMethod = "BH")
  gsea
}, simplify = FALSE)

# Write tables
lapply(names(gsea2), function(i){
  df <- as.data.frame(gsea2[[i]])
  df <- df[, c("Description", "p.adjust")]
  df$p.adjust <- format(df$p.adjust, digits = 3, scientific = TRUE)
  names(df) <- c("Pathway", "P-value")
  write.table(df, file = paste0("output/GSEA_plsmodel2_", gsub("Comp ", "comp", i), ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
})

# Emapplots
options(stringsAsFactors = TRUE)
pdf("output/GSEA_plsmodel2.pdf", 9, 8)
emapplot(gsea2$`Comp 1`, color = "pvalue")
emapplot(gsea2$`Comp 2`, color = "pvalue")
dev.off()
options(stringsAsFactors = FALSE)

# Heatmap of pathways for components of PLS model-2
lapply(names(gsea2)[1:2], function(i){
  
  # Get genesets of significant pathways and  gene weights
  pathways <- gsea2[[i]]@geneSets[gsea2[[i]]@result$ID]
  names(pathways) <- gsea2[[i]]@result$Description
  pathways_weight <- lapply(pathways, function(g){
    gene_weights2[[i]][intersect(ahba.genes(), g)]
  })

  # Average expresssion of genes in each significant pathways
  exprPathways <- sapply(pathways, function(g){
    g <- intersect(ahba.genes(), g)
    e <- roi_expr[, g]
    apply(e, 1, mean)
  })
  exprPathways <- t(scale(exprPathways))
  colnames(exprPathways) <- gsub("lh_", "", rois_lh)
  col_order <- order(plsmodel2_scores_y[, i])
  exprPathways <- exprPathways[, col_order]
  
  # Heatmap of pathways for PLS component i
  # Split table if too many rows
  nrow <- nrow(exprPathways)
  maxrow <- 80
  if (nrow > maxrow){
    pdf(paste0("output/heatmap_plsmodel2_", gsub("Comp ", "comp", i), ".pdf"), 13, 11)
    lapply(1:ceiling(nrow/maxrow), function(j){
      print(paste(i, j))
      n <- min(j*maxrow, nrow)
      rows <- ((j*maxrow-maxrow)+1):n
      pls.heatmap(exprPathways[rows, ], pathways_weight[rows], plsmodel2_scores_y[col_order, i], 
                  'gene weight', paste0('PLS component-', gsub("Comp ", "", i), ' response score'))
    })
    dev.off()
  } else {
    pdf(paste0("output/heatmap_plsmodel2_", gsub("Comp ", "comp", i), ".pdf"), 13, 11)
    pls.heatmap(exprPathways, pathways_weight, plsmodel2_scores_y[col_order, i], 
                'gene weight', paste0('PLS component-', gsub("Comp ", "", i), ' response score'))
    dev.off()
  }
  
  # Heatmap of top30 pathways
  pdf(paste0("output/heatmap_plsmodel2_", gsub("Comp ", "comp", i), "_top30.pdf"), 12, 5.8)
  hm <- pls.heatmap(exprPathways[1:30, ], pathways_weight[1:30], plsmodel2_scores_y[col_order, i], 
                    'gene weight', paste0('PLS component-', gsub("Comp ", "", i), ' response score'))
  dev.off()
  
})

# Scatter plots to show correlation between component scores of the predictor variables and each response variable
p <- lapply(c(1:2), function(i){
  comp <- paste0('component-', i)
  p <- lapply(colnames(y), function(yi){
    df <- data.frame(name = rownames(y), x = pls_model2$scores[,i], y = y[, yi])
    df$name[which(!((df$x %in% range(df$x) | (df$y %in% range(df$y)))))] <- ""
    r <- round(cor(df$x, df$y), digits = 2)
    l <- format(pls_model2$Yloadings[yi, i], digits = 2, scientific = TRUE)
    p <- ggplot(df, aes(x, y)) +
      geom_point(color = "blue") + geom_smooth(method = "lm") +
      geom_text_repel(aes(label = name), size = 5) +
      labs(x = bquote('PLS'~italic(.(comp))~"of gene expression"),
           y = bquote(atop(italic('t')*'-statistic of relationship between', 'CT and '~.(yi)))) +
      geom_hline(yintercept=0, linetype="dashed", color = "gray") +
      geom_vline(xintercept=0, linetype="dashed", color = "gray") +
      ggtitle(bquote(italic('r')~"="~.(r)~'y-loading ='~.(l))) +
      theme_classic()
    p
  })
  p[["text"]] <- textGrob(bquote('PLS'~italic(.(comp))), gp = gpar(fontsize = 24))
  p <- p[c(10,1:9)]
  p <- grid.arrange(grobs = p, layout_matrix = matrix(c(1:10), ncol = 5, byrow = TRUE))
  p
})
pdf("output/scatterplots_plsmodel2_responses.pdf", 20, 16)
grid.arrange(grobs = p, layout_matrix = matrix(c(1:2), ncol = 1, byrow = TRUE))
grid.lines(x = unit(c(0,1), "npc"), y = unit(c(.5,.5), "npc"))
dev.off()
