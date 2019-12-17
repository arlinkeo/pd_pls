# Correlation coefficent of clinical scores
library(ggpubr)
library(plyr)
library(abind)
library(reshape2)

# Relation between ct and cs for at most 123 PD patients
coef <- lapply(clinicalScores, function(cs){
  apply(ct_lh[patientIDs, ], 2, function(ct){ # For each region with ct
    res <- lm(ct ~ cs)
    res  <- summary(res)
    res$coefficients[2, c(1,4)]# coefficient + p-value
  })
})
coef <- simplify2array(coef) # 3D-array: measures (estimate, p-value) x regions x clinical features
dimnames(coef)[[2]] <- gsub("lh_", "", dimnames(coef)[[2]])
dimnames(coef)[[3]] <- gsub("_", " ", dimnames(coef)[[3]])
df <- melt(coef["Pr(>|t|)", , ])
df$value <- p.adjust(df$value, method = "BH") # BH-correct for regions and clinical features
bh <- dcast(df, Var1 ~ Var2)[, -1]
coef <- abind(coef, 'BH' = bh, along = 1) # BH added to coef
range(coef["Estimate", , ])
sum(coef["BH", , ] < 0.05)

# Heatmaps for coefficients and BH-corrected P-values
t <- t(coef["Estimate", , ])
hm1 <- Heatmap(t, name = "Beta",
              cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_side = c("left"),
              row_title_rot = 0,
              column_names_side = c("top"),
              column_names_gp = gpar(fontsize = 10),
              width = unit(ncol(t)*.8, "lines"), 
              height = unit(nrow(t)*.8, "lines"),
              heatmap_legend_param = list(title_position = "topcenter"),
              cell_fun = function(i, j, x, y, width, height, fill) {
                if(t(coef["BH", i, j]) < 0.05) {
                  grid.text("*", x = x, y = y)
                } else {
                  grid.text("", x = x, y = y)
                }
              }
)
hm1
pdf("output/coefficients_clinical_scores.pdf", 7.8, 4)
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

# PLS with al features
coef.r <- coef["Estimate", , ]
pls2 <-  plsr(roi_expr ~ coef.r, ncomp = 10, scale = TRUE, validation = "LOO")
explVar <- explvar(pls2)
paste(gsub("Comp ", "PLS", names(explVar)), paste0(round(explVar, digits = 2), "%"), sep = ": ", collapse = "; ")


# Plot
p <- lapply(c(1:3), function(i){
  df <- data.frame(
    roi = parcel_lh$id, 
    name = gsub("lh_", "", parcel_lh$name), 
    p.pls = pls2$scores[,i], # predictor PLS1 scores
    r.pls = pls2$Yscores[,i] # response PLS1 scores
  )
  ggplot(df, aes(p.pls, r.pls)) +
    geom_point(color = "blue") + geom_smooth(method = "lm") +
    geom_text_repel(aes(label = name), size = 3) +
    labs(x = paste0("PLS", i, " X"), y = paste0("PLS", i, " Y")) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", color = "gray") +
    ggtitle(paste("r =", round(cor(df$p.pls, df$r.pls), digits = 2))) +
    theme_classic()
})
pdf("output/pls_clinical_scores.pdf", 18, 5)
ggarrange(plotlist = p, nrow = 1)
dev.off()

# Heatmap top rank genes
hm <- lapply(c(1:3), function(i){
  pls2_coef <- sort(pls2$coefficients[,1,1], decreasing = TRUE) # sort genes by pls1 coefficients
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
})
# Heatmap top genes for PLS1 of X
