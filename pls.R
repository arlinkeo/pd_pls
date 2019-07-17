# Partial Least Squares

library(mixOmics)
library(ggplot2)
library(ggrepel)

# Parcel info
parcel_lh <- parcels[grep("lh-", parcels$Name), ]
parcel_lh$Name <- gsub("ctx-", "", parcel_lh$Name)
parcel_lh$Name <- gsub("-", "_", parcel_lh$Name)

# X: gene expression in cortical regions
roi_expr <- lapply(donorNames, function(d){
  e <- brainExpr[[d]]
  s <- samples_roi[[d]]
  e <- sapply(parcel_lh$ROI, function(r){ # Mean expression within ROI
    cols <- which(s$fs_roi == r)
    if (length(cols) > 1) {
      e <- e[, cols]
      apply(e, 1, mean)
    } else if (length(cols) == 1){
      e[, cols]
    } else {
      rep(NA, nrow(e))
    }
  }) # genes x samples
  colnames(e) <- parcel_lh$ROI
  t(e)
})
roi_expr <- apply(simplify2array(roi_expr), c(1,2), function(x) mean(x, na.rm = TRUE)) # Mean across donors
  
# Y: t-stats of comparing cortical thickness between PD and control,
tstats <- eq2[, c("Group", "t")]
rows <- match(parcel_lh$Name, tstats$Group)
tstats <- tstats[rows, ]
rownames(tstats) <- parcel_lh$ROI
tstats$t <- as.numeric(tstats$t)

# PLS
res <- pls(roi_expr, tstats$t, ncomp = 10, mode = "regression")
# tune <- perf(res, validation = "Mfold", folds = 10)
# plot(tune$Q2.total)
# abline(h = 0.0975)

# Plot
df <- data.frame(
  roi = parcel_lh$ROI, 
  name = gsub("lh_", "", parcel_lh$Name), 
  pls1 = res$variates$X[,1], 
  dct = tstats$t)
p <- ggplot(df, aes(pls1, dct)) +
  geom_point(color = "blue") + geom_smooth(method = "lm") +
  # geom_text(aes(label = name), hjust = -.1, size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  labs(x = "PLS1", y = "t-statistics") +
  ggtitle(paste("r =", round(cor(df$pls1, df$dct), digits = 2))) +
  theme_classic()
pdf("output/PLS1_tstats.pdf", 8, 6)
p
dev.off()

# Top ranked genes per pls component
loadings <- res$loadings$X[,1]
loadings <- sort(loadings)
topgenes <- names(loadings)[1:100]
write.table(topgenes, "output/topgenes.txt", quote = FALSE, row.names = FALSE)
