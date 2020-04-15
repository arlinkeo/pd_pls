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
hm1 <- Heatmap(t, name = "t-statistic\nof beta",
               col = col_fun,
               cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10),
               row_names_side = c("left"),
               row_title_rot = 0,
               column_names_side = c("top"),
               column_names_gp = gpar(fontsize = 10),
               width = unit(ncol(t)*.8, "lines"), 
               height = unit(nrow(t)*.8, "lines"),
               cell_fun = function(i, j, x, y, width, height, fill) {
                 if(t(coef_eq1["BH", i, j]) < 0.05) {
                   grid.text("*", x = x, y = y)
                 } else {
                   grid.text("", x = x, y = y)
                 }
               }
)
t <- t(coef_eq1["Estimate", , ])
q1 <- max(quantile(abs(t), 0.9))
col_fun <- colorRamp2(c(-q1, 0, q1), c("blue", "#EEEEEE", "red"))# Heat colors centered around 0
hm2 <- Heatmap(t, name = "beta",
               col = col_fun,
               cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10),
               row_names_side = c("left"),
               row_title_rot = 0,
               column_names_side = c("top"),
               column_names_gp = gpar(fontsize = 10),
               width = unit(ncol(t)*.8, "lines"),
               height = unit(nrow(t)*.8, "lines"),
               cell_fun = function(i, j, x, y, width, height, fill) {
                 if(t(coef_eq1["BH", i, j]) < 0.05) {
                   grid.text("*", x = x, y = y)
                 } else {
                   grid.text("", x = x, y = y)
                 }
               }
)
t <- scale(t)
q1 <- max(quantile(abs(t), 0.9))
col_fun <- colorRamp2(c(-q1, 0, q1), c("blue", "#EEEEEE", "red"))# Heat colors centered around 0
hm3 <- Heatmap(t, name = "scaled beta\n(for each column)",
               col = col_fun,
               cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10),
               row_names_side = c("left"),
               row_title_rot = 0,
               column_names_side = c("top"),
               column_names_gp = gpar(fontsize = 10),
               width = unit(ncol(t)*.8, "lines"),
               height = unit(nrow(t)*.8, "lines"),
               cell_fun = function(i, j, x, y, width, height, fill) {
                 if(t(coef_eq1["BH", i, j]) < 0.05) {
                   grid.text("*", x = x, y = y)
                 } else {
                   grid.text("", x = x, y = y)
                 }
               }
)

t <- t(coef_eq1["Estimate", , ])
t <- t(scale(t(t)))
q1 <- max(quantile(abs(t), 0.9))
col_fun <- colorRamp2(c(-q1, 0, q1), c("blue", "#EEEEEE", "red"))# Heat colors centered around 0
hm4 <- Heatmap(t, name = "scaled beta\n(for each row)",
               col = col_fun,
               cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10),
               row_names_side = c("left"),
               row_title_rot = 0,
               column_names_side = c("top"),
               column_names_gp = gpar(fontsize = 10),
               width = unit(ncol(t)*.8, "lines"),
               height = unit(nrow(t)*.8, "lines"),
               cell_fun = function(i, j, x, y, width, height, fill) {
                 if(t(coef_eq1["BH", i, j]) < 0.05) {
                   grid.text("*", x = x, y = y)
                 } else {
                   grid.text("", x = x, y = y)
                 }
               }
)
t <- scale(coef_eq1["Estimate", , ])
t <- scale(t(t))
q1 <- max(quantile(abs(t), 0.9))
col_fun <- colorRamp2(c(-q1, 0, q1), c("blue", "#EEEEEE", "red"))# Heat colors centered around 0
hm5 <- Heatmap(t, name = "scaled beta\n(for each row\nand column)",
               col = col_fun,
               cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10),
               row_names_side = c("left"),
               row_title_rot = 0,
               column_names_side = c("top"),
               column_names_gp = gpar(fontsize = 10),
               width = unit(ncol(t)*.8, "lines"), 
               height = unit(nrow(t)*.8, "lines"),
               cell_fun = function(i, j, x, y, width, height, fill) {
                 if(t(coef_eq1["BH", i, j]) < 0.05) {
                   grid.text("*", x = x, y = y)
                 } else {
                   grid.text("", x = x, y = y)
                 }
               }
)
pdf("output/heatmap_clinical_scores.pdf", 8.5, 3.8)#7.8, 3.8)
hm1
hm2
hm3
hm4
hm5
dev.off()