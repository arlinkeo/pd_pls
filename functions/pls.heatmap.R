# Function to plot heatmap of PLS results
pls.heatmap <- function(expr, row_feature, column_feature, row_name, column_name){
  # column annotation (t-statistics of Delta CT)
  ha_col <- HeatmapAnnotation(a = anno_barplot(column_feature),
                              height = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 8))
  names(ha_col) <- column_name
  # row annotation (gene weights within pathway)
  ha_row <- rowAnnotation(b = anno_boxplot(row_feature), 
                          width = unit(3, "cm"), annotation_name_gp = gpar(fontsize = 8))
  names(ha_row) <- row_name
  pathway_median <- sapply(row_feature, median)
  # Plot heatmap with annotation along the axes
  hm <- Heatmap(expr, name = 'Z-Score expression',
                row_split = ifelse(pathway_median > 0, "Positively\ncorrelated pathways", "Negatively\ncorrelated pathways"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 6),
                row_names_side = c("right"),
                row_title_rot = 90,
                row_title_gp = gpar(fontsize = 8),
                row_title_side = "left",
                column_names_side = c("top"),
                column_names_gp = gpar(fontsize = 6),
                column_names_rot = 45,
                column_title_gp = gpar(fontsize = 8),
                column_title_side = "bottom",
                width = unit(ncol(expr)*.5, "lines"),
                height = unit(nrow(expr)*.5, "lines"),
                top_annotation = ha_col,
                right_annotation = ha_row,
                heatmap_legend_param = list(direction = "horizontal", 
                                            legend_width = unit(3, "cm"), 
                                            title_position = "topcenter",
                                            labels_gp = gpar(fontsize = 8), 
                                            title_gp = gpar(fontsize = 8), fontface = "bold")
  )
  draw(hm, heatmap_legend_side = "top")
  # add vertical line at 0 in row annotation plot
  positive_pathways <- sum(pathway_median>0)
  negative_pathways <- sum(pathway_median<0)
  if (positive_pathways > 0 && negative_pathways > 0){
    decorate_annotation(row_name, {
      grid.lines(c(0, 0), c(-positive_pathways, negative_pathways+.5), default.units = "native",
                 gp = gpar(lty = 1, col = "red"))
    })
  } else if (positive_pathways > 0){
    decorate_annotation(row_name, {
      grid.lines(c(0, 0), c(0, positive_pathways+.5), default.units = "native",
                 gp = gpar(lty = 1, col = "red"))
    })
  } else {
    decorate_annotation(row_name, {
      grid.lines(c(0, 0), c(0, negative_pathways+.5), default.units = "native",
                 gp = gpar(lty = 1, col = "red"))
    })
  }
  
}
