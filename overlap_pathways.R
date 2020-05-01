# Overlapping pathways

pl <- list(
  p2.1 = gsea2$`Comp 1`@result$Description,
  p2.2 = gsea2$`Comp 2`@result$Description,
  p1 = gsea1@result$Description
)
lengths(pl)
names(pl) <-  c("model-2\ncomponent-1", "model-2\ncomponent-2", "model-1\ncomponent-1")
venn(pl)

# Table
all_pathways <- unique(unlist(pl))
table <- t(sapply(all_pathways, function(p){
  sapply(pl, function(l){
    (as.integer(p %in% l))
  })
}))
pdf("output/overlapping_pathways.pdf", 6, 8)
Heatmap(table, 
        name = "enriched pathway",
        col = c("lightgray", "darkgreen"),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 2),
        column_names_rot = 45,
        width = unit(ncol(table)*3, "lines"),
        height = unit(nrow(table)*.1, "lines")
        )
dev.off()

