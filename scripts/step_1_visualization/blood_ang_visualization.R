######################################################################################################################################################################################
######################################################################################################################################################################################
"Heatmap Visualization" 
######################################################################################################################################################################################
######################################################################################################################################################################################

# Extract the transformed matrix VST - The method removes the relationship between the mean & variance of count numbers
vst_matrix <- assay(vst_angus_blood_dds)

gene_sd <- apply(vst_matrix, 1, sd)
gene_sd <- sort(gene_sd[is.finite(gene_sd) & gene_sd > 0], decreasing = TRUE)

top_genes <- names(head(gene_sd, 1000))
heatmap_matrix_top <- vst_matrix[top_genes, , drop = FALSE]

all_genes <- names(gene_sd)
heatmap_matrix <- vst_matrix[all_genes, , drop = FALSE]


metadata_heatmap <- metadata[colnames(heatmap_matrix_top), , drop = FALSE]
metadata_heatmap_all <- metadata[colnames(heatmap_matrix), , drop = FALSE]

# meta_heatmap$Breed <- factor(meta_heatmap$Breed,
#                              levels = c("B1", "B2"))
metadata_heatmap$Treatment <- factor(metadata_heatmap$Treatment,
                                 levels = c("EXP", "POST"))
metadata_heatmap_all$Treatment <- factor(metadata_heatmap_all$Treatment,
                                     levels = c("EXP", "POST"))

sample_order <- with(metadata_heatmap, order(Treatment, IID))
metadata_heatmap <- metadata_heatmap[sample_order, , drop = FALSE]

sample_order_all <- order(metadata_heatmap_all$Treatment, metadata_heatmap_all$IID)
metadata_heatmap_all <- metadata_heatmap_all[sample_order_all, , drop = FALSE]
heatmap_matrix <- heatmap_matrix[ , rownames(metadata_heatmap_all), drop = FALSE]

heatmap_matrix_top <- heatmap_matrix_top[, rownames(metadata_heatmap), drop = FALSE]

rownames(metadata_heatmap) == metadata_heatmap$SampleID # T

annotation_col <- metadata_heatmap[, c("Treatment"), drop = FALSE] #annotation_col <- meta_heatmap[, c("Breed", "Treatment"), drop = FALSE]
groups <- paste(metadata_heatmap$Treatment)
gaps <- head(cumsum(rle(groups)$lengths), -1)

annotation_col_all <- metadata_heatmap_all[, "Treatment", drop = FALSE]
groups_all <- paste(metadata_heatmap_all$Treatment)
gaps_all <- sum(metadata_heatmap_all$Treatment == "EXP")

plot_1000_no_cluster <- pheatmap(
  heatmap_matrix_top,
  scale = "row",
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  gaps_col = gaps,
  annotation_col = annotation_col,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 90,
  main = "Top 1000 Variable Genes - VST (Not Clustered - Forcing orders)"
)


plot_1000_cluster <- pheatmap(
  heatmap_matrix_top,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 90,
  main = "Top 1000 Variable Genes - VST (Clustered)"
)

plot_all_no_cluster <- pheatmap(
  heatmap_matrix,
  scale = "row",
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  gaps_col = gaps_all,
  annotation_col = annotation_col_all,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 90,
  main = "All Variable Genes - VST (Not Clustered - Forcing orders)"
)

plot_all_cluster <- pheatmap(
  heatmap_matrix,
  scale = "row",
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  cluster_cols = TRUE,
  gaps_col = gaps_all,
  annotation_col = annotation_col_all,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 90,
  main = "All Variable Genes - VST (Clustered)"
)

heatmap_plots <- list(
  "Heatmap_Top1000_Ordered.pdf" = plot_1000_no_cluster,
  "Heatmap_Top1000_Clustered.pdf" = plot_1000_cluster,
  "Heatmap_All_Ordered.pdf" = plot_all_no_cluster,
  "Heatmap_All_Clustered.pdf" = plot_all_cluster)


for (filename in names(heatmap_plots)) {
  ggsave(filename = filename, plot = heatmap_plots[[filename]]$gtable, width = 11,height = 8,units = "in")
}
