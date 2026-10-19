setwd("/home/adam/Desktop/PhD_work_DATA_ANALYSIS_161026/PhD_work/Systems_Biology/SystemsBiology_Brasile/RNA_Seq/step_1_data_analysis/blood/")

# Libraries :
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))



library(dplyr)
library(tidyverse)
library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify) #install.packages("ggfortify")
library(stats)
library(sva)
library(pheatmap)

######################################################################################################################################################################################
######################################################################################################################################################################################
# Data Processing
######################################################################################################################################################################################
######################################################################################################################################################################################

# ============================================================
# Original gene-count matrix - Annotation set
# ============================================================

gene_count_angus_blood <- read.table("gene_count_blood/gene_count_blood_angus.txt", header = T, sep = "\t", stringsAsFactors = F)
angus_samples_ids<- grepl("^RNA_B_A", colnames(gene_count_angus_blood))
angus_samples_ids_names <- colnames(gene_count_angus_blood[,-c(1,2)])
gene_names <- gene_count_angus_blood[,1:2] 
unique(isUnique(gene_names$gene_id)) # [1] TRUE, else FALSE

matrix_counts <- gene_count_angus_blood[, angus_samples_ids]
matrix_counts <- as.matrix(matrix_counts, drop=FALSE)
rownames(matrix_counts) <- gene_names$gene_id

# Checking for . or , patterns:

num_of_dots <- grep(".", matrix_counts, fixed = TRUE, value = TRUE)
num_of_commas <- grep(",", matrix_counts, value = TRUE)
matrix_counts <- as.matrix(matrix_counts, drop=FALSE)

matrix_counts[] <- gsub(",",".", matrix_counts, fixed = TRUE) #First do: matrix_counts[] <- gsub(".","", matrix_counts, fixed = TRUE) then apply this line
storage.mode(matrix_counts) <- "numeric"

matrix_counts[] <- round(matrix_counts, digits = 0)
dim(matrix_counts)
str(matrix_counts)

### QC checks
typeof(matrix_counts) # double
sum(is.na(matrix_counts)) # checking na's
all(matrix_counts >= 0) # detect negative counts
range(matrix_counts)
quantile(matrix_counts, probs = seq(0,1,0.1)) # % of 0's in the data "mean(matrix_counts == 0) * 100"

# Zero Inflated Values (The commented was a filter on ZERO INFLATED VALUES)
# num_0_per_gene <- rowSums(matrix_counts == 0); View(gene_rows_0)
# percentage_0_gene <- rowMeans(matrix_counts == 0) * 100; View(percentage_0_gene)
# summary(percentage_0_gene)
# 
# 
# empty_genes <- rownames(matrix_counts[num_0_per_gene == ncol(matrix_counts),])
# genes_0_25pct <- rownames(matrix_counts[percentage_0_gene[] >= 75,])
# inflated_zeros_gene <- unique(c(empty_genes, genes_0_25pct))
# 
# filtered_matrix_counts <- matrix_counts[!rownames(matrix_counts) %in% inflated_zeros_gene,] ## GENE CLEANUP
# 
# num_0_per_sample <- colSums(filtered_matrix_counts == 0); View(num_0_per_sample)
# percentage_0_sample <- colMeans(filtered_matrix_counts == 0) * 100; View(percentage_0_sample)
# summary(percentage_0_sample)

# Low Count Removal:

low_count_threshold <- 3
percentage_low_sample <- 0.8 # It means I chose the threshold of removing a gene were the observation of 3 is at least observed in 80% of my samples

num_low_per_gene <- rowSums(matrix_counts <= low_count_threshold) # Using head(), a table is stdout to explain you each variable with its count of samples considered below the threshold.
percentage_low_per_gene <- rowMeans(matrix_counts <= low_count_threshold) * 100

low_expression_genes <- rownames(matrix_counts)[percentage_low_per_gene > percentage_low_sample*100] #14569
filtered_matrix_counts <- matrix_counts[!rownames(matrix_counts) %in% low_expression_genes,]
nrow(filtered_matrix_counts) #17523 

# Metadata preparation

sample_info <- strcapture("^RNA_B_([AN])([0-9]+)_([^_]+)$", angus_samples_ids_names, proto = list(
  BreedCode = character(),
  AnimalNumber = character(),
  Treatment = character()
)
)

metadata <- data.frame(
  SampleID = angus_samples_ids_names,
  IID = paste0(sample_info$BreedCode, sample_info$AnimalNumber),
  Treatment = sample_info$Treatment,
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$SampleID


# ============================================================
# All variables datasets
# ============================================================
filtered_matrix_counts
gene_names
metadata


######################################################################################################################################################################################
######################################################################################################################################################################################
"Data Processing - set Factors cols in metadata and DDS object creation"
######################################################################################################################################################################################
######################################################################################################################################################################################

# ============================================================
# Factor levels
# ============================================================

metadata$Treatment <- gsub("D4", "EXP", metadata$Treatment)
metadata$Treatment <- factor(metadata$Treatment, levels = c("EXP", "POST"))
metadata$IID <- factor(metadata$IID)

# ============================================================
# Creating DDS object
# ============================================================

# DDS Objects:
angus_blood_dds <- DESeqDataSetFromMatrix(countData = filtered_matrix_counts,
                                          colData = metadata,
                                          design = ~ IID + Treatment
)

angus_blood_dds$Treatment <- relevel(angus_blood_dds$Treatment, ref = "POST") # Right of Volcano's are UP of HS, and left of Volcano's are UP of RECOVERY

######################################################################################################################################################################################
######################################################################################################################################################################################
"DEG analysis"
######################################################################################################################################################################################
######################################################################################################################################################################################

angus_blood_dds_deseq <- DESeq(angus_blood_dds, parallel = TRUE) 
results_angus_blood_dds <- results(angus_blood_dds_deseq)
results_angus_blood_dds_df <- as.data.frame(results_angus_blood_dds)


######################################################################################################################################################################################
######################################################################################################################################################################################
"PCA analysis - PCA df" 
######################################################################################################################################################################################
######################################################################################################################################################################################

shape_column <- "Treatment"
label <- TRUE
label_size <- 3.5 # How big is the labels

# Variance Stabilizing Transformation
vst_angus_blood_dds <- vst(angus_blood_dds, blind = TRUE)
t_filtered_matrix_counts <- t(assay(vst_angus_blood_dds))
metadata_pca <- metadata[rownames(t_filtered_matrix_counts), ,drop = FALSE]

# Computing Principal Components
pca_angus_blood_dds <- prcomp(t_filtered_matrix_counts, scale. = F, center = T)

# PLotting Principal Components

pca_plot <- autoplot(
  pca_angus_blood_dds,
  data = metadata_pca,
  shape = shape_column,
  label = FALSE
) +
  geom_text_repel(
    aes(label = SampleID),
    size = label_size
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "PCA - Rounded Counts (VST)",
    shape = "Treatment"
  )
ggsave(filename = "PCA_angus_blood.pdf",plot = pca_plot,width = 8,height = 6)

######################################################################################################################################################################################
######################################################################################################################################################################################
"PCA analysis - PCA df - test 2" 
######################################################################################################################################################################################
######################################################################################################################################################################################

shape_column <- "Treatment"
label_size <- 3.5

# Variance Stabilizing Transformation
vst_angus_blood_dds <- vst(angus_blood_dds, blind = TRUE)

t_filtered_matrix_counts <- t(assay(vst_angus_blood_dds))

metadata_pca <- metadata[
  rownames(t_filtered_matrix_counts), , drop = FALSE
]

# Validate sample alignment
stopifnot(
  identical(rownames(metadata_pca),
            rownames(t_filtered_matrix_counts)),
  !anyNA(metadata_pca$SampleID),
  !anyNA(metadata_pca$IID),
  !anyNA(metadata_pca$Treatment)
)


# Computing Principal Components
pca_angus_blood_dds <- prcomp(t_filtered_matrix_counts, scale. = FALSE, center = TRUE)

# PCA scores + metadata
pca_df <- cbind(as.data.frame(pca_angus_blood_dds$x[, 1:nrow(pca_angus_blood_dds$x), drop = FALSE]),metadata_pca)
pca_df$Treatment <- factor(pca_df$Treatment)
pca_df$IID <- factor(pca_df$IID)

animal_shapes <- setNames(c(16, 17, 15, 18, 3, 4, 7, 8),levels(pca_df$IID))
treatment_colors <- c("EXP" = "#0072B2","POST" = "#D55E00")


# Percentage of variance explained
variance_explained <- 100 *
  pca_angus_blood_dds$sdev^2 /
  sum(pca_angus_blood_dds$sdev^2)

# Create variance-explained dataframe
variance_df <- data.frame(
  PC = factor(
    paste0("PC", seq_along(variance_explained)),
    levels = paste0("PC", seq_along(variance_explained))
  ),
  Variance = variance_explained
)

# PCA plot
plot_pca <- function(pc_x, pc_y) {
  
  ggplot(
    pca_df,
    aes(
      x = .data[[pc_x]],
      y = .data[[pc_y]],
      color = Treatment,
      shape = IID
    )
  ) +
    geom_point(size = 2.5, stroke = 1.2) +
    geom_text_repel(
      aes(label = SampleID),
      size = 3.5,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = treatment_colors) +
    scale_shape_manual(values = animal_shapes) +
    labs(
      title = "PCA - Angus Blood (VST)",
      x = paste0(
        pc_x, " (",
        round(variance_explained[as.integer(sub("PC", "", pc_x))], 2),
        "%)"
      ),
      y = paste0(
        pc_y, " (",
        round(variance_explained[as.integer(sub("PC", "", pc_y))], 2),
        "%)"
      ),
      color = "Treatment",
      shape = "Animal ID"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 17, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11)
    )
}

pca_plot_12 <- plot_pca("PC1", "PC2")
pca_plot_34 <- plot_pca("PC3", "PC4")

print(pca_plot_12)
print(pca_plot_34)

ggsave(filename = "PCA_angus_blood_PC1_PC2.pdf",plot = pca_plot_12,width = 11,height = 8,units = "in")
ggsave(filename = "PCA_angus_blood_PC3_PC4.pdf",plot = pca_plot_34,width = 11,height = 8,units = "in")



# Plot variance explained
variance_plot <- ggplot(
  variance_df,
  aes(x = PC, y = Variance)
) +
  geom_col(fill = "orange", width = 0.7) +
  geom_text(
    aes(label = paste0(round(Variance, 1), "%")),
    vjust = -0.4,
    size = 3.2
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Explained Variance - Angus Blood",
    x = "Principal Components",
    y = "Variance Explained (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(variance_plot)

ggsave(
  filename = "PCA_variance_angus_blood.pdf",
  plot = variance_plot,
  width = 10,
  height = 6,
  units = "in"
)


######################################################################################################################################################################################
######################################################################################################################################################################################
# Investiguating Why weird PCA's
######################################################################################################################################################################################
######################################################################################################################################################################################

## Create QC dataframe, and we check for the CORR(PC1/2 vs Library SIZE)
# if (is.null(sizeFactors(angus_blood_dds))) {
#   angus_blood_dds <- estimateSizeFactors(angus_blood_dds)
# }
# 
# qc_samples <- data.frame(
#          SampleID = colnames(angus_blood_dds),
#          LibrarySize = colSums(counts(angus_blood_dds)),
#          SizeFactor = sizeFactors(angus_blood_dds))
# print(qc_samples)
# 
# pca_scores <- pca_angus_blood_dds$x[qc_samples$SampleID, c("PC1","PC2")]
# correlation_libsize_PC1 <- cor(pca_scores[,"PC1"], log10(qc_samples$LibrarySize)) # -0.18
# correlation_libsize_PC2 <- cor(pca_scores[,"PC2"], log10(qc_samples$LibrarySize)) # -0.48
# 
# # SECOND CHECK for the A19 sample:
# detected_genes <- colSums(counts(angus_blood_dds) > 0)
# qc_detected <- data.frame(
#   SampleID = qc_samples$SampleID,
#   Detected_Genes = detected_genes[qc_samples$SampleID],
#   PC2 = pca_angus_blood_dds$x[qc_samples$SampleID, "PC2"]
# )
# 
# qc_detected$A19 <- grepl("_A19_", qc_detected$SampleID)
# qc_detected[qc_detected$A19, ]
# 
# correlation_detected_PC2 <- cor(
#   qc_detected$Detected_Genes,
#   qc_detected$PC2
# )
# 
# correlation_without_A19 <- cor(
#   qc_detected$Detected_Genes[!qc_detected$A19],
#   qc_detected$PC2[!qc_detected$A19]
# )
# 
# correlation_detected_PC2 # -0.64
# correlation_without_A19 # -0.48

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

sample_order_all <- with(metadata_heatmap_all, order(Treatment, IID))
metadata_heatmap_all <- metadata_heatmap_all[sample_order_all, , drop = FALSE]

rownames(metadata_heatmap) == metadata_heatmap$SampleID # T

annotation_col <- metadata_heatmap[, c("Treatment"), drop = FALSE] #annotation_col <- meta_heatmap[, c("Breed", "Treatment"), drop = FALSE]
groups <- paste(metadata_heatmap$Treatment)
gaps <- head(cumsum(rle(groups)$lengths), -1)

annotation_col_all <- metadata_heatmap_all[, c("Treatment"), drop = FALSE] #annotation_col <- meta_heatmap[, c("Breed", "Treatment"), drop = FALSE]
groups_all <- paste(metadata_heatmap_all$Treatment)
gaps_all <- head(cumsum(rle(groups_all)$lengths), -1)

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
