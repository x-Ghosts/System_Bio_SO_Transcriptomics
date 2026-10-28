setwd("/home/adam/Desktop/PhD_work_DATA_ANALYSIS_161026/PhD_work/Systems_Biology/SystemsBiology_Brasile/RNA_Seq/step_2_DEG_analysis/blood/")


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
library(apeglm)

######################################################################################################################################################################################
######################################################################################################################################################################################
# Data Processing
######################################################################################################################################################################################
######################################################################################################################################################################################

# ============================================================
# Original gene-count matrix - Annotation set
# ============================================================

gene_count_stress_blood <- read.delim("gene_count_blood/gene_count_blood_ang_nell_stress.txt", header = T, stringsAsFactors = F)
stress_samples_ids<- grepl("^RNA_B_[AN]", colnames(gene_count_stress_blood))
stress_samples_ids_names <- colnames(gene_count_stress_blood[,-c(1,2)])
gene_names <- gene_count_stress_blood[,1:2] 
unique(isUnique(gene_names$gene_id)) # [1] TRUE, else FALSE

matrix_counts <- gene_count_stress_blood[, stress_samples_ids]
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

low_expression_genes <- rownames(matrix_counts)[percentage_low_per_gene > percentage_low_sample*100] #14825
filtered_matrix_counts <- matrix_counts[!rownames(matrix_counts) %in% low_expression_genes,]
nrow(filtered_matrix_counts) #17267 

# Metadata preparation

sample_info <- strcapture("^RNA_B_([AN])([0-9]+)_([^_]+)$", stress_samples_ids_names, proto = list(
  BreedCode = character(),
  AnimalNumber = character(),
  Treatment = character()
)
)

metadata <- data.frame(
  SampleID = stress_samples_ids_names,
  IID = paste0(sample_info$BreedCode, sample_info$AnimalNumber),
  Treatment = sample_info$Treatment,
  Breed = sample_info$BreedCode,
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
metadata$Breed <- factor(metadata$Breed, levels = c("A","N"))


# ============================================================
# Creating DDS object
# ============================================================

# DDS Objects:
stress_blood_dds <- DESeqDataSetFromMatrix(countData = filtered_matrix_counts,
                                           colData = metadata,
                                           design = ~ Breed
)

#stress_blood_dds$Treatment <- relevel(stress_blood_dds$Treatment, ref = "POST") # Right of Volcano's are UP of HS, and left of Volcano's are UP of RECOVERY
stress_blood_dds$Breed <- relevel(stress_blood_dds$Breed, ref = "N") # Right of Volcano's are UP of N, and left of Volcano's are UP of A

######################################################################################################################################################################################
######################################################################################################################################################################################
"DEG analysis"
######################################################################################################################################################################################
######################################################################################################################################################################################
FDR <- 0.05
lfc <- 1

stress_blood_dds_deseq <- DESeq(stress_blood_dds, parallel = TRUE) 
results_stress_blood_dds <- results(stress_blood_dds_deseq)
results_stress_blood_dds_df <- as.data.frame(results_stress_blood_dds)

padjusted_results_stress_blood_dds_df <- results_stress_blood_dds_df[!is.na(results_stress_blood_dds_df$padj),,drop = FALSE] # 17267
names_to_test <- rownames(padjusted_results_stress_blood_dds_df)
test_matrix <- filtered_matrix_counts[rownames(filtered_matrix_counts) %in% names_to_test,]
padjusted_results_stress_blood_dds_df <- padjusted_results_stress_blood_dds_df[padjusted_results_stress_blood_dds_df$padj < FDR & 
                                                                               abs(padjusted_results_stress_blood_dds_df$log2FoldChange) > lfc,] # 1225

######################################################################################################################################################################################
######################################################################################################################################################################################
"Volcanot plotting"
######################################################################################################################################################################################
######################################################################################################################################################################################

fdr_cutoff <- 0.05 #0.1
lfc_cutoff <- 1 #0.58
n_labels <- 35

resultsNames(stress_blood_dds_deseq)
coef_name <- "Breed_A_vs_N" # The statistical comparison of the Design Model ~ Breed

res_stress <- results(stress_blood_dds_deseq, name = coef_name, alpha = fdr_cutoff)
res_stress_shrunk <- lfcShrink(stress_blood_dds_deseq, coef = coef_name, res = res_stress, type = "apeglm")

volcano_df <- as.data.frame(res_stress_shrunk)
volcano_df$gene <- rownames(volcano_df)

volcano_df <- volcano_df[!is.na(volcano_df$padj) &is.finite(volcano_df$log2FoldChange),, drop = FALSE]
volcano_df$neglog10padj <- -log10(pmax(volcano_df$padj, .Machine$double.xmin)) #-log10(FDR)


volcano_df$status <- as.character("Not significant")
sig <- volcano_df$padj < fdr_cutoff

volcano_df$status[sig] <- "Significant (small FC)"
volcano_df$status[sig & volcano_df$log2FoldChange >= lfc_cutoff] <- "UP in EXP"
volcano_df$status[sig & volcano_df$log2FoldChange <= -lfc_cutoff] <- "UP in POST"
volcano_df$status <- factor(volcano_df$status, levels = c("UP in EXP", "UP in POST", "Significant (small FC)", "Not significant"))

volcano_labels <- volcano_df[volcano_df$status %in% c("UP in EXP", "UP in POST"), , drop = FALSE]
volcano_labels <- head( volcano_labels[order(volcano_labels$padj), , drop = FALSE], n_labels)

status_colors <- c(
  "UP in POST" = "#4575B4",
  "Significant (small FC)" = "#E6AB02",
  "Not significant" = "grey75"
)
status_colors["UP in EXP"] <- "#D73027"


volcano_plot_stress <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = neglog10padj)
) +
  geom_point(
    aes(color = status),
    alpha = 0.65,
    size = 1.5
  ) +
  
  scale_color_manual(
    values = status_colors,
    drop = FALSE
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_text_repel(
    data = volcano_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot - stress Blood",
    subtitle = "EXP vs POST | DESeq2",
    x = "log2FoldChange (EXP vs POST)",
    y = "-log10(adjusted p-value)",
    color = "Gene classification"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(volcano_plot_stress)

ggsave("Volcano_stress_Blood.pdf",plot = volcano_plot_stress,width = 10,height = 7)




######################################################################################################################################################################################
######################################################################################################################################################################################
"Heatmap plotting and PCA plotting of the Significant Markers - FC 1 & FDR < 0.05"
######################################################################################################################################################################################
######################################################################################################################################################################################

heatmap_significant_df <- padjusted_results_stress_blood_dds_df
pca_significant_df <- padjusted_results_stress_blood_dds_df

# Heatmap:

gene_names_significant <- rownames(heatmap_significant_df)
normalized_matrix_count <- counts(stress_blood_dds_deseq, normalized = TRUE)
gene_sd <- apply(normalized_matrix_count, 1,sd)
gene_sd_significant <- gene_sd[names(gene_sd) %in% rownames(padjusted_results_stress_blood_dds_df)]
gene_sd_significant <- sort(gene_sd_significant, decreasing = TRUE)

heatmap_matrix_normalized_sig <- log2(
  normalized_matrix_count[names(gene_sd_significant), , drop = FALSE] + 1
)

metadata_heatmap <- metadata[colnames(heatmap_matrix_normalized_sig), , drop = FALSE]
metadata_heatmap$Treatment <- factor(
  metadata_heatmap$Treatment,
  levels = c("EXP", "POST"))
metadata_heatmap$Breed <- factor(
  metadata_heatmap$Breed,
  levels = c("A", "N"))

sample_order <- with(metadata_heatmap, order(Breed, IID))
metadata_heatmap <- metadata_heatmap[sample_order, , drop = FALSE]

heatmap_matrix_normalized_sig <- heatmap_matrix_normalized_sig[, rownames(metadata_heatmap), drop = FALSE]
annotation_col <- metadata_heatmap[, "Breed", drop = FALSE]
gaps <- sum(metadata_heatmap$Breed == "A")


plot_deg_ordered <- pheatmap(
  heatmap_matrix_normalized_sig,
  scale = "row",
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  gaps_col = gaps,
  annotation_col = annotation_col,
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  border_color = NA,
  show_rownames = F,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 90,
  main = "stress Blood - Significant DEGs (EXP vs POST)"
)

plot_deg_clustered <- pheatmap(
  heatmap_matrix_normalized_sig,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  border_color = NA,
  show_rownames = F,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 90,
  main = "stress Blood - Significant DEGs (Clustered)"
)

heatmap_plots <- list("stress_DEG_Ordered.pdf" = plot_deg_ordered, "stress_DEG_Clustered.pdf" = plot_deg_clustered)
for (filename in names(heatmap_plots)) {
  ggsave(filename, plot = heatmap_plots[[filename]]$gtable, width = 11, height = 8)
}



# PCA:

label_size <- 3.5

pca_significant_df <- padjusted_results_stress_blood_dds_df
gene_names_significant <- rownames(heatmap_significant_df)
normalized_matrix_count <- counts(stress_blood_dds_deseq, normalized = TRUE)
pca_matrix_normalized_sig <- log2(normalized_matrix_count[gene_names_significant,] + 1)

t_pca_matrix_deg <- t(pca_matrix_normalized_sig)
pca_stress_deg <- prcomp(t_pca_matrix_deg, center = T, scale. = F)

metadata_pca <- metadata[rownames(t_pca_matrix_deg), , drop = FALSE]
pca_df <- cbind(as.data.frame(pca_stress_deg$x[, 1:4, drop = FALSE]), metadata_pca)
pca_df$Treatment <- factor(pca_df$Treatment,levels = c("EXP", "POST"))
pca_df$IID <- factor(pca_df$IID)
animal_shapes <- setNames(c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6, 9, 10, 11),levels(pca_df$IID))

breed_colors <- c("A" = "#0072B2","N" = "#D55E00")

variance_explained <- 100 * pca_stress_deg$sdev^2 / sum(pca_stress_deg$sdev^2)

variance_df <- data.frame(PC = factor(paste0("PC", seq_along(variance_explained)),levels = paste0("PC", seq_along(variance_explained))),
                          Variance = variance_explained)

plot_pca_deg <- function(pc_x, pc_y) {
  
  ggplot(
    pca_df,
    aes(
      x = .data[[pc_x]],
      y = .data[[pc_y]],
      color = Breed,
      shape = IID
    )
  ) +
    geom_point(size = 2.5, stroke = 1.2) +
    geom_text_repel(
      aes(label = SampleID),
      size = label_size,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = breed_colors) +
    scale_shape_manual(values = animal_shapes) +
    labs(
      title = "PCA - stress Blood (Significant DEGs)",
      subtitle = "FDR < 0.05 | Absolute log2FC >= 1",
      x = paste0(
        pc_x, " (",
        round(variance_explained[
          as.integer(sub("PC", "", pc_x))
        ], 2), "%)"
      ),
      y = paste0(
        pc_y, " (",
        round(variance_explained[
          as.integer(sub("PC", "", pc_y))
        ], 2), "%)"
      ),
      color = "Breed",
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

pca_plot_12 <- plot_pca_deg("PC1", "PC2"); pca_plot_12
pca_plot_34 <- plot_pca_deg("PC3", "PC4"); pca_plot_34

ggsave("PCA_stress_DEGs_PC1_PC2.pdf", pca_plot_12,width = 11, height = 8)
ggsave("PCA_stress_DEGs_PC3_PC4.pdf", pca_plot_34,width = 11, height = 8)


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
    title = "Explained Variance - stress Significant DEGs",
    x = "Principal Components",
    y = "Variance Explained (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(variance_plot)
ggsave("PCA_stress_variance.pdf", variance_plot,width = 11, height = 8)


######################################################################################################################################################################################
######################################################################################################################################################################################
"Heatmap plotting and PCA plotting of the Significant Markers - FC 1 & FDR < 0.05"
######################################################################################################################################################################################
######################################################################################################################################################################################

write.table(padjusted_results_stress_blood_dds_df, "significant_blood_stress_genes.txt", quote = F, sep = "\t", row.names = T)
write.table(results_stress_blood_dds_df, "results_df_blood_stress_genes.txt", quote = F, sep = "\t", row.names = T)
