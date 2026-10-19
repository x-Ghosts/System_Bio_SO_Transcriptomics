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

angus_blood_dds_deseq <- DESeq(angus_blood_dds, parallel = TRUE) # 17523
results_angus_blood_dds <- results(angus_blood_dds_deseq)
results_angus_blood_dds_df <- as.data.frame(results_angus_blood_dds)

padjusted_results_angus_blood_dds_df <- results_angus_blood_dds_df[!is.na(results_angus_blood_dds_df$padj),,drop = FALSE] # 11068
padjusted_results_angus_blood_dds_df <- padjusted_results_angus_blood_dds_df[padjusted_results_angus_blood_dds_df$padj < 0.05,] # 43

######################################################################################################################################################################################
######################################################################################################################################################################################
"Volcanot plotting"
######################################################################################################################################################################################
######################################################################################################################################################################################

fdr_cutoff <- 0.05
lfc_cutoff <- 1
n_labels <- 25

resultsNames(angus_blood_dds_deseq)
coef_name <- "Treatment_EXP_vs_POST" # The statistical comparison of the Design Model ~ IID + Treatment

res_angus <- results(angus_blood_dds_deseq, name = coef_name, alpha = fdr_cutoff)
res_angus_shrunk <- lfcShrink(angus_blood_dds_deseq, coef = coef_name, res = res_angus, type = "apeglm")

volcano_df <- as.data.frame(res_angus_shrunk)
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

volcano_plot_angus <- ggplot(
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
    title = "Volcano Plot - Angus Blood",
    subtitle = paste("EXP vs POST | DESeq2"),
    x = paste0("Shrunken log2FoldChange ( EXP vs POST)"),
    y = "-log10(adjusted p-value)",
    color = "Gene classification"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(volcano_plot_angus)

ggsave("Volcano_Angus_Blood.pdf",plot = volcano_plot_angus,width = 10,height = 7)


# Check the normalized expression
ido1_counts <- plotCounts(
  angus_blood_dds_deseq,
  gene = "IDO1",
  intgroup = c("IID", "Treatment"),
  returnData = TRUE
)

ggplot(
  ido1_counts,
  aes(x = Treatment, y = count, group = IID)
) +
  geom_line(color = "white") +
  geom_point(aes(color = IID), size = 3) +
  scale_y_log10() +
  theme_classic(base_size = 13) +
  labs(
    title = "IDO1 - Normalized Expression",
    x = "Treatment",
    y = "Normalized Counts (log10)"
  )


