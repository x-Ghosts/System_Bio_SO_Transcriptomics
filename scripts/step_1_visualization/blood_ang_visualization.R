setwd("/home/adam/Desktop/PhD_work/Systems_Biology/SystemsBiology_Brasile/RNA_Seq/step_1_data_analysis/blood/angus/")

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

# Zero Inflated Values
num_0_per_gene <- rowSums(matrix_counts == 0); View(gene_rows_0)
percentage_0_gene <- rowMeans(matrix_counts == 0) * 100; View(percentage_0_gene)
summary(percentage_0_gene)


empty_genes <- rownames(matrix_counts[num_0_per_gene == ncol(matrix_counts),])
genes_0_25pct <- rownames(matrix_counts[percentage_0_gene[] >= 75,])
inflated_zeros_gene <- unique(c(empty_genes, genes_0_25pct))

filtered_matrix_counts <- matrix_counts[!rownames(matrix_counts) %in% inflated_zeros_gene,] ## GENE CLEANUP

num_0_per_sample <- colSums(filtered_matrix_counts == 0); View(num_0_per_sample)
percentage_0_sample <- colMeans(filtered_matrix_counts == 0) * 100; View(percentage_0_sample)
summary(percentage_0_sample)


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
                                          colData = meta,
                                          design = ~ IID + Treatment
)

dds_rounded$Treatment <- relevel(dds_rounded$Treatment, ref = "POST")
