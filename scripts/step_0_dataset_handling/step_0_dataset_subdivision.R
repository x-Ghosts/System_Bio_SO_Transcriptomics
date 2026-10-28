setwd("/home/adam/Desktop/PhD_work_DATA_ANALYSIS_161026/PhD_work/Systems_Biology/SystemsBiology_Brasile/RNA_Seq/step_0_dataset_subdivision/")

#X_data <- read.delim("transcript_count_blood/transcripts_liver_angus.txt")

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


# Loading datasets:

transcripts <- read.delim("transcript_counts_S.tsv", header = T, stringsAsFactors = F)
gene <- read.delim("gene_counts_S.tsv", header = T, stringsAsFactors = F)

#########################################################################################################################################################################
#########################################################################################################################################################################
# GENE COUNT
##########################################################################################################################################################################
#########################################################################################################################################################################

# Blood only:
blood_cols <- grepl("^RNA_B_", colnames(gene))
blood_cols_gene <- grepl("^RNA_B_", colnames(gene))
transcripts_blood <- transcripts[,c("tx", "gene_id", colnames(transcripts)[blood_cols])]
gene_count_blood <- gene[,c("gene_id", "gene_name", colnames(gene)[blood_cols])]

# Liver only:
liver_cols <- grepl("^RNA_L_", colnames(transcripts))
liver_cols_gene <- grepl("^RNA_L_", colnames(gene))
transcripts_liver <- transcripts[,c("tx", "gene_id", colnames(transcripts)[blood_cols])]
gene_count_liver <- gene[,c("gene_id", "gene_name", colnames(gene)[blood_cols])]



# Blood Angus:
blood_ang_cols <- grepl("^RNA_B_A", colnames(gene))
gene_count_blood_angus <- gene[,c("gene_id", "gene_name", colnames(gene)[blood_ang_cols])]

# Blood Nellore:
blood_nel_cols <- grepl("^RNA_B_N", colnames(gene))
gene_count_blood_nellore <- gene[,c("gene_id", "gene_name", colnames(gene)[blood_nel_cols])]

# Liver Angus:
liver_ang_cols <- grepl("^RNA_L_A", colnames(gene))
gene_count_liver_angus <- gene[,c("gene_id", "gene_name", colnames(gene)[liver_ang_cols])]

# Liver Nellore:
liver_nel_cols <- grepl("^RNA_L_N", colnames(gene))
gene_count_liver_nellore <- gene[,c("gene_id", "gene_name", colnames(gene)[liver_nel_cols])]

# Stressed Nellore & Stressed Angus together "BLOOD"
stress_blood_cols <- grepl("^RNA_B_.*_D4$", colnames(gene))
gene_count_stress_only_blood <- gene[,c("gene_id", "gene_name", colnames(gene)[stress_blood_cols])]

# Non Stressed Nellore & Stressed Angus together "BLOOD"
no_stress_blood_cols <- grepl("^RNA_B_.*_POST$", colnames(gene))
gene_count_no_stress_only_blood <- gene[,c("gene_id", "gene_name", colnames(gene)[no_stress_blood_cols])]


# Stressed Nellore & Stressed Angus together "LIVER"
stress_liver_cols <- grepl("^RNA_L_.*_D4$", colnames(gene))
gene_count_stress_only_liver <- gene[,c("gene_id", "gene_name", colnames(gene)[stress_liver_cols])]


# Non Stressed Nellore & Stressed Angus together "LIVER"
no_stress_liver_cols <- grepl("^RNA_L_.*_POST$", colnames(gene))
gene_count_no_stress_only_liver <- gene[,c("gene_id", "gene_name", colnames(gene)[no_stress_liver_cols])]


write.table(gene_count_blood, "gene_count_blood_all_animals.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_liver, "gene_count_liver_all_animals.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_blood_angus, "gene_count_blood/gene_count_blood_angus.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_blood_nellore, "gene_count_blood/gene_count_blood_nellore.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_liver_angus, "gene_count_liver/gene_count_liver_angus.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_liver_nellore, "gene_count_liver/gene_count_liver_nellore.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_stress_only_blood, "gene_count_blood/gene_count_blood_ang_nell_stress.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_stress_only_liver, "gene_count_liver/gene_count_liver_ang_nell_stress.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_no_stress_only_blood, "gene_count_blood/gene_count_blood_ang_nell_NOstress.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene_count_no_stress_only_liver, "gene_count_liver/gene_count_liver_ang_nell_NOstress.txt", quote = F, sep = "\t", row.names = F, col.names = T)


#########################################################################################################################################################################
#########################################################################################################################################################################
# TRANSCRIPT COUNT
##########################################################################################################################################################################
#########################################################################################################################################################################

# Blood only:
blood_cols_transcripts <- grepl("^RNA_B_", colnames(transcripts))
transcripts_blood <- transcripts[,c("tx", "gene_id", colnames(transcripts)[blood_cols_transcripts])]

# Liver only:
liver_cols_transcripts <- grepl("^RNA_L_", colnames(transcripts))
transcripts_liver <- transcripts[,c("tx", "gene_id", colnames(transcripts)[liver_cols_transcripts])]

# Blood Angus Transcripts:
blood_ang_transcripts_cols <- grepl("^RNA_B_A", colnames(transcripts))
transcripts_blood_angus <- transcripts[,c("tx", "gene_id", colnames(transcripts)[blood_ang_transcripts_cols])]

# Blood Nellore Transcripts:
blood_nel_transcripts_cols <- grepl("^RNA_B_N", colnames(transcripts))
transcripts_blood_nellore <- transcripts[,c("tx", "gene_id", colnames(transcripts)[blood_nel_transcripts_cols])]

# Liver Angus Transcripts:
liver_ang_transcripts_cols <- grepl("^RNA_L_A", colnames(transcripts))
transcripts_liver_angus <- transcripts[,c("tx", "gene_id", colnames(transcripts)[liver_ang_transcripts_cols])]

# Liver Nellore Transcripts:
liver_nel_transcripts_cols <- grepl("^RNA_L_N", colnames(transcripts))
transcripts_liver_nellore <- transcripts[,c("tx", "gene_id", colnames(transcripts)[liver_nel_transcripts_cols])]

# Stressed Nellore & Stressed Angus together "BLOOD"
stress_blood_cols <- grepl("^RNA_B_.*_D4$", colnames(transcripts))
transcripts_stress_only_blood <- transcripts[,c("tx", "gene_id", colnames(transcripts)[stress_blood_cols])]

# Non Stressed Nellore & Stressed Angus together "BLOOD"
no_stress_blood_cols <- grepl("^RNA_B_.*_POST$", colnames(transcripts))
transcripts_no_stress_only_blood <- transcripts[,c("tx", "gene_id", colnames(transcripts)[no_stress_blood_cols])]


# Stressed Nellore & Stressed Angus together "LIVER"
stress_liver_cols <- grepl("^RNA_L_.*_D4$", colnames(transcripts))
transcripts_stress_only_liver <- transcripts[,c("tx", "gene_id", colnames(transcripts)[stress_liver_cols])]

# Non Stressed Nellore & Stressed Angus together "LIVER"
no_stress_liver_cols <- grepl("^RNA_L_.*_POST$", colnames(transcripts))
transcripts_no_stress_only_liver <- transcripts[,c("tx", "gene_id", colnames(transcripts)[no_stress_liver_cols])]



write.table(transcripts_blood, "transcripts_blood_all_animals.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_liver, "transcripts_liver_all_animals.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_blood_angus, "transcript_count_blood/transcripts_liver_angus.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_blood_nellore, "transcript_count_blood/transcripts_blood_nellore.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_liver_angus, "transcript_count_liver/transcripts_liver_angus.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_liver_nellore, "transcript_count_liver/transcripts_liver_nellore.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_stress_only_blood, "transcript_count_blood/transcripts_blood_ang_nell_stress.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_stress_only_liver, "transcript_count_liver/transcripts_liver_ang_nell_stress.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_no_stress_only_blood, "transcript_count_blood/transcripts_blood_ang_nell_NO_stress.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(transcripts_no_stress_only_liver, "transcript_count_liver/transcripts_liver_ang_nell_NO_stress.txt", quote = F, sep = "\t", row.names = F, col.names = T)

