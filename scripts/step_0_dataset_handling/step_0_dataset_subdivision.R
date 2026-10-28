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
library(tidyr)
library(stringr)
library(patchwork)


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

#########################################################################################################################################################################
#########################################################################################################################################################################
# RNA material quality
##########################################################################################################################################################################
#########################################################################################################################################################################

blood <- read.delim("RNA_blood.csv", header = T, stringsAsFactors = F)
liver <- read.delim("RNA_liver.csv", header = T, stringsAsFactors = F)
muscle <- read.delim("RNA_Muscle.csv", header = T, stringsAsFactors = F)
skin <- read.delim("RNA_Skin.csv", header = T, stringsAsFactors = F)


all_animals <- blood %>%
  mutate(
    Sample = str_squish(as.character(Sample)),
    Animal = sub("\\s+(D4|POST)$", "", Sample)
  ) %>%
  pull(Animal) %>%
  unique()

plot_RIN <- function(df, tissue_name) {
  
  # Extract animal and treatment
  rin_df <- df %>%
    mutate(
      Sample = str_squish(as.character(Sample)),
      Animal = sub("\\s+(D4|POST)$", "", Sample),
      Period = sub("^.*\\s+(D4|POST)$", "\\1", Sample),
      RIN = suppressWarnings(as.numeric(as.character(RIN)))
    ) %>%
    select(Animal, Period, RIN)
  
  # Include ALL animals and BOTH periods
  rin_complete <- expand_grid(
    Animal = all_animals,
    Period = c("D4", "POST")
  ) %>%
    left_join(rin_df, by = c("Animal", "Period")) %>%
    mutate(
      Animal = factor(Animal, levels = all_animals),
      Period = factor(Period, levels = c("D4", "POST")),
      Missing = is.na(RIN),
      Bar_RIN = coalesce(RIN, 0),
      Label = if_else(Missing, "N/A", sprintf("%.1f", RIN)),
      Label_Y = if_else(Missing, 0.4, Bar_RIN + 0.25)
    )
  
  ggplot(rin_complete, aes(x = Animal, y = Bar_RIN, fill = Period)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black"
    ) +
    geom_text(
      aes(y = Label_Y, label = Label),
      position = position_dodge(width = 0.8),
      size = 3
    ) +
    scale_fill_manual(
      values = c("D4" = "steelblue", "POST" = "tomato")
    ) +
    scale_y_continuous(
      limits = c(0, 11),
      breaks = seq(0, 10, 2)
    ) +
    labs(
      title = paste("RIN Distribution -", tissue_name),
      x = "Animal",
      y = "RNA Integrity Number (RIN)",
      fill = "Period"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

p1 <- plot_RIN(blood, "Blood")
p2 <- plot_RIN(liver, "Liver")
p3 <- plot_RIN(skin, "Skin")
p4 <- plot_RIN(muscle, "Muscle")
combined_RIN <- (p1 | p2) / (p3 | p4)

# Display
combined_RIN

ggsave(filename = "RIN_All_Tissues.pdf",
  plot = combined_RIN,
  width = 20,
  height = 12,
  bg = "white"
)
