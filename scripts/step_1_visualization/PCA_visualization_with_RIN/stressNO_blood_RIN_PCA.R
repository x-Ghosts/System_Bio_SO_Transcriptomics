library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggnewscale)


treatment_colors <- c(
  "POST" = "#0072B2",  
  "EXP"   = "#D55E00"   
)

pca_results <- readRDS(
  "PCA_Nostress_Blood_Results.rds"
)

blood <- read.delim(
  "../../step_0_dataset_subdivision/Wetlab_RNA_quality$/RNA_blood.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "N/A", "")
)
pca_df <- pca_results$pca_df
variance_explained <- pca_results$variance_explained


# Using your existing blood dataframe
rna_quality_blood <- blood %>%
  mutate(
    Sample = str_squish(as.character(Sample)),
    SampleID = paste0(
      "RNA_B_",
      gsub("\\s+", "_", Sample)
    ),
    RIN = as.numeric(as.character(RIN)),
    ng.ul = as.numeric(as.character(ng.ul))
  ) %>%
  select(SampleID, RIN, ng.ul)

# Join RNA quality with PCA
pca_quality_df <- pca_df %>%
  left_join(
    rna_quality_blood,
    by = "SampleID"
  )

# Verify matching
pca_quality_df %>%
  select(SampleID, IID, Treatment, RIN, ng.ul)


animal_shapes <- setNames(
  c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6, 9, 10, 11),
  levels(factor(pca_quality_df$IID))
)

plot_pca_quality <- function(variable, pc_x = "PC1", pc_y = "PC2") {
  
  ggplot(
    pca_quality_df,
    aes(
      x = .data[[pc_x]],
      y = .data[[pc_y]],
      color = .data[[variable]],
      shape = IID
    )
  ) +
    geom_point(size = 3, stroke = 1.2) +
    geom_text_repel(
      aes(label = SampleID),
      size = 3.5,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = animal_shapes) +
    scale_color_viridis_c(
      option = "plasma",
      na.value = "grey70",
      name = variable
    ) +
    labs(
      title = paste("PCA - Nostress Blood -", variable),
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
      shape = "Animal ID"
    ) +
    theme_classic(base_size = 14)
}

# plot_pca_quality <- function(variable, pc_x = "PC1", pc_y = "PC2") {
#   
#   ggplot(
#     pca_quality_df,
#     aes(
#       x = .data[[pc_x]],
#       y = .data[[pc_y]]
#     )
#   ) +
#     
#     # POINTS = gradient color (RIN or ng.ul)
#     geom_point(
#       aes(
#         color = .data[[variable]],
#         shape = IID
#       ),
#       size = 3,
#       stroke = 1.2
#     ) +
#     scale_shape_manual(values = animal_shapes) +
#     scale_color_viridis_c(
#       option = "plasma",
#       na.value = "grey70",
#       name = variable
#     ) +
#     
#     # reset color scale for text
#     ggnewscale::new_scale_color() +
#     
#     # TEXT = treatment colors
#     geom_text_repel(
#       aes(
#         label = SampleID,
#         color = Treatment
#       ),
#       size = 3.5,
#       max.overlaps = Inf,
#       show.legend = FALSE
#     ) +
#     scale_color_manual(
#       values = treatment_colors
#     ) +
#     
#     labs(
#       title = paste("PCA - Nostress Blood -", variable),
#       x = paste0(
#         pc_x, " (",
#         round(variance_explained[as.integer(sub("PC", "", pc_x))], 2),
#         "%)"
#       ),
#       y = paste0(
#         pc_y, " (",
#         round(variance_explained[as.integer(sub("PC", "", pc_y))], 2),
#         "%)"
#       ),
#       shape = "Animal ID"
#     ) +
#     theme_classic(base_size = 14)
# }

pca_RIN <- plot_pca_quality("RIN")
pca_ngul <- plot_pca_quality("ng.ul")

print(pca_RIN)
print(pca_ngul)

ggsave("PCA_Nostress_Blood_RIN.pdf",pca_RIN,width = 11,height = 8)
ggsave("PCA_Nostress_Blood_ngul.pdf",pca_ngul,width = 11,height = 8)
