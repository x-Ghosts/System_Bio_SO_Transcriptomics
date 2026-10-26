setwd("/home/adam/Desktop/PhD_work_DATA_ANALYSIS_161026/PhD_work/Systems_Biology/SystemsBiology_Brasile/RNA_Seq/step_3_GO_analysis/blood/")
options(timeout = 1200)

BiocManager::install(c("clusterProfiler", "org.Bt.eg.db", "enrichplot"))



library(dplyr)
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Bt.eg.db)
library(enrichplot)
library(patchwork)
library(ggplot2)
library(pathview)
library(KEGGREST)


######################################################################################################################################################################################
######################################################################################################################################################################################
# Data Loading
######################################################################################################################################################################################
######################################################################################################################################################################################

# Loading the dataset containing DEGS
results_angus_blood_dds_df <- read.delim("results_df_blood_angus_genes.txt", header = TRUE, stringsAsFactors = FALSE)
significants_blood_genes <- read.delim("significant_blood_angus_genes.txt", header = TRUE, stringsAsFactors = FALSE)
results_angus_blood_dds_df$gene <- rownames(results_angus_blood_dds_df)
significants_blood_genes$gene <- rownames(significants_blood_genes)
significant_gene_ids <- significants_blood_genes$gene

#Extracting the annotation db of Bos Taurus - org.Bt.eg.db + Mapping
valid_symbols <- AnnotationDbi::keys(org.Bt.eg.db,keytype = "SYMBOL")
mapped_genes <- intersect(significant_gene_ids, valid_symbols)
unmapped_genes <- setdiff(significant_gene_ids, valid_symbols)

background_genes <- results_angus_blood_dds_df %>%
  filter(!is.na(padj)) %>%
  pull(gene)
background_genes <- unique(background_genes) #11068 num of genes

mapped_background <- intersect(background_genes, valid_symbols)


######################################################################################################################################################################################
######################################################################################################################################################################################
# GO Enrichment
######################################################################################################################################################################################
######################################################################################################################################################################################

# Biological Process
go_bp_angus <- enrichGO(
  gene          = mapped_genes,
  universe      = mapped_background,
  OrgDb         = org.Bt.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  readable      = TRUE
)
# Molecular Function
go_mf_angus <- enrichGO(
  gene          = mapped_genes,
  universe      = mapped_background,
  OrgDb         = org.Bt.eg.db,
  keyType       = "SYMBOL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  readable      = TRUE
)
# Cellular Component
go_cc_angus <- enrichGO(
  gene          = mapped_genes,
  universe      = mapped_background,
  OrgDb         = org.Bt.eg.db,
  keyType       = "SYMBOL",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  readable      = TRUE
)
go_bp_results <- as.data.frame(go_bp_angus)
go_mf_results <- as.data.frame(go_mf_angus)
go_cc_results <- as.data.frame(go_cc_angus)

head(go_bp_results)
head(go_mf_results)
head(go_cc_results)


######################################################################################################################################################################################
######################################################################################################################################################################################
# GO Plotting
######################################################################################################################################################################################
######################################################################################################################################################################################

plot_GO_bp_angus <- dotplot(go_bp_angus, showCategory = 10,title = "GO Biological Process - Angus Blood")
plot_GO_mf_angus <- dotplot(go_mf_angus, showCategory = 10,title = "GO Molecular Function - Angus Blood")
plot_GO_cc_angus <- dotplot(go_cc_angus, showCategory = 10,title = "GO Cellular Component - Angus Blood")

print(plot_GO_bp_angus)
print(plot_GO_mf_angus)
print(plot_GO_cc_angus)

plot_all_GO_angus_blood <- (plot_GO_bp_angus / plot_GO_mf_angus / plot_GO_cc_angus)
print(plot_all_GO_angus_blood)

ggsave("GO_BP_Angus_Blood.pdf", plot = plot_GO_bp_angus, width = 10, height = 6)
ggsave("GO_MF_Angus_Blood.pdf", plot = plot_GO_mf_angus, width = 10, height = 6)
ggsave("GO_CC_Angus_Blood.pdf", plot = plot_GO_cc_angus, width = 10, height = 6)
ggsave("GO_ALL_Angus_Blood.pdf", plot = plot_all_GO_angus_blood, width = 17, height = 12)


write.csv(go_bp_results,"GO_BP_Angus_Blood_Results.csv",row.names = FALSE)
write.csv(go_mf_results,"GO_MF_Angus_Blood_Results.csv",row.names = FALSE)
write.csv(go_cc_results,"GO_CC_Angus_Blood_Results.csv",row.names = FALSE)

GO_angus_results <- bind_rows(as.data.frame(go_bp_angus) %>% mutate(Ontology = "BP"), as.data.frame(go_mf_angus) %>% mutate(Ontology = "MF"))
GO_angus_results <- GO_angus_results %>%
  dplyr::select(Ontology, ID, Description, GeneRatio, FoldEnrichment, pvalue, p.adjust, geneID, Count) %>%
  relocate(Ontology, .after = Description)
write.csv(GO_angus_results,"GO_ALL_Angus_Blood_Results.csv",row.names = FALSE)

# Useful if adding the LOG2FC & FDR from DEG to GO genes
# GO_angus_genes <- GO_angus_results %>%
#   dplyr::select(Ontology, Description, geneID) %>%
#   separate_rows(geneID, sep = "/") %>%
#   left_join(
#     significants_blood_genes %>%
#       dplyr::select(gene, log2FoldChange, padj),
#     by = c("geneID" = "gene")
#   )
# 
# View(GO_angus_genes)

######################################################################################################################################################################################
######################################################################################################################################################################################
# KEGG (Pathways) Enrichment
######################################################################################################################################################################################
######################################################################################################################################################################################

kegg_angus <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Bt.eg.db
)

kegg_significant_genes <- unique(kegg_angus$ENTREZID[kegg_angus$SYMBOL %in% significant_gene_ids])
kegg_background_genes <- unique(kegg_angus$ENTREZID)

kegg_angus_jp <- enrichKEGG(
  gene          = kegg_significant_genes,
  universe      = kegg_background_genes,
  organism      = "bta",
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  minGSSize     = 10
)

kegg_angus_jp <- as.data.frame(kegg_angus_jp)

kegg_angus_jp$gene_SYMBOL <- vapply(
  strsplit(kegg_angus_jp$geneID, "/", fixed = TRUE),
  function(ids) {
    symbols <- kegg_angus$SYMBOL[
      match(ids, kegg_angus$ENTREZID)
    ]
    paste(symbols, collapse = "/")
  },
  character(1)
)


kegg_bubble <- ggplot(
  kegg_angus_jp,
  aes(
    x = RichFactor,
    y = reorder(Description, RichFactor),
    size = Count,
    color = p.adjust
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    name = "Adjusted p-value"
  ) +
  scale_size_continuous(range = c(5, 12)) +
  labs(
    title = "KEGG Pathway Enrichment - Angus Blood",
    subtitle = "EXP vs POST | Significant DEGs",
    x = "Rich Factor",
    y = NULL,
    size = "Gene Count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

print(kegg_bubble)

ggsave("KEGG_Angus_Blood_Bubble.pdf",plot = kegg_bubble,width = 11,height = 6)


#####################################################################
# KEGG Pathview - Native KEGG + Blue/Red log2FC
#####################################################################

#####################################################################
# Prepare log2FC values for KEGG Pathview
#####################################################################

# 1. Match log2FC values to ENTREZ IDs
gene_fc <- results_angus_blood_dds_df$log2FoldChange[
  match(kegg_angus$SYMBOL, results_angus_blood_dds_df$gene)
]

names(gene_fc) <- kegg_angus$ENTREZID

# 2. Remove missing values and duplicated ENTREZ IDs
gene_fc <- gene_fc[is.finite(gene_fc)]
gene_fc <- gene_fc[!duplicated(names(gene_fc))]

pv_out <- list()
kegg_node_tables <- list()

symbol_map <- setNames(
  kegg_angus$SYMBOL,
  kegg_angus$ENTREZID
)

for (i in seq_len(nrow(kegg_angus_jp))) {
  
  pathway_full_id <- kegg_angus_jp$ID[i]
  pathway_id <- sub("^bta", "", pathway_full_id)
  
  # Significant genes in this pathway only
  pathway_genes <- strsplit(
    kegg_angus_jp$geneID[i], "/"
  )[[1]]
  
  gene_fc_pathway <- gene_fc[pathway_genes]
  gene_fc_pathway <- gene_fc_pathway[
    is.finite(gene_fc_pathway)
  ]
  
  if (length(gene_fc_pathway) == 0) next
  
  suffix <- paste0("Angus_", pathway_full_id, "_CyanOrange")
  
  # Preserve original KEGG labels and shapes
  pv_out[[pathway_full_id]] <- pathview(
    gene.data = gene_fc_pathway,
    pathway.id = pathway_id,
    species = "bta",
    gene.idtype = "entrez",
    kegg.native = TRUE,
    same.layer = FALSE,
    map.symbol = FALSE,
    node.sum = "mean",
    limit = list(gene = 1, cpd = 1),
    low = list(gene = "cyan", cpd = "cyan"),
    mid = list(gene = "grey90", cpd = "grey90"),
    high = list(gene = "orange", cpd = "orange"),
    both.dirs = list(gene = TRUE, cpd = TRUE),
    na.col = "transparent",
    out.suffix = suffix, kegg.dir = "KEGG_pathways/angus/"
  )
  
  # Retrieve mapped KEGG nodes
  pv_nodes <- pv_out[[pathway_full_id]]$plot.data.gene
  
  # Identify exactly which significant genes mapped to each node
  mapped_ids <- lapply(pv_nodes$all.mapped, function(x) {
    
    if (is.na(x) || !nzchar(x)) return(character(0))
    
    ids <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    
    intersect(ids, names(gene_fc_pathway))
  })
  
  # Convert Entrez IDs to symbols
  symbol_labels <- vapply(mapped_ids, function(ids) {
    
    if (length(ids) == 0) return("")
    
    symbols <- unname(symbol_map[ids])
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    
    paste(unique(symbols), collapse = "\n")
    
  }, character(1))
  
  # Store mapping information
  kegg_node_tables[[pathway_full_id]] <- data.frame(
    KEGG_Node = rownames(pv_nodes),
    KEGG_KGML_Label = pv_nodes$labels,
    Significant_Symbols = gsub("\n", "/", symbol_labels),
    stringsAsFactors = FALSE
  )
  
  ###################################################################
  # Optional: Add significant gene lists ABOVE shared KEGG nodes
  ###################################################################
  
  shared_nodes <- which(lengths(mapped_ids) > 1)
  
  if (length(shared_nodes) > 0) {
    
    input_png <- paste0(
      pathway_full_id, ".", suffix, ".png"
    )
    
    output_png <- paste0(
      pathway_full_id, ".", suffix, "_Labeled.png"
    )
    
    img <- png::readPNG(input_png)
    
    h <- dim(img)[1]
    w <- dim(img)[2]
    
    png(output_png, width = w, height = h)
    
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(
      xlim = c(0, w),
      ylim = c(0, h),
      xaxs = "i",
      yaxs = "i"
    )
    
    rasterImage(img, 0, 0, w, h, interpolate = FALSE)
    
    # Add names above nodes, without replacing original KEGG labels
    text(
      x = pv_nodes$x[shared_nodes],
      y = h - pv_nodes$y[shared_nodes] +
        pv_nodes$height[shared_nodes] / 2 + 8,
      labels = symbol_labels[shared_nodes],
      cex = 0.55,
      adj = c(0.5, 0),
      col = "black",
      font = 2
    )
    
    dev.off()
  }
}

#####################################################################
# Replace Pathview Labels with My RNA-seq Gene Symbols
#####################################################################

for (pathway_full_id in names(pv_out)) {
  
  pv_nodes <- pv_out[[pathway_full_id]]$plot.data.gene
  
  # Significant genes belonging to this pathway
  pathway_genes <- strsplit(
    kegg_angus_jp$geneID[
      match(pathway_full_id, kegg_angus_jp$ID)
    ],
    "/"
  )[[1]]
  
  # Extract significant Entrez IDs mapped to each node
  mapped_ids <- lapply(pv_nodes$all.mapped, function(x) {
    
    if (is.na(x) || !nzchar(x)) return(character(0))
    
    ids <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    
    intersect(ids, pathway_genes)
  })
  
  # Convert Entrez IDs using my own dataset
  gene_labels <- vapply(mapped_ids, function(ids) {
    
    if (length(ids) == 0) return("")
    
    symbols <- kegg_angus$SYMBOL[
      match(ids, kegg_angus$ENTREZID)
    ]
    
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    
    paste(unique(symbols), collapse = "\n")
    
  }, character(1))
  
  # Keep only significant mapped nodes
  keep <- nzchar(gene_labels) & is.finite(pv_nodes$mol.data)
  
  pv_nodes <- pv_nodes[keep, , drop = FALSE]
  pv_nodes$labels <- gene_labels[keep]
  
  if (nrow(pv_nodes) == 0) next
  
  # Redraw using my own gene symbols
  pathview::keggview.native(
    plot.data.gene = pv_nodes, kegg.dir = "KEGG_pathways/angus/",
    cols.ts.gene = matrix(pv_nodes$mol.col, ncol = 1),
    node.data = pathview::node.info(
      paste0(pathway_full_id, ".xml")
    ),
    pathway.name = pathway_full_id,
    same.layer = FALSE,
    low = list(gene = "cyan", cpd = "cyan"),
    mid = list(gene = "grey90", cpd = "grey90"),
    high = list(gene = "orange", cpd = "orange"),
    limit = list(gene = 1, cpd = 1),
    out.suffix = paste0(
      "Angus_", pathway_full_id, "_GeneSymbols"
    )
  )
}

# Requires a fix
# #####################################################################
# # KEGG - Original Nodes vs Pathview Mapped Genes
# #####################################################################
# 
# pathview_mapping <- dplyr::bind_rows(
#   lapply(names(pv_out), function(pathway_id) {
#     
#     # Pathview mapped nodes
#     pv_nodes <- pv_out[[pathway_id]]$plot.data.gene
#     
#     # Original KEGG node information
#     kegg_nodes <- pathview::node.info(
#       paste0(pathway_id, ".xml"),
#       short.name = FALSE
#     )
#     
#     node_ids <- rownames(pv_nodes)
#     
#     # Build node mapping table
#     data.frame(
#       ID = pathway_id,
#       KEGG_Node = node_ids,
#       Original_KEGG_Label = unname(kegg_nodes$labels[node_ids]),
#       Original_Node_Genes = vapply(
#         node_ids,
#         function(id) {
#           paste(
#             sub("^bta:", "", kegg_nodes$kegg.names[[id]]),
#             collapse = "/"
#           )
#         },
#         character(1)
#       ),
#       Pathview_Label = pv_nodes$labels,
#       EntrezID = pv_nodes$all.mapped,
#       Node_log2FC = pv_nodes$mol.data,
#       stringsAsFactors = FALSE
#     ) %>%
#       dplyr::filter(!is.na(EntrezID), EntrezID != "") %>%
#       tidyr::separate_rows(EntrezID, sep = ",") %>%
#       dplyr::mutate(EntrezID = trimws(EntrezID))
#   })
# )
# 
# # Compare with the original enrichment genes
# kegg_mapping_comparison <- kegg_angus_pathway_genes %>%
#   dplyr::transmute(
#     ID,
#     Description,
#     EntrezID = as.character(geneID),
#     Gene_Symbol = SYMBOL
#   ) %>%
#   dplyr::left_join(
#     pathview_mapping,
#     by = c("ID", "EntrezID")
#   ) %>%
#   dplyr::mutate(
#     Mapping_Status = ifelse(
#       is.na(KEGG_Node),
#       "Not mapped",
#       "Mapped"
#     )
#   )
# 
# View(kegg_mapping_comparison)
# 
# 
