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

######################################################################################################################################################################################
######################################################################################################################################################################################
# Data Loading
######################################################################################################################################################################################
######################################################################################################################################################################################

# Loading the dataset containing DEGS
results_angus_blood_dds_df <- read.delim("results_df_blood_angus_genes.txt", header = TRUE, stringsAsFactors = FALSE)
significants_blood_genes <- read.delim("significant_blood_angus_genes.txt", header = TRUE, stringsAsFactors = FALSE)
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
#genes_to_check <- c("TULP2", "KDM4B", "B4GALT2", "PRR14L", "HDAC4", "ATG16L2", "PPARD", "TENM4", "WNK2", "FMO2", "STYXL1", "TRIO", "RAB27B", "DBP", "LZTR1", "SEC14L5", "GRIA1", "XPR1", "SIGIRR", "NR2E1", "CCL5", "WWOX", "ATF6B", "ADD1", "DDX23", "CAPZB", "DLG5", "SNRNP200", "RUNX1", "IHH", "FARSA", "DKKL1", "ACVR2B", "GDNF", "MARK2", "FAM131C", "FBLN2", "HSD17B14", "RABIF", "UBN1", "EIF2AK1", "ZNF133", "CH25H", "CSNK1D", "GJB3", "FBXL13", "TIMELESS", "TMEM159", "MEF2D", "MACROD1", "MTFMT", "F9", "CACNA1C", "SUGP1", "GAS6", "bta-mir-10b", "RNF26", "WNT5A", "ERI3", "EPHA2", "RALGAPA2", "ISCA2", "PHF19", "KIAA1324L", "CFAP61", "KRT72", "STK11IP", "RMDN3", "MRPL40", "EFCAB6", "MAPK8IP3", "CAD", "PIGZ", "AGAP3", "FGFR2", "FKBP10", "KCTD10", "HPS1", "GTPBP1", "INPP5E", "RANBP10", "NRG1", "GTF2IRD1", "SULT2B1", "GRM4", "PANK4", "PNOC", "C1QTNF6", "CACNA1G", "UBTD1")
#angus_methylation_genes <- c("TBC1D22A", "MYH9", "ERGIC1", "MED27", "FBXL19", "ATG7", "ACADL", "GALNT9", "SETD1A", "MFSD14B", "PHKA1", "PARN", "PLPP3", "PARVB", "KIAA1324L", "PHLDB1", "CRB1", "CDH10", "MOB3B", "WWOX", "SMG6", "CRTC1", "MLH1", "UBXN6", "WDFY4", "CCDC36", "MMD2", "USP13", "ONECUT1", "POR", "THOC3", "MAP4K2", "PDZRN3", "EPB41L4B", "ZNF827", "MAP2K1", "SEPT14", "TMC4", "HUNK", "BCR", "MATN2", "INSC", "PDE4DIP", "PACRG", "TTC7B", "TDRD9", "IMPDH1", "ATP13A1", "CDK17", "MGLL", "U6", "LIMK2", "ANKRD29", "TRAPPC9", "PPARA", "HADH", "GAS7", "TEKT5", "BPI", "PTBP1", "FHAD1", "ERF", "B4GALT1", "MYBPC2", "PNMA1", "PSD4", "XPO4", "DLK1", "FMNL2", "POLR3B", "MPPED2", "SCARB1", "NRP1", "PGM1", "EVC2", "NBEA", "CRYL1", "CCDC85A", "CLPTM1", "RBFOX1", "ZFAND3", "OPTN", "CTNNA3", "DMTN", "GRHL2", "ASTE1", "GRIN2B", "KIF1C", "COL20A1", "SLC25A21", "MELK", "GSG1L", "SOX10", "KAT2B", "PPP1R10", "PNPLA7", "EML1", "ACVR1B", "SRC", "CCDC158", "DHX34", "MLLT6", "RGMA", "KLHL40", "TWIST2", "TMEM266", "GRK3", "ATP1A3", "GNA15", "KANSL2", "PHACTR1", "LHFPL1", "PTBP3", "ACOX2", "TNRC18", "DNM2", "KHNYN", "NLRX1", "ZFAND1", "EFTUD2", "NDRG1", "ZNF692", "YIF1B", "EZR", "CIB4", "SIRT7", "ZNF594", "AGPAT4", "PTTG1IP", "KIAA0586", "PTK2B", "CHMP6", "DAGLB", "EPHA2", "CSAD", "NDE1", "CSNK1G1", "ME3", "IL1R1", "DYNLT1", "ARHGAP35", "ABLIM2", "UBN2", "GRAMD1B", "RPRD1B", "ARHGAP25", "CTNNA1", "IGF2R", "GADD45A", "APBB1IP", "ITPRIPL2", "PID1", "GGCX", "LCT", "DLGAP4", "RAP1GAP", "RNF43", "STAC", "PROM1", "TMED3", "PTPRO", "CLCN2", "DRD3", "TECPR2", "ZNF532", "FAM196A", "RECK", "PCBP3", "KRT85", "CSK", "E2F1", "PNMA2", "L3MBTL2", "ADCK1", "TYROBP", "WBSCR17", "COL16A1", "HS3ST3A1", "TUBB1", "CCR4", "LIG1", "HDAC4", "POLA2", "SCUBE1", "REXO1", "MTNR1A", "FAM178B", "CHKA", "CGN", "KCNMA1", "NRP2", "DTX4", "SLC11A1", "KIRREL", "PERM1", "SLC25A37", "TBXA2R", "HTRA1", "LBH", "PIGF", "UHRF1BP1", "LRRC27", "HIP1R", "ATP11A", "NKX2-5", "AP2A2", "STUM", "FLYWCH2", "SLC43A2", "NAA25", "EVX1", "SYNE2", "CAMTA2", "IQCD", "BCL2", "SERGEF", "AIM1L", "TRPV6", "SEC63", "IP6K3", "SEC14L5", "CALD1", "ETS1", "TGFB3", "PARP11", "GRIP2", "GLDN", "KCND2", "ELAC2", "SLC45A1", "CPT1A", "ERBB3", "MAP3K14", "OSBP2", "SLC6A6", "C19orf68", "WNK2", "ACIN1", "UBAC1", "HLTF", "ABCA12", "PLCXD3", "MACROD1", "TCF4", "CTNNBL1", "TTYH2", "NFATC1", "HERC2", "ABCA13", "EVL", "EPHA1", "PAPSS2", "POLN", "CAPN5", "ZNF423", "NR2C2", "WWC1", "KLHL31", "GATB", "TOM1L2", "IFT88", "MAN1C1", "VWA5B1", "CC2D1A", "FNIP2", "TRANK1", "NEK7", "ACAD11", "ERI3", "TRAPPC5", "SMC1B", "COL5A1", "KCND3", "DCP1B", "RERE", "MLXIP", "LRP5", "CTIF", "LMX1A", "STX8", "CFAP77", "ITFG1", "EFCAB6", "HEG1", "ANAPC5", "CLIC5", "SERPINE1", "EXT2", "RHEB", "LRRC30", "MFSD4A", "UNC45A", "LONRF2", "PAX7", "EBF3", "LAMA1", "NPFF", "MYO1E", "WWC2", "PSAP", "FARS2", "CERK", "DNAH17", "PIGG", "MAP2K4", "LRP1", "MIPEP", "MINK1", "SNRPA", "EXOC3", "NFIA", "TMEM94", "TBC1D14", "RUNX3", "EPN2", "TXNRD1", "WDR1", "IGSF21", "EVI5L", "SIL1", "NFATC2", "TBCD", "DXO", "CELF2", "DUOX1", "TRMT1", "MAD1L1", "RNF157")

genes_found <- intersect(angus_methylation_genes, kegg_significant_genes)
genes_not_found <- setdiff(angus_methylation_genes, kegg_significant_genes)

