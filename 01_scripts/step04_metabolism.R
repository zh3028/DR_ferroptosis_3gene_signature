#!/usr/bin/env Rscript
# Step 04: 代谢通路分析（GSVA）
suppressPackageStartupMessages({
  library(Seurat); library(GSVA); library(GSEABase)
  library(dplyr); library(tidyr); library(ggplot2)
})
log_file <- "/home/zm3028/research_platform/ZJ1/DR_project_v2/05_logs/step04_metabolism_log.txt"
sink(log_file, append = TRUE, split = TRUE)
mac <- readRDS("/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step00_mac_with_celltype.rds")
metabolism_pathways <- list(
  Glycolysis = c("HK1", "HK2", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA", "ALDOB", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB"),
  Oxidative_Phosphorylation = c("NDUFA1", "NDUFA2", "NDUFB1", "NDUFB2", "SDHA", "SDHB", "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "ATP5A1", "ATP5B"),
  Fatty_Acid_Oxidation = c("CPT1A", "CPT1B", "CPT2", "ACADVL", "ACADM", "ACADS", "HADHA", "HADHB", "ACAA1", "ACAA2"),
  Glutamine_Metabolism = c("GLS", "GLUD1", "GLUD2", "GOT1", "GOT2", "GPT", "GPT2"),
  Pentose_Phosphate = c("G6PD", "PGD", "RPE", "RPIA", "TKT", "TALDO1"),
  Iron_Metabolism = c("FTH1", "FTL", "TFRC", "SLC40A1", "HMOX1", "HMOX2", "IREB2"),
  Lipid_Synthesis = c("FASN", "ACACA", "ACACB", "SCD", "ELOVL6", "DGAT1", "DGAT2"),
  Amino_Acid_Metabolism = c("GOT1", "GOT2", "GPT", "GPT2", "GLUD1", "GLUD2", "GLS")
)
metabolism_pathways_filtered <- lapply(metabolism_pathways, function(x) x[x %in% rownames(mac@assays$RNA@data)])
expr_matrix <- as.matrix(mac@assays$RNA@data)
gsva_result <- gsva(gsvaParam(exprData = expr_matrix, geneSets = metabolism_pathways_filtered, kcdf = "Gaussian"), verbose = TRUE)
saveRDS(gsva_result, "/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step04_metabolism_gsva.rds")
metabolism_scores <- as.data.frame(t(gsva_result))
metabolism_scores$Subpopulation <- mac$cell_type
metabolism_by_subtype <- metabolism_scores %>% group_by(Subpopulation) %>% summarise_all(mean)
write.csv(metabolism_by_subtype, "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step04_metabolism.csv", row.names = FALSE)
p_metabolism <- metabolism_by_subtype %>% pivot_longer(-Subpopulation, names_to = "Pathway", values_to = "Score") %>% ggplot(aes(x = Subpopulation, y = Pathway, fill = Score)) + geom_tile() + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 10)) + labs(title = "Metabolic pathway activity across macrophage subpopulations", x = "Subpopulation", y = "Metabolic pathway") + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("/home/zm3028/research_platform/ZJ1/DR_project_v2/03_figures/step04_metabolism_heatmap.pdf", p_metabolism, width = 12, height = 6)
pathogenic_subtypes <- c("Ferroptosis_Mac", "Pro_fibrotic_Mac", "Pro_fibrotic_Mac2")
metabolism_pathogenic <- metabolism_by_subtype %>% filter(Subpopulation %in% pathogenic_subtypes) %>% dplyr::select(-Subpopulation)
top_pathways <- sort(colMeans(metabolism_pathogenic), decreasing = TRUE)
write.csv(data.frame(Pathway = names(top_pathways), Activity = top_pathways), "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step04_top_pathways.csv", row.names = FALSE)
sink()
