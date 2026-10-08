#!/usr/bin/env Rscript
# Step 01: GSVA 铁死亡通路评分
suppressPackageStartupMessages({
  library(Seurat); library(GSVA); library(GSEABase)
  library(dplyr); library(tidyr); library(ggplot2)
})
log_file <- "/home/zm3028/research_platform/ZJ1/DR_project_v2/05_logs/step01_GSVA_log.txt"
sink(log_file, append = TRUE, split = TRUE)
mac <- readRDS("/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step00_mac_with_celltype.rds")
ferroptosis_genes <- list(
  Ferroptosis_Driver = c("ACSL4", "LPCAT3", "ALOX15", "ALOX12", "POR", "NCOA4",
                          "TFRC", "SLC39A14", "HMOX1", "ATP6V1G2", "MAP1LC3A",
                          "GABARAPL1", "PEX3", "VDAC2", "VDAC3", "RPL8", "CS",
                          "ACO1", "ATF3", "ATF4", "IREB2"),
  Ferroptosis_Suppressor = c("GPX4", "SLC7A11", "SLC3A2", "GCLC", "GCLM", "GSR",
                              "NFS1", "AIFM2", "DHODH", "GCH1", "FTH1", "FTL",
                              "SLC40A1", "HSPB1", "NQO1", "NFE2L2", "MT1G",
                              "FANCD2", "CDKN1A", "TP53")
)
ferroptosis_genes_filtered <- lapply(ferroptosis_genes, function(x) {
  x[x %in% rownames(mac@assays$RNA@data)]
})
expr_matrix <- as.matrix(mac@assays$RNA@data)
gsva_result <- gsva(gsvaParam(exprData = expr_matrix, geneSets = ferroptosis_genes_filtered, kcdf = "Gaussian"), verbose = TRUE)
saveRDS(gsva_result, "/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step01_gsva_result.rds")
gsva_scores <- as.data.frame(t(gsva_result))
gsva_scores$Subpopulation <- mac$cell_type
write.csv(gsva_scores, "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step01_gsva_scores.csv", row.names = FALSE)
gsva_long <- gsva_scores %>% pivot_longer(cols = c("Ferroptosis_Driver", "Ferroptosis_Suppressor"), names_to = "Pathway", values_to = "Score")
p_gsva <- ggplot(gsva_long, aes(x = Subpopulation, y = Score, fill = Pathway)) + geom_violin(scale = "width", trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Ferroptosis pathway activity across macrophage subpopulations", x = "Subpopulation", y = "GSVA score") + scale_fill_manual(values = c("Ferroptosis_Driver" = "#E41A1C", "Ferroptosis_Suppressor" = "#377EB8"))
ggsave("/home/zm3028/research_platform/ZJ1/DR_project_v2/03_figures/step01_GSVA_Ferroptosis.pdf", p_gsva, width = 12, height = 6)
pathogenic_subtypes <- c("Ferroptosis_Mac", "Pro_fibrotic_Mac", "Pro_fibrotic_Mac2")
non_pathogenic_subtypes <- c("Resident_Mac", "Antigen_Mac", "M2_Mac", "Lipid_Mac")
pathogenic_scores <- gsva_scores %>% filter(Subpopulation %in% pathogenic_subtypes) %>% pull(Ferroptosis_Driver)
non_pathogenic_scores <- gsva_scores %>% filter(Subpopulation %in% non_pathogenic_subtypes) %>% pull(Ferroptosis_Driver)
wilcox_result <- wilcox.test(pathogenic_scores, non_pathogenic_scores)
write.csv(data.frame(Comparison = "Pathogenic vs Non-pathogenic", W = wilcox_result$statistic, P_value = wilcox_result$p.value), "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step01_gsva_wilcox.csv", row.names = FALSE)
sink()
