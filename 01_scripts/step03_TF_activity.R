#!/usr/bin/env Rscript
# Step 03: DoRothEA + VIPER 转录因子活性分析
suppressPackageStartupMessages({
  library(Seurat); library(decoupleR); library(dorothea)
  library(dplyr); library(tidyr); library(ggplot2)
})
log_file <- "/home/zm3028/research_platform/ZJ1/DR_project_v2/05_logs/step03_TF_activity_log.txt"
sink(log_file, append = TRUE, split = TRUE)
mac <- readRDS("/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step00_mac_with_celltype.rds")
expr_matrix <- as.matrix(mac@assays$RNA@data)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
tf_activity <- decoupleR::run_viper(mat = expr_matrix, network = regulons, .source = "tf", .target = "target", .mor = "mor", minsize = 5)
saveRDS(tf_activity, "/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step03_tf_activity.rds")
tf_wide <- tf_activity %>% select(source, condition, score) %>% pivot_wider(id_cols = condition, names_from = source, values_from = score)
tf_by_subtype <- tf_wide %>% mutate(Subpopulation = mac$cell_type) %>% group_by(Subpopulation) %>% summarise_all(mean)
write.csv(tf_by_subtype, "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step03_TF_activity.csv", row.names = FALSE)
p_tf <- tf_by_subtype %>% pivot_longer(-Subpopulation, names_to = "TF", values_to = "Activity") %>% ggplot(aes(x = Subpopulation, y = TF, fill = Activity)) + geom_tile() + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8)) + labs(title = "Transcription factor activity across macrophage subpopulations", x = "Subpopulation", y = "Transcription factor") + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("/home/zm3028/research_platform/ZJ1/DR_project_v2/03_figures/step03_TF_activity_heatmap.pdf", p_tf, width = 14, height = 10)
pathogenic_subtypes <- c("Ferroptosis_Mac", "Pro_fibrotic_Mac", "Pro_fibrotic_Mac2")
tf_pathogenic <- tf_by_subtype %>% filter(Subpopulation %in% pathogenic_subtypes) %>% dplyr::select(-Subpopulation)
top_tfs <- sort(colMeans(tf_pathogenic), decreasing = TRUE)[1:20]
write.csv(data.frame(TF = names(top_tfs), Activity = top_tfs), "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step03_top_TFs.csv", row.names = FALSE)
sink()
