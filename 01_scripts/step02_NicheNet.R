#!/usr/bin/env Rscript
# Step 02: NicheNet 配体-受体-靶基因网络
suppressPackageStartupMessages({
  library(Seurat); library(nichenetr); library(dplyr); library(tidyr); library(ggplot2)
})
log_file <- "/home/zm3028/research_platform/ZJ1/DR_project_v2/05_logs/step02_NicheNet_log.txt"
sink(log_file, append = TRUE, split = TRUE)
seu <- readRDS("/home/zm3028/research_platform/ZJ1/DR_project/seu_combined_annotated.rds")
seu$cell_type <- as.character(Idents(seu))
ligand_target_matrix <- readRDS("/home/zm3028/research_platform/ZJ1/nichenet_db/ligand_target_matrix_nsga2r_final.rds")
lr_network <- readRDS("/home/zm3028/research_platform/ZJ1/nichenet_db/lr_network_human_21122021.rds")
sender <- "Macro_SPP1+"
receivers <- c("Fibroblast", "Endothelial", "Pericyte")
background_genes <- rownames(seu@assays$RNA@data)
all_ligand_activities <- list()
all_de_genes <- list()
for (recv in receivers) {
  recv_de <- FindMarkers(seu, ident.1 = recv, ident.2 = NULL)
  recv_de_genes <- rownames(recv_de)[recv_de$avg_log2FC > 0.5 & recv_de$p_val_adj < 0.05]
  all_de_genes[[recv]] <- recv_de_genes
  ligand_activities <- predict_ligand_activities(geneset = recv_de_genes, background_expressed_genes = background_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = lr_network$from)
  ligand_activities$receiver <- recv
  all_ligand_activities[[recv]] <- ligand_activities
}
combined_ligand_activities <- do.call(rbind, all_ligand_activities)
top_ligands <- combined_ligand_activities %>% group_by(receiver) %>% arrange(desc(aupr_corrected)) %>% distinct(test_ligand, .keep_all = TRUE) %>% slice_head(n = 20) %>% ungroup()
write.csv(top_ligands, "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step02_NicheNet_top_ligands_clean.csv", row.names = FALSE)
all_ligand_target_links <- list()
for (recv in receivers) {
  recv_top_ligands <- top_ligands %>% filter(receiver == recv) %>% pull(test_ligand)
  active_ligand_target_links <- get_weighted_ligand_target_links(ligand = recv_top_ligands, geneset = all_de_genes[[recv]], ligand_target_matrix = ligand_target_matrix, n = 250)
  active_ligand_target_links$receiver <- recv
  all_ligand_target_links[[recv]] <- active_ligand_target_links
}
combined_ligand_target_links <- do.call(rbind, all_ligand_target_links)
write.csv(combined_ligand_target_links, "/home/zm3028/research_platform/ZJ1/DR_project_v2/04_tables/step02_NicheNet_ligand_target_links_clean.csv", row.names = FALSE)
saveRDS(combined_ligand_activities, "/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step02_ligand_activities.rds")
saveRDS(combined_ligand_target_links, "/home/zm3028/research_platform/ZJ1/DR_project_v2/02_rds/step02_ligand_target_links.rds")
sink()
