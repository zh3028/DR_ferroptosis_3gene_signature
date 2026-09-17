library(Seurat)
macro_subset <- readRDS('macro_subset_annotated_final.rds')
DimPlot(macro_subset, reduction = 'umap', label = TRUE)
ggsave('Figure2_Pathogenic_vs_NonPathogenic_Heatmap.pdf')
