# ========================================================
# 三叶青治疗DR单细胞全流程分析代码 (由云端自动生成备份)
# ========================================================

# 1. 加载包
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(CellChat)
library(slingshot)

# 2. 设置工作目录
setwd('~/DR_project')

# 3. 读取人源GSE165784数据（6个样本）
# 注意：前4个是tsv.gz，后2个是标准10x文件夹
expr_1 <- fread('GSM5049904_RRD-ERM1_matrix.tsv.gz', data.table = FALSE); rownames(expr_1) <- expr_1[,1]; expr_1 <- expr_1[,-1]
expr_2 <- fread('GSM5049905_PDR-ERM-51_matrix.tsv.gz', data.table = FALSE); rownames(expr_2) <- expr_2[,1]; expr_2 <- expr_2[,-1]
expr_3 <- fread('GSM5049906_PDR-ERM-59_matrix.tsv.gz', data.table = FALSE); rownames(expr_3) <- expr_3[,1]; expr_3 <- expr_3[,-1]
expr_4 <- fread('GSM5277737_F47-PDR-ERM_matrix.tsv.gz', data.table = FALSE); rownames(expr_4) <- expr_4[,1]; expr_4 <- expr_4[,-1]
expr_5 <- Read10X(data.dir = 'PDR-FM-0609_matrix_10X/')
expr_6 <- Read10X(data.dir = 'PDR-ERM-210630_matrix_10X/')

# 4. 构建并合并人源对象
seu1 <- CreateSeuratObject(counts = expr_1, project = 'PDR_1', min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts = expr_2, project = 'PDR_2', min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts = expr_3, project = 'PDR_3', min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts = expr_4, project = 'PDR_4', min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts = expr_5, project = 'PDR_5', min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts = expr_6, project = 'PDR_6', min.cells = 3, min.features = 200)
seu_combined <- merge(seu1, y = c(seu2, seu3, seu4, seu5, seu6), add.cell.ids = c('S1','S2','S3','S4','S5','S6'), project = 'PDR')

# 5. 人源数据质控、降维、聚类
seu_combined[['percent.mt']] <- PercentageFeatureSet(seu_combined, pattern = '^MT-')
seu_combined <- subset(seu_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
seu_combined <- NormalizeData(seu_combined)
seu_combined <- FindVariableFeatures(seu_combined, selection.method = 'vst', nfeatures = 2000)
seu_combined <- ScaleData(seu_combined)
seu_combined <- RunPCA(seu_combined, features = VariableFeatures(seu_combined))
seu_combined <- RunHarmony(seu_combined, group.by.vars = 'orig.ident')
seu_combined <- RunUMAP(seu_combined, reduction = 'harmony', dims = 1:30)
seu_combined <- FindNeighbors(seu_combined, reduction = 'harmony', dims = 1:30)
seu_combined <- FindClusters(seu_combined, resolution = 0.8)
saveRDS(seu_combined, file = 'seu_combined_scored.rds')

# 6. 人源细胞注释（16个群）
new.cluster.ids <- c('Macro_C3+', 'Macro_SPP1+', 'Macro_GPNMB+', 'Macro_MRC1+', 'Mono_FCN1+', 'Mono_EREG+', 'Fibroblast', 'Proliferating', 'T cell', 'DC', 'Endothelial', 'Pericyte', 'Macro_IFN+', 'Proliferating_2', 'Plasma cell', 'NK cell')
names(new.cluster.ids) <- levels(seu_combined)
seu_combined <- RenameIdents(seu_combined, new.cluster.ids)

# 7. 提取并深挖致病巨噬细胞
macro_subset <- subset(seu_combined, idents = c('Macro_C3+', 'Macro_SPP1+', 'Macro_GPNMB+', 'Macro_MRC1+', 'Mono_FCN1+', 'Mono_EREG+', 'Macro_IFN+'))
macro_subset <- NormalizeData(macro_subset)
macro_subset <- FindVariableFeatures(macro_subset, selection.method = 'vst', nfeatures = 2000)
macro_subset <- ScaleData(macro_subset)
macro_subset <- RunPCA(macro_subset, features = VariableFeatures(macro_subset))
macro_subset <- RunUMAP(macro_subset, dims = 1:20)
macro_subset <- FindNeighbors(macro_subset, dims = 1:20)
macro_subset <- FindClusters(macro_subset, resolution = 0.5)
saveRDS(macro_subset, file = 'macro_subset_annotated_final.rds')

# 8. 致病打分与差异分析
macro_subset <- AddModuleScore(macro_subset, features = list(c('LDHA', 'SLC2A1', 'HK2', 'PKM')), name = 'Glycolysis_Score')
macro_subset <- AddModuleScore(macro_subset, features = list(c('GPX4', 'SLC7A11', 'ACSL4', 'HMOX1', 'NCOA4')), name = 'Ferroptosis_Score')
macro_subset <- AddModuleScore(macro_subset, features = list(c('EP300', 'HDAC1', 'HDAC2')), name = 'PTM_Enzyme_Score')
patho_markers <- FindMarkers(macro_subset, ident.1 = 'Ferroptosis_Mac', ident.2 = 'Resident_Mac', only.pos = TRUE)
write.csv(patho_markers, file = 'pathogenic_vs_non_pathogenic_markers.csv')

# 9. 细胞通讯分析 (CellChat)
cellchat <- createCellChat(object = GetAssayData(seu_combined, assay = 'RNA', slot = 'data'), meta = data.frame(labels = Idents(seu_combined), row.names = names(Idents(seu_combined))), group.by = 'labels')
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = 'triMean')
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = 'cellchat_DR_results.rds')

# 10. 拟时序轨迹分析 (Slingshot)
emb <- Embeddings(macro_subset, reduction = 'umap')
clusters <- macro_subset$seurat_clusters
sds <- slingshot(data = emb, clusterLabels = clusters, start.clus = '0')
saveRDS(sds, file = 'macro_trajectory_slingshot.rds')

# 11. 读取大鼠GSE209872数据并合并
rat_seu1 <- CreateSeuratObject(counts = counts_1, project = 'WT1_0wk', min.cells = 3, min.features = 200)
rat_seu2 <- CreateSeuratObject(counts = counts_2, project = 'WT2_0wk', min.cells = 3, min.features = 200)
rat_seu3 <- CreateSeuratObject(counts = counts_3, project = 'DR1_2wk', min.cells = 3, min.features = 200)
rat_seu4 <- CreateSeuratObject(counts = counts_4, project = 'DR2_4wk', min.cells = 3, min.features = 200)
rat_seu5 <- CreateSeuratObject(counts = counts_5, project = 'DR3_8wk', min.cells = 3, min.features = 200)
rat_seu <- merge(rat_seu1, y = c(rat_seu2, rat_seu3, rat_seu4, rat_seu5), add.cell.ids = c('WT1','WT2','DR1','DR2','DR3'), project = 'Rat_DR')

# 12. 大鼠数据降维、聚类、跨物种验证
rat_seu[['percent.mt']] <- PercentageFeatureSet(rat_seu, pattern = '^Mt-')
rat_seu <- subset(rat_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
rat_seu <- NormalizeData(rat_seu)
rat_seu <- FindVariableFeatures(rat_seu, selection.method = 'vst', nfeatures = 2000)
rat_seu <- ScaleData(rat_seu)
rat_seu <- RunPCA(rat_seu, features = VariableFeatures(rat_seu))
rat_seu <- RunHarmony(rat_seu, group.by.vars = 'orig.ident')
rat_seu <- RunUMAP(rat_seu, reduction = 'harmony', dims = 1:30)
rat_seu <- FindNeighbors(rat_seu, reduction = 'harmony', dims = 1:30)
rat_seu <- FindClusters(rat_seu, resolution = 0.8)
saveRDS(rat_seu, file = 'rat_seu_processed.rds')

# 13. 大鼠致病打分验证
rat_core <- c('Hmox1', 'Ctsd', 'Ctsb', 'Fth1', 'Cd68', 'Ldha', 'Gpnmb')
rat_seu <- AddModuleScore(rat_seu, features = list(rat_core), name = 'Pathogenic_Score')
saveRDS(rat_seu, file = 'rat_seu_final_scored.rds')
