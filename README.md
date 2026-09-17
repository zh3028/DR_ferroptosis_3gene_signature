# Data and Code for "Ferroptosis-related three-gene diagnostic signature for diabetic retinopathy"

## ScienceDB DOI
https://doi.org/10.57760/sciencedb.0134m

## Study Title
Integrative single-cell and machine learning analysis identifies a 
ferroptosis-related three-gene diagnostic signature for diabetic retinopathy

## Directory Structure
- `01_scripts/`: R scripts used for all analyses
- `02_processed_data/`: Processed Seurat objects, ML models, expression matrices
- `03_figures/`: All main and supplementary figures
- `04_tables/`: Tables used in the manuscript
- `05_session_info/`: R session information for reproducibility

## Data Sources
- GSE165784: Human PDR fibrovascular membranes (scRNA-seq)
- GSE209872: Rat STZ-induced DR retina (scRNA-seq)
- GSE102485: Human PDR bulk RNA-seq (training cohort)
- GSE160306: Human DR bulk RNA-seq (external validation)
- GSE60436: Human PDR microarray (external validation)

## Key Findings
1. Identified 11 macrophage subpopulations in human PDR
2. SPP1+/HMOX1+ macrophages exhibit ferroptosis activation (HMOX1, CTSD, FTH1)
3. Three-gene diagnostic signature (FTH1, LDHA, CTSZ) with AUC=0.960 
   (training), 0.889 (GSE60436 after Z-score), 0.705 (FTH1 alone in GSE160306)
4. FTH1 correlates with monocyte/macrophage infiltration (r=0.316, P=0.015)

## Contact
Meng Zhang, Department of Clinical Laboratory, Tongde Hospital of Zhejiang Province
Email: zh3028@zuaa.zju.edu.cn