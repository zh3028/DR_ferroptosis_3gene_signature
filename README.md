# Data and Code for Two Related Studies on Diabetic Retinopathy

## Overview
This repository contains analysis scripts for two related studies on diabetic retinopathy (DR):

- **Study 1 (in submission):** Integrative single-cell and machine learning analysis identifies a ferroptosis-related three-gene diagnostic signature for diabetic retinopathy
- **Study 2 (in preparation):** Therapeutic mechanism of *Tetrastigma hemsleyanum* in diabetic retinopathy: targeting HMOX1/CTSD-mediated ferroptosis in pathogenic macrophages

## ScienceDB DOI
All processed data (including Seurat objects, machine learning models, expression matrices, and ROC objects) are available at ScienceDB: https://doi.org/10.57760/sciencedb.0134m

## Repository Contents

### Scripts
| File | Description |
| :--- | :--- |
| `step01_GSVA.R` | GSVA analysis for ferroptosis pathways |
| `step02_NicheNet.R` | NicheNet ligand-target interaction analysis |
| `step03_TF_activity.R` | DoRothEA + VIPER transcription factor activity analysis |
| `step04_metabolism.R` | GSVA analysis for metabolic pathways |

### Processed Data (CSV)
| File | Description |
| :--- | :--- |
| `step01_gsva_scores.csv` | GSVA scores for ferroptosis pathways across macrophage subpopulations |
| `step01_gsva_wilcox.csv` | Wilcoxon test results for ferroptosis pathway activity |
| `step02_NicheNet_top_ligands.csv` | Top 20 predicted ligands from SPP1+ macrophages to stromal cells |
| `step02_NicheNet_ligand_target_links.csv` | Ligand-target gene interactions |
| `step03_TF_activity.csv` | Transcription factor activity scores across macrophage subpopulations |
| `step03_top_TFs.csv` | Top 20 transcription factors in pathogenic subpopulations |
| `step04_metabolism.csv` | Metabolic pathway activity across macrophage subpopulations |
| `step04_top_pathways.csv` | Top metabolic pathways in pathogenic subpopulations |

### Figures
| File | Description |
| :--- | :--- |
| `step03_TF_activity_heatmap.pdf` | Heatmap of transcription factor activity |
| `step04_metabolism_heatmap.pdf` | Heatmap of metabolic pathway activity |
| `Figure6_Combined_ABCDEF_222.png` | Combined Figure 6 (GSVA, NicheNet, VIPER, metabolism) |

### Session Info
- `sessionInfo.txt`: R session information for reproducibility

## Data Sources
- **GSE165784**: Human PDR fibrovascular membranes (scRNA-seq)
- **GSE209872**: Rat STZ-induced DR retina (scRNA-seq)
- **GSE102485**: Human PDR bulk RNA-seq (training cohort)
- **GSE160306**: Human DR bulk RNA-seq (external validation)
- **GSE60436**: Human PDR microarray (external validation)

## Key Findings

### Study 1
- Identified 11 macrophage subpopulations in human PDR
- SPP1+/HMOX1+ macrophages exhibit ferroptosis activation (HMOX1, CTSD, FTH1)
- Three-gene diagnostic signature (FTH1, LDHA, CTSZ) with 10-fold CV AUC=0.800 (training), 0.889 (GSE60436 after Z-score), 0.705 (FTH1 alone in GSE160306)
- FTH1 correlates with monocyte/macrophage infiltration (r=0.316, P=0.015)
- **GSVA analysis**: Ferroptosis_Mac exhibited simultaneous elevation of both driver and suppressor pathways (P = 5.82e-40)
- **NicheNet analysis**: SPP1+ macrophages communicate with fibroblasts through BMP2 and TGFB1, with endothelial cells through VEGFA and ANGPT1, and with pericytes through BMP2 and BMP4
- **DoRothEA + VIPER analysis**: MYC, SPI1, RFX5, HIF1A, ESR2, LYL1, NFKB1, and FOXL2 were top transcriptional regulators
- **Metabolic pathway analysis**: Ferroptosis_Mac exhibited the highest iron metabolism (+0.293) and pentose phosphate pathway (+0.223) activity

### Study 2
- CellChat analysis reveals SPP1_Mac communicates with pericytes (0.907), fibroblasts (0.819), and endothelial cells (0.658)
- Molecular docking shows active compounds of *Tetrastigma hemsleyanum* bind to HMOX1, CTSD, and SPP1
- Cross-species conservation of pathogenic signatures in rat STZ-DR model
- Wet-lab validation (in progress): HPLC quantification of ocular penetration, cell and animal experiments

## Requirements
- R v4.3.1
- Seurat v5.0.3, Harmony v1.2.0
- GSVA v2.2.1
- nichenetr, decoupleR, dorothea
- ggplot2, dplyr, tidyr

## Contact
Meng Zhang, Department of Clinical Laboratory, Tongde Hospital of Zhejiang Province
Email: zh3028@zuaa.zju.edu.cn

## License
This project is licensed under the MIT License.
