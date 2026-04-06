# BRCA_Splicing

This repository contains custom R scripts used for the analysis of alternative splicing in breast cancer (BRCA), including DEAS detection, subtype classification, robustness evaluation, and immune-related analyses.

---

## Overview

This project focuses on:

- Identification of differentially expressed alternative splicing (DEAS) events
- Characterization of splicing landscapes in breast cancer
- Construction of splicing-based molecular subtypes
- Evaluation of subtype robustness and reproducibility
- Development of classification models
- Integration of immune infiltration and genomic features

---

## Methods Summary

### Alternative Splicing Quantification

- **SpliceSeq** was used to calculate PSI values based on splice junction reads.
- **SUPPA2** was applied to transcript-level quantification data for cross-tool validation.

Both analyses were conducted following official documentation:

- SpliceSeq: Ryan et al., *Nat Methods* (2012)  
- SUPPA2: Trincado et al., *Genome Biology* (2018)

---

### Differential Splicing Analysis

DEAS events were identified based on:

- PSI filtering thresholds (e.g., PSI ≥ 0.05, detection ≥ 75%)
- Paired statistical testing (Mann–Whitney U test)
- Multiple testing correction (BH-adjusted P < 0.05)

---

### Clustering

- Consensus clustering was applied to PSI matrices
- Robustness was assessed via:
  - repeated subsampling
  - adjusted Rand index (ARI)
  - patient-level assignment stability

---

### Classification

- **XGBoost model** was used for subtype prediction (internal validation in TCGA)
- **Centroid-based classification** was applied for cross-cohort subtype projection

---

### Validation

- Independent in-house RNA-seq dataset
- External cohort validation (SCAN-B)
- Perturbation dataset (SF3A3 knockdown, GSE147505)
- RT-PCR/qPCR validation for selected events

---

## Data Availability

The datasets analyzed in this study are publicly available:

- **TCGA-BRCA**: Genomic Data Commons (GDC)  
  https://portal.gdc.cancer.gov/

- **SCAN-B cohort**: GEO accession **GSE96058**

- **SF3A3 perturbation dataset**: GEO accession **GSE147505**

- **In-house RNA-seq dataset**: Genome Sequence Archive for Human  
  accession **HRA008406**  
  https://ngdc.cncb.ac.cn/gsa-human/s/13vB8KaC

---

## Scripts

### `BRCA_cluster_robustness_and_model.R`
Main script for clustering robustness analysis and classification modeling, including internal validation and performance evaluation.

### `DEAS_SF_cor.R`
Performs Pearson correlation analysis between splicing factor expression and DEAS PSI values to identify significant regulatory associations.

### `Detect_DEAS.R`
Identifies differentially expressed alternative splicing (DEAS) events between tumor and normal samples using PSI matrices and statistical testing.

### `centroid_based.R`
Implements centroid-based classification for subtype assignment and cross-cohort validation.

### `filter_psi.R`
Filters PSI matrices based on detection rate and PSI thresholds to retain high-confidence splicing events.

### `immune_analysis.R`
Performs immune infiltration analysis, including ESTIMATE score calculation and integration of CIBERSORT results for subtype comparisons.

### `robustness_consensuscluster_final.R`
Performs consensus clustering with resampling to assess stability of subtype classification.

---

## Reproducibility

All scripts are provided with example workflows.  
Users should adjust file paths, parameters, and input formats according to their local environment.

Note:

- SpliceSeq and SUPPA2 analyses were conducted following official pipelines.
- CIBERSORT was performed using the web platform (LM22, 100 permutations), and downstream analyses are included here.

---

## Requirements

- R (≥ 4.1)

### Key R packages:

- data.table  
- ggplot2  
- survival  
- ComplexHeatmap  
- ConsensusClusterPlus  
- xgboost  

---

## Contact

**Qiang Sun, Ph.D.**  
Email: qsun95@zju.edu.cn

For questions or issues related to this repository, please feel free to contact.
