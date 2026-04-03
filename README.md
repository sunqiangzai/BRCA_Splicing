## Methods Summary

### Alternative Splicing Quantification

- **SpliceSeq** was used to calculate PSI values based on splice junction reads.
- **SUPPA2** was applied to transcript-level quantification data for cross-tool validation.

Both analyses were conducted following official documentation:

- SpliceSeq: Ryan et al., Nat Methods (2012)
- SUPPA2: Trincado et al., Genome Biology (2018)

### Differential Splicing Analysis

DEAS events were identified based on:

- ΔPSI threshold
- Statistical significance (as described in the manuscript)

### Clustering

- Consensus clustering was applied to PSI matrices
- Robustness assessed via resampling and ARI

### Classification

- XGBoost model for subtype prediction (TCGA internal validation)
- Centroid-based classification for cross-cohort projection

### Validation

- Independent in-house RNA-seq dataset
- Perturbation dataset (SF3A3 knockdown)
- RT-PCR validation for selected events

## Reproducibility

All scripts are provided with example workflows.  
Users should adjust file paths and input formats based on their local environment.

## Requirements

- R (≥ 4.1)
- Key packages:
  - data.table
  - ggplot2
  - survival
  - ComplexHeatmap
  - ConsensusClusterPlus
  - xgboost

## Contact

For questions or issues, please contact the authors.
