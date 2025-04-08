# In Silico Purification for Methylation-Based Tumor Classification

This repository accompanies our manuscript, *"In silico purification improves DNA methylation-based classification rates of pediatric low-grade gliomas"*.

## Background
DNA methylation-based classification has transformed diagnostics for CNS tumors, but many samples remain non-classifiable due to low tumor purity and infiltration by non-malignant cells. This repository provides an R-based framework for in silico purification, enabling confident classification by computationally removing epigenetic signals from non-malignant cells (e.g., neurons, immune cells). The method systematically subtracts reference methylation signatures of sorted non-malignant cell types from tumor profiles at varying proportions. Each adjusted profile is classified individually, and the most confident result is retained.

While developed for pediatric low-grade gliomas and demonstrated here using the Heidelberg Classifier v11b4 (Capper et al., 2018), the approach is flexible and broadly applicable, regardless of tumor type, classifier, or reference signatures.


## Repository Overview & Usage Instructions
### Included Resources
- **Reference signatures** (`ref.csv`) contain median beta-values for non-malignant cell types (microglia, monocytes, neutrophils, T cells, and neurons), derived from in vitro enriched methylation data from post-mortem brain tissue (GSE234520, GSE191200; Hannon et al., 2024; de Witte et al., 2022) and peripheral blood (GSE122244).

- **Example data** (`idats_example.tar.gz`) includes IDAT files for three representative samples from the Sturm et al. dataset (GSE215240), provided to demonstrate and test the framework.

- **Color mappings** (`colors_pur.Rds`) define a consistent color scheme for visualization, compatible with the methylation class (family) names of Heidelberg classifier v11b4.

### Requirements
- **Classifier:** A classifier is required to run this framework. The workflow is set up for the Heidelberg classifier v11b4 but can be easily adapted to work with any other methylation-based tumor classifier.

- **Functions:** The workflow uses MNPpreprocessIllumina and MNPpredict_betas for preprocessing and prediction (Capper et al., 2018). If these functions are not available, equivalent alternatives (e.g., minfi::preprocessIllumina, custom predict functions) can be used.

### Set up & Installations
Install required packages if not already installed. 
```R
install.packages("BiocManager")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
install.packages("randomForest")
install.packages("glmnet")
```

Clone the repository.
```bash
git clone https://github.com/livjuergensen/InSilicoPurification
```

### Preprocessing
```markdown
1_preprocessing.R
```
Uses raw IDAT files and generates a beta-value matrix `beta_values.rds`. Reference signatures `ref.csv` provided in this repository were processed using the same pipeline.

**Recommendation:** Apply identical preprocessing to both tumor and reference samples.

### Purification & Classification
```markdown
2_purification.R
```
Subtracts non-malignant reference signatures in increments (default: 0 to 0.98 in steps of 0.02) and classifies each adjusted profile using a random forest model. Generates classification probabilities across all reference cell types and subtraction fractions for methylation classes (`prob.RDS`) and methylation class families (`prob_mcf.RDS`), along with summary tables reporting the most confident classification per sample (`max_pur.csv`) and per sample–cell type combination (`max_ct_pur.csv`). Configured for the Heidelberg Classifier v11b4; adapt as needed for other classifiers.

### Visualization
```markdown
3_visualization.R
```
Generates stacked bar plots of calibrated classification probabilities (y-axis) across subtraction fractions (x-axis), grouped by reference cell type.

### Customization
- **Reference signatures:** Expand or replace `ref.csv` with custom cell-type profiles.
- **Classifier:**  Any methylation-based classifier can be used. Replace the Heidelberg Classifier loading and probability calculation steps with your own routine.
- **Fraction grid:** Adjust `seq(0, 0.98, 0.02)`.


## Citation
