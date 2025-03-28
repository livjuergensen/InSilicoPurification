# In Silico Purification for Methylation-Based Tumor Classification

This is the accompanying repository for our (not yet published) manuscript, "In silico purification improves DNA methylation-based classification rates of pediatric low-grade gliomas". 

## Background
Methylation-based CNS tumor classification has greatly improved brain tumor diagnostics. However, many samples fail to reach the classifier’s confidence threshold, often due to low tumor purity and infiltration by non-malignant cells. This repository provides an R-based framework for performing **in silico purification** of tumor DNA methylation profiles to enable confident classification of previously non-classifiable tumor samples. Therefore, we systematically  subtract reference methylation profiles of sorted non-malignant cell types from tumor profiles across a range of proportions. All adjusted profiles are classified individually, and the most confident result is retained. While developed for pediatric low-grade gliomas and shown here using the Heidelberg Classifier v11b4 (Capper et al., 2018), the approach is independent of classifier model, non-malignant reference and tumor type, and may provide advantages for broader integration. 

## Repository Overview

01_preprocessing.R
Demonstrates preprocessing of Illumina methylation array data as performed in our study, illustrated using example data from GEO (GSExx; Sturm et al., 2023). The script can be easily adapted for your own data. We used MNPpreprocessIllumina (Capper et al., 2018); if unavailable, preprocessIllumina from the minfi package yields highly similar results. Identical preprocessing steps are recommended for tumor and reference datasets. The script outputs beta-value matrices (.rds files) ready for downstream analysis.

02_in_silico_purification.R
Implements the core in silico purification method. The script systematically subtracts reference methylation signatures of non-malignant cell types from tumor profiles in small increments. Each adjusted profile is classified individually (here demonstrated with the Heidelberg Classifier v11b4). While the purification framework is classifier-independent, you must supply your own classification model and may need to adjust the provided code accordingly. The script outputs classification probabilities (.RDS files) for visualization and CSV summaries indicating the optimal (most confident) classification for each sample.
