# In Silico Purification for Methylation-Based Tumor Classification

This repository provides an R-based framework for performing **in silico purification** of tumor DNA methylation profiles by computationally removing non-malignant cell-type signatures. This approach is designed to enhance classification rates of methylation-based classifiers, such as the **Heidelberg Brain Tumor Classifier**, especially in samples with low tumor purity. While developed for pediatric low-grade gliomas, the approach is independent of classifier model, non-malignant reference and tumor type. 

Methylation-based classifiers rely on the tumor-specific signal within a sample. In low-purity tumors, this signal may be diluted by non-malignant cell contributions, leading to non-classifiable results. Here, we systematically subtract reference methylation profiles of sorted non-malignant cell types from tumor profiles across a range of proportions. All adjusted profiles are classified individually, and the most confident result is retained.

## Overview
This pipeline includes:
- Preprocessing of raw methylation data
- Fractional removal of cell-type-specific methylation signatures
- Reclassification of all adjusted profiles
- Summary outputs of optimal classification
- Visualization of classification confidence across purified profiles

## Installation

Clone the repository:

```bash
git clone https://github.com/livjuergensen/in-silico-purification.git
cd in-silico-purification
