# In Silico Purification for Methylation-Based Tumor Classification

This repository provides an R-based framework for performing **in silico purification** of tumor DNA methylation profiles by computationally removing non-malignant cell-type signatures. This approach is designed to enhance the performance of methylation-based classifiers, such as the **Heidelberg Brain Tumor Classifier**, especially in samples with low tumor purity.

While developed for pediatric brain tumors, the approach is broadly applicable and classifier-agnostic.

## Overview

Methylation-based classifiers rely on the tumor-specific signal within a sample. In low-purity tumors, this signal may be diluted by immune or stromal cell contributions, leading to ambiguous or non-classifiable results. Here, we systematically subtract reference methylation profiles of sorted non-malignant cell types from tumor profiles across a range of proportions. All adjusted profiles are classified individually, and the most confident result is retained.

This pipeline includes:
- Fractional subtraction of cell-type-specific methylation signatures
- Reclassification of all adjusted profiles
- Summary outputs per sample and per cell type
- Visualization of classification confidence across purifications

## Installation

Clone the repository:

```bash
git clone https://github.com/livjuergensen/in-silico-purification.git
cd in-silico-purification
