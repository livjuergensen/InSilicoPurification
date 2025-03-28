# purify
in silico purification
# 🧬 In Silico Purification of Tumor Methylation Profiles

This repository provides the code and framework for performing **in silico purification** of tumor DNA methylation profiles by computationally removing non-malignant cell-type signatures. The goal is to improve the classification of low-purity pediatric CNS tumors using the **Heidelberg Brain Tumor Classifier** (Capper et al., 2018).

---

## 📌 Overview

> 💡 This approach is **model-agnostic** and could be adapted for any methylation-based classifier.  
> In this project, we use the Heidelberg Classifier v11b4 to assess classification improvements following computational purification.

---

## 🚀 Workflow

1. **Input**: Preprocessed beta value matrix and reference cell-type methylation signatures
2. **Purification**: Subtraction of cell-type-specific methylation profiles across a range of proportions
3. **Classification**: Use of the Heidelberg Classifier to evaluate adjusted profiles
4. **Result Aggregation**: Selection of the most confident classification per sample or per sample–cell-type combination
5. **Visualization**: Barplots of classification scores across purification levels

---

## 📁 Repository Structure
