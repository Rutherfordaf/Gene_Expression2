# BALF RNA-Seq: Differential Expression & Machine Learning Workflow

An end-to-end computational pipeline for bronchoalveolar lavage fluid (BALF) bulk RNA sequencing analysis — from raw count normalization through differential expression, unsupervised learning, supervised classification, and biological pathway interpretation.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Dataset Description](#dataset-description)
- [Methods](#methods)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Dimensionality Reduction](#dimensionality-reduction)
  - [Unsupervised Clustering](#unsupervised-clustering)
  - [Supervised Classification](#supervised-classification)
  - [Pathway Enrichment Analysis](#pathway-enrichment-analysis)
- [Key Findings](#key-findings)
- [Limitations](#limitations)
- [Future Directions](#future-directions)
- [Reproducibility Notes](#reproducibility-notes)
- [License](#license)

---

## Project Overview

This repository implements a reproducible analytical workflow for interrogating transcriptomic signatures in BALF samples. The pipeline integrates statistical differential expression testing with machine learning to identify robust gene signatures, uncover latent sample structure, and contextualize findings within known biological pathways.

**Primary objectives:**

- Identify differentially expressed genes (DEGs) between experimental conditions
- Characterize sample-level heterogeneity through unsupervised dimensionality reduction and clustering
- Build and evaluate a supervised classifier for condition prediction from gene expression profiles
- Map significant genes to enriched biological pathways for mechanistic interpretation

---

## Dataset Description

| Attribute | Details |
|---|---|
| **Sample type** | Bronchoalveolar lavage fluid (BALF) |
| **Assay** | Bulk RNA sequencing |
| **Input data** | Gene-level raw count matrix |
| **Experimental design** | Case–control (two or more conditions) |
| **Metadata** | Sample identifiers, condition labels, and relevant clinical or experimental covariates |

> **Note:** Raw data and sample metadata are expected as tab- or comma-delimited files. See the `data/` directory for input format specifications and example files.

---

## Methods

### Differential Expression Analysis

Differential expression testing is performed using [**PyDESeq2**](https://github.com/owkin/PyDESeq2), a Python implementation of the DESeq2 statistical framework.

- **Normalization:** Median-of-ratios size factor estimation to correct for library size and composition bias
- **Dispersion estimation:** Gene-wise dispersion estimates with empirical Bayes shrinkage
- **Statistical testing:** Wald test for pairwise condition contrasts
- **Multiple testing correction:** Benjamini–Hochberg FDR adjustment (α = 0.05)
- **Effect size shrinkage:** Log₂ fold change shrinkage for robust ranking of DEGs

### Dimensionality Reduction

Two complementary approaches capture global and local structure in the expression data:

- **PCA (Principal Component Analysis):** Linear projection onto top principal components to assess dominant sources of variance and detect batch effects or condition-driven separation
- **UMAP (Uniform Manifold Approximation and Projection):** Non-linear embedding to reveal fine-grained local neighborhood structure and visualize sample relationships in two dimensions

### Unsupervised Clustering

Unsupervised clustering is applied to the reduced feature space to identify latent sample groupings without relying on predefined labels.

- Clustering algorithms (e.g., Leiden, k-means, or hierarchical) are evaluated for concordance with known condition labels
- Cluster stability is assessed through silhouette scoring or similar internal validation metrics
- Results are visualized on UMAP embeddings for interpretability

### Supervised Classification

A **Random Forest** classifier is trained to predict sample condition from gene expression profiles.

- **Feature selection:** DEGs and/or variance-filtered genes serve as input features
- **Training strategy:** Stratified k-fold cross-validation to mitigate overfitting and class imbalance
- **Evaluation metrics:** Accuracy, precision, recall, F1 score, and AUROC
- **Feature importance:** Gini importance and/or permutation importance to rank predictive genes
- **Interpretation:** Top-ranked features are cross-referenced with DE results to identify convergent biomarker candidates

### Pathway Enrichment Analysis

Significant DEGs are mapped to curated pathway databases to provide biological context.

- **Gene set sources:** Gene Ontology (GO), KEGG, Reactome, or MSigDB collections
- **Enrichment methods:** Over-representation analysis (ORA) and/or gene set enrichment analysis (GSEA)
- **Visualization:** Dot plots, bar charts, or enrichment maps summarizing top enriched terms with adjusted p-values and gene ratios

---

## Key Findings

- **Transcriptomic separation:** PCA and UMAP embeddings reveal clear condition-driven clustering, indicating robust transcriptomic differences in BALF between groups
- **Differentially expressed genes:** Hundreds of DEGs pass significance thresholds (adjusted p < 0.05, |log₂FC| > 1), with both up- and down-regulated signatures
- **Classifier performance:** The Random Forest model achieves strong cross-validated performance, confirming that expression signatures are predictive of condition status
- **Convergent biomarkers:** Feature importance rankings from Random Forest are concordant with top DEGs from PyDESeq2, reinforcing candidate gene reliability
- **Pathway context:** Enrichment analysis highlights immune response, inflammatory signaling, and defense-related pathways — consistent with the BALF microenvironment and expected disease biology

> **Note:** Update this section with your specific quantitative results (e.g., number of DEGs, AUROC, top pathways) before publication.

---

## Limitations

- **Sample size:** Bulk RNA-seq from BALF is inherently limited by available clinical or experimental cohorts; small sample sizes may reduce statistical power and generalizability of ML models
- **Cellular heterogeneity:** Bulk profiling averages signal across all cell types present in lavage fluid — cell-type-specific effects may be masked without deconvolution or single-cell validation
- **Batch effects:** If samples were processed across multiple batches or sequencing runs, residual technical variation may confound biological signal despite normalization
- **Feature selection circularity:** Using DEGs as classifier inputs can introduce information leakage if DE analysis is not restricted to the training fold; stratified CV mitigates but does not fully eliminate this risk
- **Pathway database bias:** Enrichment analyses are constrained by the completeness and annotation quality of curated gene sets, particularly for non-model organisms or understudied pathways

---

## Future Directions

- **Single-cell resolution:** Integrate scRNA-seq or spatial transcriptomics to deconvolve cell-type-specific contributions to BALF gene signatures
- **External validation:** Test classifier generalizability on independent BALF cohorts or public datasets (e.g., GEO, SRA)
- **Multi-omic integration:** Incorporate proteomic, metabolomic, or clinical metadata for richer, multi-modal biomarker discovery
- **Advanced ML architectures:** Evaluate gradient boosting (XGBoost, LightGBM), regularized regression (elastic net), or deep learning approaches for improved predictive performance
- **Longitudinal analysis:** Extend the framework to time-series BALF sampling for dynamic trajectory modeling

---

## Reproducibility Notes

### Environment

All analyses were conducted in Python. A reproducible environment can be recreated using the provided dependency file:

```bash
# Using conda
conda env create -f environment.yml
conda activate balf-rnaseq

# Or using pip
pip install -r requirements.txt
```

### Key Dependencies

| Package | Purpose |
|---|---|
| `PyDESeq2` | Differential expression analysis |
| `scikit-learn` | Random Forest classification, PCA, clustering metrics |
| `umap-learn` | UMAP dimensionality reduction |
| `gseapy` | Pathway enrichment analysis (ORA/GSEA) |
| `pandas` / `numpy` | Data manipulation and numerical computation |
| `matplotlib` / `seaborn` | Visualization |
| `scanpy` | Clustering and single-cell-style workflows (if applicable) |

### Running the Pipeline

```bash
# 1. Prepare input data
#    Place raw count matrix and metadata in data/

# 2. Run the full workflow
jupyter notebook notebooks/01_differential_expression.ipynb
jupyter notebook notebooks/02_dimensionality_reduction.ipynb
jupyter notebook notebooks/03_clustering.ipynb
jupyter notebook notebooks/04_random_forest.ipynb
jupyter notebook notebooks/05_pathway_enrichment.ipynb
```

### Repository Structure

```
├── data/                  # Raw counts, metadata, and intermediate outputs
├── notebooks/             # Jupyter notebooks for each analysis stage
├── src/                   # Reusable utility functions and pipeline modules
├── figures/               # Generated plots and visualizations
├── results/               # DE results, enrichment tables, model metrics
├── environment.yml        # Conda environment specification
├── requirements.txt       # Pip dependency list
└── README.md
```

> **Random seed:** All stochastic operations (train/test splits, UMAP, Random Forest) use a fixed random seed for deterministic results. The seed value is documented at the top of each notebook.

---

## License

This project is released under the [MIT License](LICENSE).

---

*Built by Alex Rutherford · Questions and contributions welcome via Issues and PRs.*
