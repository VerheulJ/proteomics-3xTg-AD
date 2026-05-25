# Proteomics Analysis of 3xTg-AD Alzheimer's Mouse Model

## Overview

This repository contains R scripts for the quantitative proteomics analysis of the **3xTg-AD** triple transgenic mouse model of Alzheimer's disease, comparing it against wild-type (WT) controls.

The pipeline covers the full workflow from raw LFQ (Label-Free Quantification) abundance data to statistical analysis, visualization, and network/enrichment analysis.

---

## Repository Structure

```
proteomics-3xTg-AD/
│
├── R/
│   ├── preprocessing/
│   │   └── 01_preprocessing.R       # Data loading, cleaning, imputation & normalization
│   ├── analysis/
│   │   └── 02_differential_expression.R  # Statistical tests & volcano plot data
│   └── visualization/
│       └── 03_visualization.R       # Volcano plots, heatmaps, network & enrichment plots
│
├── data/
│   ├── raw/                         # Input files (LFQ abundances, sample metadata)
│   └── processed/                   # Intermediate outputs (ajustado.xlsx, entera_3xtg.xlsx)
│
├── results/
│   ├── figures/                     # Volcano plots, heatmaps, network graphs, enrichment
│   └── tables/                      # Differentially expressed proteins, full volcano table
│
├── docs/
│   └── pipeline_overview.md         # Detailed description of each analysis step
│
├── .gitignore
└── README.md
```

---

## Analysis Pipeline

### 1. Preprocessing (`01_preprocessing.R`)
- Loads raw LFQ abundance data from Proteome Discoverer output
- Fills missing gene symbols using a manually curated dictionary
- Separates samples by genotype (WT / 3xTg) and sex (male / female)
- **Outlier removal** per protein per group using a 3×IQR rule
- Filters proteins with fewer than 3 valid values per group
- **Median imputation** within groups
- **Normalization** by a correction factor based on column sums

### 2. Differential Expression (`02_differential_expression.R`)
- Compares 3xTg vs WT across all animals
- Selects test per protein based on Shapiro-Wilk normality:
  - Normal → Welch's t-test
  - Non-normal → Wilcoxon rank-sum test
- Computes log₂(Fold Change) and –log₁₀(p-value)
- Classifies proteins as **up-regulated**, **down-regulated**, or unchanged
  - Thresholds: p < 0.05 and |log₂FC| > log₂(1.25)

### 3. Visualization (`03_visualization.R`)
- **Volcano plot** with labeled top candidates
- **Heatmap** of differentially expressed proteins (log₁₀ abundances, hierarchical clustering)
- **Protein–protein interaction network** (STRING-derived interactions, ggraph)
- **GO and KEGG enrichment** bar plots (top 20 terms per category)

---

## Input Data

| File | Description |
|------|-------------|
| `data/raw/230414-LFQ_Filtradas.xlsx` | LFQ protein abundances from Proteome Discoverer |
| `data/raw/sample-metadata.xlsx` | Sample metadata (ID, genotype, sex) |
| `data/processed/Enrichment_Results_3xtg_005.xlsx` | GO/KEGG enrichment results (from STRING) |
| `data/processed/RED_3xtg_005.xlsx` | Protein–protein interaction edges (from STRING) |

> ⚠️ Raw data files are not included in this repository due to size/privacy constraints. Contact the authors to request access.

---

## Dependencies

All analyses were performed in **R**. Required packages:

```r
# Data handling
library(readxl); library(openxlsx); library(dplyr); library(tidyverse); library(data.table)

# Statistics
library(stats); library(effsize); library(car)

# Visualization
library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer); library(gplots)

# Dimensionality reduction
library(FactoMineR); library(factoextra); library(missMDA); library(VIM)

# Network analysis
library(igraph); library(ggraph)

# Enrichment / NGS utilities
library(edgeR); library(tidytext); library(tm)
```

Install all at once:
```r
install.packages(c("readxl", "openxlsx", "dplyr", "tidyverse", "data.table",
                   "effsize", "car", "ggplot2", "ggrepel", "pheatmap",
                   "RColorBrewer", "gplots", "FactoMineR", "factoextra",
                   "missMDA", "VIM", "igraph", "ggraph", "edgeR",
                   "tidytext", "tm", "formattable"))
```

---

## Usage

```r
# Step 1 — Preprocessing
source("R/preprocessing/01_preprocessing.R")

# Step 2 — Differential expression
source("R/analysis/02_differential_expression.R")

# Step 3 — Visualization
source("R/visualization/03_visualization.R")
```

Set your working directory to the root of this repository before running.

---

## Results Summary

| Comparison | Up-regulated | Down-regulated |
|------------|-------------|----------------|
| 3xTg vs WT | — | — |

> Fill in after running the analysis.

---

## Citation

If you use this code, please cite the corresponding publication (in preparation).

---

## License

MIT License — see `LICENSE` for details.
