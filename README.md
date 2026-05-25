# Proteomics Analysis of the 3×Tg-AD Alzheimer's Mouse Model

## Publication

> **Early Oral Administration of D-Chiro-Inositol Reverses Hippocampal Insulin and Glutamate Signaling Deficits in the 3×Tg Humanized Mouse Model of Alzheimer's Disease**
>
> Pacheco-Sánchez B., Verheul-Campos J., Vargas A., Tovar R., Rodríguez-Pozo M., Navarro J.A., López-Gambero A.J., Baixeras E., Serrano-Castro P.J., Suárez J., Sanjuan C., Rivera P.\*, Rodríguez de Fonseca F.\*
>
> *Nutrients* 2025, **17**, 3024. https://doi.org/10.3390/nu17183024

---

## Project Overview

This repository contains the R scripts used for the **shotgun proteomics and phosphoproteomics analysis** of the dorsal hippocampus in the 3×Tg-AD triple transgenic mouse model of Alzheimer's disease (AD).

The proteomics data were generated to characterize the molecular phenotype of 12-month-old 3×Tg-AD mice and served as the basis for selecting early intervention targets (insulin signaling and glutamatergic synaptic markers) subsequently validated in 6-month-old animals treated with **D-chiro-inositol (DCI, Gyneos®)**.

The pipeline covers the full workflow from raw LFQ (Label-Free Quantification) abundance data exported from Proteome Discoverer to statistical analysis, visualization, and protein network/enrichment analysis.

---

## Biological Context

The **3×Tg-AD** mouse carries three human familial AD mutations (*APP Swedish*, *MAPT P301L*, *PSEN1 M146V*) and progressively develops amyloid-β plaques and neurofibrillary tau tangles. Proteomics of 12-month-old hippocampus identified:

- Dysregulation of **mitochondrial oxidative phosphorylation** (Cyc1, Cycs, Sdhb, Cox5a, Uqcrq, Ndufa7, Atp5d)
- Disruption of **glutamatergic synaptic proteins** involved in vesicle trafficking, cytoskeletal dynamics, and local translation
- Upregulation of **AKT2**, indicating insulin signaling pathway disruption
- Phosphoproteomic alterations in **AMPA receptor trafficking** (Dlg4, Epb41l1, Camk4, Prkcb) and presynaptic release machinery

These findings motivated the DCI intervention study reported in the paper.

---

## Animals

| Cohort | Age | n | Purpose |
|--------|-----|---|---------|
| Batch 1 | 10–12 months | 5F + 5M WT; 5F + 4M 3×Tg-AD | Shotgun proteomics & phosphoproteomics |
| Batch 2 | 6 months | 12 WT + 12 3×Tg-AD (mixed sex) | Plasma DCI monitoring |
| Batch 3 | 3 → 6 months | 12F + 12M WT; 12F + 12M 3×Tg-AD | DCI treatment (3 months) |

All experiments were approved by the Local Ethical Committee for Animal Research of the University of Málaga (CEUMA No. 203-2023-A, 1 April 2024), and conducted in compliance with EU Directive 2010/63/EU and Spanish Royal Decree 53/2013.

---

## Repository Structure

```
proteomics-3xTg-AD/
│
├── R/
│   ├── preprocessing/
│   │   └── 01_preprocessing.R            # Data loading, cleaning, imputation & normalization
│   ├── analysis/
│   │   └── 02_differential_expression.R  # Statistical tests, log2FC, volcano table
│   └── visualization/
│       └── 03_visualization.R            # Volcano, heatmap, PCA, network, enrichment plots
│
├── data/
│   ├── raw/                              # Input files (not tracked — see below)
│   └── processed/                        # Intermediate outputs
│
├── results/
│   ├── figures/                          # All generated figures (tif, svg, pdf, png)
│   └── tables/                           # Differentially expressed proteins, full stats table
│
├── docs/
│   └── pipeline_overview.md              # Step-by-step methodology description
│
├── .gitignore
└── README.md
```

---

## Analysis Pipeline

### 1. Preprocessing (`01_preprocessing.R`)
- Loads raw LFQ abundance data from Proteome Discoverer output
- Fills missing gene symbols using a manually curated dictionary
- Separates samples by genotype (WT / 3×Tg) and sex (male / female)
- **Outlier removal** per protein per group using a 3×IQR rule
- Filters proteins requiring ≥3 valid values per group
- Retains only proteins detected across **all four** sex × genotype groups
- **Median imputation** within groups
- **Normalization** by a down-scaling correction factor (column mean / minimum column mean), applied independently per sex group

### 2. Differential Expression (`02_differential_expression.R`)
- Compares 3×Tg vs WT across all animals
- Normality assessed per protein with **Shapiro-Wilk test**:
  - Both groups normal → **Welch's t-test**
  - At least one non-normal → **Wilcoxon rank-sum test**
- Computes log₂(Fold Change) and −log₁₀(p-value)
- Classifies proteins as **up-regulated** or **down-regulated**
  - Thresholds: p < 0.05 and |log₂FC| > log₂(1.25) [≈ ±0.32 / −0.41]

### 3. Visualization (`03_visualization.R`)
- **PCA** of normalized abundances (colored by genotype, shaped by sex)
- **Volcano plot** with labeled top candidates
- **Heatmap** of differentially expressed proteins (log₁₀ abundances, Euclidean + complete linkage clustering)
- **Protein–protein interaction network** (STRING-derived edges, ggraph / Fruchterman–Reingold layout)
- **GO and KEGG enrichment** bar plots (top 20 terms per category, FDR-colored)

All figures are saved in TIFF (600 dpi), SVG, PDF, and PNG formats.

---

## Input Data

| File | Location | Description |
|------|----------|-------------|
| `230414-LFQ_Filtradas.xlsx` | `data/raw/` | LFQ protein abundances from Proteome Discoverer (Mus musculus UniProt 2023-03-01) |
| `sample-metadata.xlsx` | `data/raw/` | Sample annotations: ID, genotype, sex, experiment |
| `Enrichment_Results_3xtg_005.xlsx` | `data/processed/` | GO/KEGG enrichment results from STRING |
| `RED_3xtg_005.xlsx` | `data/processed/` | Protein–protein interaction edges from STRING |

> ⚠️ Raw data files are not tracked in this repository. Raw data are available upon reasonable request from the corresponding authors: [riveragonzalezpatricia@uma.es](mailto:riveragonzalezpatricia@uma.es) or [fernando.rodriguez@ibima.eu](mailto:fernando.rodriguez@ibima.eu)

---

## Proteomics Methods Summary

- **Instrument:** Easy nLC 1200 UHPLC coupled to Q Exactive HF-X Hybrid Quadrupole-Orbitrap (Thermo Fisher Scientific)
- **Search engine:** Sequest HT via Proteome Discoverer 2.5; 1% FDR (Percolator)
- **Database:** *Mus musculus* UniProt (2023-03-01)
- **Precursor tolerance:** 10 ppm | **Fragment tolerance:** 0.02 Da
- **Fixed modification:** Carbamidomethylation (Cys +57.021 Da)
- **Variable modifications:** Methionine oxidation (+15.996 Da); Phosphorylation (Ser/Thr/Tyr, +79.966 Da)
- **Quantification:** Label-free, precursor ion intensities

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
# Set working directory to repo root
setwd("path/to/proteomics-3xTg-AD")

# Step 1 — Preprocessing
source("R/preprocessing/01_preprocessing.R")

# Step 2 — Differential expression
source("R/analysis/02_differential_expression.R")

# Step 3 — Visualization
source("R/visualization/03_visualization.R")
```

---

## Results Summary

| Comparison | Proteins detected | Up-regulated | Down-regulated |
|------------|:-----------------:|:------------:|:--------------:|
| 3×Tg vs WT (12 months) | 3,257 | 100 | 111 |
| Phosphoproteome (12 months) | 539 | 30 | 84 |

---

## Authors

| Role | Name | Affiliation |
|------|------|-------------|
| Proteomics analysis | **Julia Verheul-Campos** | IBIMA Plataforma BIONAND, Málaga |
| Animal experiments & sampling | Beatriz Pacheco-Sánchez, A.J. López-Gambero, J.A. Navarro, R. Tovar | IBIMA, Málaga |
| Biochemical analyses | Beatriz Pacheco-Sánchez, A. Vargas, J.A. Navarro, P. Rivera | IBIMA, Málaga |
| Study concept & design | Patricia Rivera, Fernando Rodríguez de Fonseca | IBIMA, Málaga |

**Correspondence:**
- Patricia Rivera — [riveragonzalezpatricia@uma.es](mailto:riveragonzalezpatricia@uma.es)
- Fernando Rodríguez de Fonseca — [fernando.rodriguez@ibima.eu](mailto:fernando.rodriguez@ibima.eu)

---

## Funding

This study was funded by:
- Instituto de Salud Carlos III (ISCIII): projects DTS22/00021 and PI22/00427 (co-funded by the European Union)
- Ministerio de Ciencia e Innovación — Agencia Estatal de Investigación: RTC2019-007329-1
- Consejería de Salud y Consumo, Junta de Andalucía: HOR2025-05 (ERA4Health, EU Horizon Europe GA No. 101095426)
- Consejería de Universidad, Investigación e Innovación, Junta de Andalucía: PI21/00291
- Grant RYC2023-044921-I (MCIU/AEI/10.13039/501100011033 and FSE+)

---

## Citation

If you use this code, please cite:

```bibtex
@article{PachecoSanchez2025,
  author    = {Pacheco-Sánchez, Beatriz and Verheul-Campos, Julia and Vargas, Antonio
               and Tovar, Rubén and Rodríguez-Pozo, Miguel and Navarro, Juan A.
               and López-Gambero, Antonio J. and Baixeras, Elena
               and Serrano-Castro, Pedro J. and Suárez, Juan and Sanjuan, Carlos
               and Rivera, Patricia and Rodríguez de Fonseca, Fernando},
  title     = {Early Oral Administration of {D}-Chiro-Inositol Reverses Hippocampal
               Insulin and Glutamate Signaling Deficits in the {3$\times$Tg} Humanized
               Mouse Model of {Alzheimer's} Disease},
  journal   = {Nutrients},
  year      = {2025},
  volume    = {17},
  pages     = {3024},
  doi       = {10.3390/nu17183024}
}
```

---

## License

This code is released under the **MIT License** — see `LICENSE` for details.

© 2025 Pacheco-Sánchez et al. Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (article); MIT (code).
