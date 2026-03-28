# Subterranean Invertebrate Trap — Data & Code

> **Associated publication:**
> Jubber W, Fuller A, Manser MB, Paniw M. *Capturing the unseen: a low-cost method for stratified subterranean sampling of soil invertebrates in drylands.* **Ecological Solutions and Evidence**, in press.

---
> **Zenodo link:**
[![DOI](https://zenodo.org/badge/1037187717.svg)](https://doi.org/10.5281/zenodo.19300403)

---

## Overview

This repository contains the data and R analysis code used in the above paper. The study introduces and evaluates two novel, low-cost subterranean pitfall trap designs for sampling soil invertebrates in dryland ecosystems. Traps were deployed across **10 sites in the Kalahari** between **October 2024 and March 2025**, yielding records of invertebrate taxonomy, abundance, and dry biomass across two trap designs and up to three soil depth layers.

---

## Repository structure

```
subterranean-invertebrate-trap/
│
├── TableS1_forAnalysis.csv              # Main invertebrate occurrence data
├── subterraneanHabitat_202502.csv       # Site-level habitat covariates
├── biomass_data/                        # Per-collection biomass Excel files
│   ├── *_25_10_2024.xlsx
│   ├── *_08_11_2024.xlsx
│   ├── *_21_12_2024.xlsx
│   ├── *_27_12_2024.xlsx
│   ├── *_28_01_2025.xlsx
│   ├── *_30_01_2025.xlsx
│   ├── *_20_02_2025.xlsx
│   ├── *_27_02_2025.xlsx
│   ├── *_26_03_2025.xlsx
│   └── *_31_03_2025.xlsx
│
├── analyses_subterranean_submission.R   # Full analysis script
├── LICENSE                              # CC0 1.0 Public Domain Dedication
└── README.md                           # This file
```

---

## Data files

### `TableS1_forAnalysis.csv`

The primary occurrence dataset. Each row is one taxonomic record from one trap on one collection date.

| Column | Description |
|---|---|
| `Ref` | Trap reference code. |
| `Timestamp` | Collection date and time (`DD/MM/YY HH:MM`). |
| `Method` | Trap design: `"double-stratified subterranean"` or `"three-stratified subterranean"`. |
| `Class` | Invertebrate class (e.g., `Insecta`). Rows with `"Unidentified"` are excluded from analyses. |
| `Order` | Taxonomic order. |
| `Family` | Taxonomic family. |
| `Subfamily` | Taxonomic subfamily (where identified). |
| `Tribe` | Taxonomic tribe (where identified). |
| `Genera` | Genus (where identified). |
| `Species` | Species epithet (where identified). |
| `Number_caught` | Count of individuals of this taxon in this trap at this collection. |

**Tip for reuse:** The `Ref` column encodes multiple variables. Use `substr(Ref, 2, 4)` for site, `substr(Ref, 5, 7)` for line, and `substr(Ref, 8, 10)` for soil layer — exactly as in the analysis script.

---

### `subterraneanHabitat_202502.csv`

Site- and line-level habitat covariates. Contains 20 rows (10 sites × 2 trap lines).

| Column | Description |
|---|---|
| `Site` | Composite site–line reference (same encoding as `Ref` above). |
| `land_use` | Macrohabitat / land-use category for the site. |
| `Sand` | Sand cover (%) at the site–line combination. Used as the continuous habitat covariate (`hab`) in models. |

---

### `biomass_data/` directory

A folder of Excel workbooks, one per collection event, named by date (`DD_MM_YYYY`). Each file has a single sheet with the columns:

| Column | Description |
|---|---|
| Column 1 | Collection date. |
| `Ref.` | Trap reference code (same encoding as above; 25 unique references per event). |
| `Weight` | Total dry biomass of all invertebrates in the trap (grams). |
| `Number` | Total individual count across all taxa. |

The analysis script reads and combines all files in this folder automatically.

---

## Sampling design

Two trap designs were deployed at each of the 10 sites:

| Line code | Design | Soil layers sampled |
|---|---|---|
| `L01` | Double-stratified | L01, L02 |
| `L05` | Three-stratified | L01, L02, L03 |

Collections were made approximately monthly, yielding five **sampling occasions**:

| Occasion | Months |
|---|---|
| 1 | October – November 2024 |
| 2 | December 2024 |
| 3 | January 2025 |
| 4 | February 2025 |
| 5 | March 2025 |

---

## Analysis script

### `analyses_subterranean_submission.R`

A self-contained R script (1,296 lines, last updated 2025-07-06) that reproduces all analyses and figures in the paper. It is divided into the following sections:

1. **Data loading and preparation** — reads all CSVs and biomass Excel files, parses dates, extracts site/line/layer from reference codes.
2. **Species richness — seasonal variation** — zero-inflated Poisson and negative binomial GLMMs comparing richness across sampling occasions, land-use types, habitat, and trap design.
3. **Species richness — by trap layer** — separate GLMMs for double- and three-stratified designs to test for depth effects on richness.
4. **Biomass — seasonal variation** — Tweedie GLMMs of total biomass across the same candidate predictors.
5. **Biomass — by trap layer** — Tweedie GLMMs testing for depth effects on biomass within each design.
6. **Abundance — multivariate Bayesian models** — three parallel-chain MCMCglmm Poisson models for Scarabaeidae and Tenebrionidae jointly, with habitat, sampling occasion, and trap design as predictors.

### Required R packages

```r
install.packages(c(
  "lubridate",   # date parsing
  "mgcv",        # GAMs (loaded but available for extension)
  "dplyr",       # data manipulation
  "ggplot2",     # plotting
  "viridis",     # colour scales
  "glmmTMB",     # GLMMs (Poisson, negative binomial, Tweedie)
  "readxl",      # reading biomass Excel files
  "MCMCglmm",    # Bayesian multivariate Poisson models
  "tidyr",       # data reshaping
  "hdrcde",      # highest-density region estimates
  "coda",        # MCMC diagnostics
  "MCMCvis",     # MCMC visualisation
  "patchwork",   # combining ggplot panels
  "scales"       # axis scaling in ggplot2
))
```

### How to run

1. Clone or download this repository.
2. Open `analyses_subterranean_submission.R` in R or RStudio.
3. Set your local working directory on the line:
   ```r
   setwd("YOUR_LOCAL_PATH")
   ```
   This should point to the root of the repository (where the CSV files live).
4. Update the biomass path:
   ```r
   path = "YOUR_PATH/biomass_data"
   ```
5. Run the script from top to bottom. Figures are saved as PDF files in the working directory.

> **Note:** The abundance analyses (Section 6) run three MCMC chains of 700,000 iterations each and may take several hours on a standard laptop. Convergence is assessed automatically via Gelman–Rubin diagnostics printed to the console.

---

## License

This dataset and code are released into the public domain under the **Creative Commons Zero (CC0 1.0) Universal** licence. You may copy, modify, and distribute this work without asking permission. See [LICENSE](LICENSE) or https://creativecommons.org/publicdomain/zero/1.0/ for details.

---

## Citation

If you use these data or code, please cite the associated paper:

> Jubber W, Fuller A, Manser MB, Paniw M. Capturing the unseen: a low-cost method for stratified subterranean sampling of soil invertebrates in drylands. *Ecological Solutions and Evidence*, in press.

---

## Contact

For questions about the data or methods, please contact wrjubber@gmail.com  or m.paniw@gmail.com.
