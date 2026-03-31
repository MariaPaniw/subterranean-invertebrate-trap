# Subterranean Invertebrate Trap — Data & Code

> **Associated publication:**
> Jubber W, Fuller A, Manser MB, Paniw M. *Capturing the unseen: a low-cost method for stratified subterranean sampling of soil invertebrates in drylands.* **Ecological Solutions and Evidence**, in press.

---
> **Zenodo link:**
[![DOI](https://zenodo.org/badge/1037187717.svg)](https://doi.org/10.5281/zenodo.19300403)

---

## Overview

This repository contains the data and R analysis code used in the above paper. The study introduces and evaluates two novel, low-cost subterranean pitfall trap designs for sampling soil invertebrates in dryland ecosystems. Traps were deployed across **10 sites in the Kalahari** between **October 2024 and March 2025**, yielding records of invertebrate taxonomy, abundance, and dry biomass across two trap designs and up to three soil depth layers. The repository also includes a companion conventional pitfall trap dataset and a descriptive analysis script that directly compares taxonomic capture profiles between the two sampling methods.

---

## Repository structure

```
subterranean-invertebrate-trap/
│
├── TableS1_forAnalysis.csv                    # Main subterranean invertebrate occurrence data
├── Pitfall_data2.csv                          # Conventional pitfall trap occurrence data
├── subterraneanHabitat_202502.csv             # Site-level habitat covariates
├── biomass_data/                              # Per-collection biomass Excel files
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
├── analyses_subterranean_submission.R         # Statistical analyses (GLMMs, MCMCglmm)
├── Subterranean_trap_descriptive_analyses.R   # Descriptive & comparative visualisations
├── LICENSE                                    # CC0 1.0 Public Domain Dedication
└── README.md                                  # This file
```

---

## Data files

### `TableS1_forAnalysis.csv`

The primary subterranean trap occurrence dataset. Each row is one taxonomic record from one trap on one collection date.

| Column | Description |
|---|---|
| `Ref` | Trap reference code. Encodes site (chars 2–4), trap line (chars 5–7: `L01` = double-stratified, `L05` = three-stratified), and soil layer (chars 8–10: `L01`–`L03`). |
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

### `Pitfall_data2.csv`

Conventional surface pitfall trap occurrence data, collected at the same Kalahari sites. Used in `Subterranean_trap_descriptive_analyses.R` to compare taxonomic capture profiles between sampling methods. Each row is one taxonomic record from one trap on one collection date (5,190 records spanning January 2024 – February 2025).

| Column | Description |
|---|---|
| `Timestamp` | Collection date and time (`YYYY/MM/DD HH:MM`). |
| `Ref` | Trap reference code encoding site and trap position (e.g., `K18L05`). |
| `Class` | Invertebrate class. |
| `Order` | Taxonomic order. |
| `Family` | Taxonomic family. |
| `Subfamily` | Taxonomic subfamily (where identified). |
| `Tribe` | Taxonomic tribe (where identified). |
| `Genera` | Genus (where identified). |
| `Species` | Species epithet (where identified). |
| `LifeStage` | Life stage of the captured individual (e.g., `Adult`, `Larva`). |
| `Method` | Always `"Pitfall"` in this file. |
| `Number_caught` | Count of individuals of this taxon in this trap at this collection. |

**Note:** The descriptive script filters to records from 2024 onwards to align with the subterranean trap deployment period. `LifeStage` is present in this file but not in `TableS1_forAnalysis.csv`; the descriptive script handles this by treating missing life-stage values as `"Unidentified"`.

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

Two subterranean trap designs were deployed at each of the 10 sites, alongside conventional surface pitfall traps:

| Method | Line / file | Soil layers sampled |
|---|---|---|
| Double-stratified subterranean | `L01` in `TableS1_forAnalysis.csv` | L01, L02 |
| Three-stratified subterranean | `L05` in `TableS1_forAnalysis.csv` | L01, L02, L03 |
| Conventional pitfall | `Pitfall_data2.csv` | Surface |

Subterranean trap collections were made approximately monthly, yielding five **sampling occasions**:

| Occasion | Months |
|---|---|
| 1 | October – November 2024 |
| 2 | December 2024 |
| 3 | January 2025 |
| 4 | February 2025 |
| 5 | March 2025 |

---

## Analysis scripts

### `analyses_subterranean_submission.R`

A self-contained R script (1,296 lines, last updated 2025-07-06) that reproduces all statistical analyses and figures reported in the paper. It is divided into the following sections:

1. **Data loading and preparation** — reads all CSVs and biomass Excel files, parses dates, extracts site/line/layer from reference codes.
2. **Species richness — seasonal variation** — zero-inflated Poisson and negative binomial GLMMs comparing richness across sampling occasions, land-use types, habitat, and trap design.
3. **Species richness — by trap layer** — separate GLMMs for double- and three-stratified designs to test for depth effects on richness.
4. **Biomass — seasonal variation** — Tweedie GLMMs of total biomass across the same candidate predictors.
5. **Biomass — by trap layer** — Tweedie GLMMs testing for depth effects on biomass within each design.
6. **Abundance — multivariate Bayesian models** — three parallel-chain MCMCglmm Poisson models for Scarabaeidae and Tenebrionidae jointly, with habitat, sampling occasion, and trap design as predictors.

#### Required R packages

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

#### How to run

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

### `Subterranean_trap_descriptive_analyses.R`

A companion R script for qualitative and descriptive comparison of the subterranean and conventional pitfall trap datasets. It reads both `TableS1_forAnalysis.csv` and `Pitfall_data2.csv`, combines them, and produces the following outputs:

1. **Qualitative statistics table** — summarises total counts and proportional contributions by taxon (Family: Genus), life stage, and sampling method; displayed as an interactive `DT` table with copy/CSV/Excel export buttons.
2. **Venn diagram** — shows the overlap in taxa (Family: Genus combinations) detected by subterranean vs. pitfall traps.
3. **Sankey diagram (by taxon)** — flow diagram linking each sampling method to its top captured taxa (top 8 per method plus top 6 shared taxa), with node labels showing counts and percentage contributions.
4. **Sankey diagram (by family and life stage)** — as above but with `TaxonID` defined as Family: LifeStage, to visualise differences in life-stage composition between methods.

#### Required R packages

```r
install.packages(c(
  "tidyverse",    # data manipulation and plotting (includes dplyr, ggplot2, tidyr, etc.)
  "lubridate",    # date parsing
  "ggpubr",       # publication-ready ggplot2 extensions
  "RColorBrewer", # colour palettes
  "ggvenn",       # Venn diagrams with ggplot2
  "patchwork",    # combining ggplot panels
  "networkD3",    # interactive network/Sankey diagrams
  "flextable",    # formatted tables
  "officer",      # export to Word/PowerPoint
  "ggiraph",      # interactive ggplot2 graphics
  "scales",       # axis and colour scaling
  "viridis",      # colour scales
  "glue",         # string interpolation
  "DT"            # interactive HTML data tables
))

# ggsankey is not on CRAN — install from GitHub:
install.packages("remotes")
remotes::install_github("davidsjoberg/ggsankey")
```

#### How to run

1. Ensure both `TableS1_forAnalysis.csv` and `Pitfall_data2.csv` are in your working directory.
2. Open `Subterranean_trap_descriptive_analyses.R` in R or RStudio.
3. Update the `ggsave()` file paths near the end of the script to point to your preferred output directory.
4. Run the script from top to bottom. The interactive `DT` table will render in the RStudio Viewer pane or a browser; Venn and Sankey figures are saved as PNG files.

---

## License

This dataset and code are released into the public domain under the **Creative Commons Zero (CC0 1.0) Universal** licence. You may copy, modify, and distribute this work without asking permission. See [LICENSE](LICENSE) or https://creativecommons.org/publicdomain/zero/1.0/ for details.

---

## Citation

If you use these data or code, please cite the associated paper:

> Jubber W, Fuller A, Manser MB, Paniw M. Capturing the unseen: a low-cost method for stratified subterranean sampling of soil invertebrates in drylands. *Ecological Solutions and Evidence*, in press.

---

## Contact

For questions about the data or methods, please contact wrjubber@gmail.com or m.paniw@gmail.com.
