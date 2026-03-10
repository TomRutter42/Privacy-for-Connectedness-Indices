# Replication Kit: "Differential Privacy for Network Connectedness Indices"

## Authors

Tom Rutter, Yuxin Liu, Amin Rahimian

## Overview

This repository contains all code required to reproduce the figures and results in the paper "Differential Privacy for Network Connectedness Indices". 

The scripts cover:
- Theoretical simulations using stochastic block models (SBM) and graphon networks, with both binary and continuous node labels.
- Empirical applications to three real-world datasets: Indian village networks (Banerjee et al. 2013), the Amazon co-purchase graph (Yang and Leskovec 2012), and the Twitch social network (Rozemberczki and Sarkar 2021).

## Directory Structure

```
(root)/                        Main scripts and data files
  00_run_all.R                 Master script — sources all scripts in order
  binary_simulations.R         Binary label simulations (SBM, MSE vs homophily/epsilon/size)
  continuous_simulations.R     Continuous label simulations (graphon, allocation, regression)
  amazon_analysis.R
  twitch_analysis.R
  karnataka_villages_analysis.R

Figures/                       All output PDF and PNG figures
Tables/                        Reserved for any table outputs
DerivedData/                   Intermediate CSV results written by scripts
```

## Data Requirements

Three external datasets must be obtained separately before running the code. For convenience, all three datasets have been bundled and made available in a [shared Dropbox folder](https://www.dropbox.com/scl/fo/0qliudow4iy3nsn5vnxgt/AOJyP-LzLSXCwg6BFKB6Lq0?rlkey=73480n7acexzb0kfw380dlzwh&st=5x0l3j6o&dl=0). Download and point `ext_data_dir` to wherever you store it. Alternatively, download each dataset from its original source:

### 1. Banerjee et al. (2013) Indian Village Networks

| Field | Detail |
|---|---|
| Dataset name | datav4.0 |
| Source | [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21538) |
| Citation | Banerjee, A., Chandrasekhar, A. G., Duflo, E., & Jackson, M. O. (2013). The diffusion of microfinance. *Science*, 341(6144). |
| Used by | `karnataka_villages_analysis.R` |

### 2. Amazon Co-purchase Graph

| Field | Detail |
|---|---|
| File names | `com-amazon.ungraph.txt` (edge list), `amazon-meta.txt` (product metadata, same archive) |
| Source | [Stanford SNAP](https://snap.stanford.edu/data/com-Amazon.html) |
| Citation | Yang, J., & Leskovec, J. (2012). Defining and evaluating network communities based on ground-truth. *ICDM*. |
| Used by | `amazon_analysis.R` |

### 3. Twitch Gamers Dataset

| Field | Detail |
|---|---|
| File names | `large_twitch_edges.csv`, `large_twitch_features.csv` |
| Source | [Stanford SNAP](https://snap.stanford.edu/data/twitch_gamers.html) |
| Citation | Rozemberczki, B., & Sarkar, R. (2021). Twitch Gamers: a Dataset for Evaluating Proximity Preserving and Structural Role-based Node Embeddings. arXiv:2101.03091. |
| Used by | `twitch_analysis.R` |

## How to Reproduce

1. Install required R packages (see [Software Requirements](#software-requirements) below).
2. Download the three external datasets and set paths as described above.
3. From the repository root, run:
   ```r
   source("00_run_all.R")
   ```
   Alternatively, run individual scripts directly. Scripts that do not require external data (all simulation scripts) can be run immediately. Scripts requiring the Karnataka village data, Amazon graph, or Twitch data must have the relevant files in place first.

4. Output figures are written to `Figures/` and intermediate results to `DerivedData/`.

## Software Requirements

R (version 4.0 or later) with the following packages:

| Package | Purpose |
|---|---|
| `data.table` | fast data manipulation |
| `ggplot2` | plotting |
| `VGAM` | Laplace distribution (`rlaplace`) |
| `igraph` | graph construction and analysis |
| `Matrix` | sparse matrix operations |
| `parallel` | parallel simulation (`mclapply`) |
| `haven` | reading Stata `.dta` files (village data) |
| `dplyr` | data wrangling |
| `scales` | axis label formatting |

Install all packages with:
```r
install.packages(c("data.table", "ggplot2", "VGAM", "igraph", "Matrix",
                   "parallel", "haven", "dplyr", "scales", "ggrepel", "tidyr"))
```

## Computational Note

Some simulation scripts are computationally intensive. In particular:

- `binary_simulations.R` — MSE vs network size runs 10,000 iterations per combination; expect several hours on a laptop.
- `continuous_simulations.R` — allocation sweep scales up to N=500,000 nodes and uses `parallel::mclapply`; run on a multi-core machine.
- `twitch_analysis.R` — processes ~168,000 nodes and runs 1,000 DP simulation iterations per epsilon level.

For quick testing, reduce `NUM_ITERATIONS` in the relevant scripts.
