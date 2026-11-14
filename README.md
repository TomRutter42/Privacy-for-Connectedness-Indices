# Replication Kit: "Differential Privacy for Network Connectedness Indices"

## Authors

Tom Rutter, Yuxin Liu, Amin Rahimian

## Overview

This repository contains the code and saved outputs for the paper "Differential Privacy for Network Connectedness Indices." It covers binary-label stochastic block model simulations, continuous-label graphon simulations, and applications to Karnataka village networks, Amazon co-purchases, and Twitch social networks.

## Directory Structure

```text
00_run_all.R                 Run selected analyses
replication_helpers.R        Shared privacy mechanisms, graph/data helpers, and configuration
binary_simulations.R         Binary connectedness simulations
continuous_simulations.R     Continuous regression and bottom-quartile friend-rank simulations
karnataka_villages_analysis.R Village connectedness, table, and figures
twitch_analysis.R            Within-language Twitch regression and figures
amazon_analysis.R            Amazon product-level regression and figures
Figures/                     Output PDF figures and existing PNG variants
Tables/                      Village summary table
DerivedData/                 Simulation results and empirical summaries
Logs/                        Clean-run logs, session details, checksums, and validation
```

Each analysis can run independently from the repository root. The shared helper implements the privacy mechanisms used across analyses.

## Data Requirements

Three external datasets must be obtained separately before running the code. For convenience, all three datasets have been bundled and made available in a [shared Dropbox folder](https://www.dropbox.com/scl/fo/0qliudow4iy3nsn5vnxgt/AOJyP-LzLSXCwg6BFKB6Lq0?rlkey=73480n7acexzb0kfw380dlzwh&st=5x0l3j6o&dl=0). Download the inputs and set `DP_DATA_DIR` to the directory containing them. Alternatively, download each dataset from its original source:

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
| File names | `com-amazon.ungraph.txt` (edge list), `amazon-meta.txt` ([separate product metadata download](https://snap.stanford.edu/data/amazon-meta.html)) |
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

1. Install the R packages listed below.
2. Download the inputs. The village files must retain their `datav4.0/datav4.0/Data/` subdirectories; the Amazon and Twitch files belong directly under `DP_DATA_DIR`.
3. From the repository root, run:

   ```sh
   DP_DATA_DIR="/path/to/NetworkPrivacyData" Rscript 00_run_all.R
   ```

   In an interactive R session, use:

   ```r
   Sys.setenv(DP_DATA_DIR = "/path/to/NetworkPrivacyData")
   source("00_run_all.R")
   ```

To run a single empirical analysis:

```sh
DP_DATA_DIR="/path/to/NetworkPrivacyData" Rscript twitch_analysis.R
```

By default, the scripts use the paper's experiment settings and write to the existing `Figures/`, `Tables/`, and `DerivedData/` directories. The master runner checks empirical inputs before starting the simulations.

The current saved outputs come from one full run completed on November 13, 2026, at 18:15 Pacific, using an empty `DP_OUTPUT_DIR`, all five analyses, and the default repetition counts. The continuous homophily sweep uses total epsilon 4 (label epsilon 2 and edge epsilon 2). All 41 outputs are installed, including `twitch_language_summary.csv`; the current scripts and outputs match the recorded SHA-256 checksums. See the [run log](Logs/revisions3_2026-11-13_174815/run.log), [run manifest](Logs/revisions3_2026-11-13_174815/run-manifest.json), and [validation results](Logs/revisions3_2026-11-13_174815/validation.json).

| Setting | Meaning | Default |
|---|---|---|
| `DP_DATA_DIR` | Directory containing the raw datasets | `.data` |
| `DP_OUTPUT_DIR` | Output directory | Repository root |
| `DP_TASKS` | Comma-separated analyses for the master runner | `binary,continuous,villages,twitch,amazon` |
| `DP_ITERATIONS` | Optional integer repetition override, at least 2 | Experiment-specific counts |

Use a separate `DP_OUTPUT_DIR` when overriding counts. The unit controlled by `DP_ITERATIONS` depends on the experiment:

| Experiment | What `DP_ITERATIONS` controls |
|---|---|
| Binary homophily, privacy-budget, and cell-composition sweeps | Independent graphs per parameter setting. Each graph retains 15 privacy draws for homophily and composition and 10 for the privacy-budget sweep. |
| Binary budget allocation | Privacy draws per allocation on one fixed graph per network size. |
| Binary network-size sweep | Independent graphs per network size, with one privacy draw per budget on each graph. |
| Continuous graphon budget allocation | Privacy draws per allocation on one fixed graph and label vector per network size. |
| Continuous graphon homophily sweep | Privacy draws on one fixed graph and label vector per homophily/size setting. |
| Village, Twitch, and Amazon analyses | Privacy draws per budget on each fixed empirical network. |

Seeds are fixed separately for each experiment, so results do not depend on which other analyses were selected. Each analysis script reruns its simulations and generates its outputs; it does not load saved simulation results.

## Software Requirements

The code was validated with R 4.6.1 and uses these packages:

| Package | Purpose |
|---|---|
| `data.table` | Data manipulation |
| `ggplot2` | Plotting |
| `igraph` | Graph simulation |
| `Matrix` | Sparse matrix operations |
| `haven` | Reading village Stata files |
| `dplyr` | Village data preparation |
| `scales` | Axis formatting |
| `ggrepel` | Twitch plot labels |

```r
install.packages(c("data.table", "ggplot2", "igraph", "Matrix",
                   "haven", "dplyr", "scales", "ggrepel"))
```
