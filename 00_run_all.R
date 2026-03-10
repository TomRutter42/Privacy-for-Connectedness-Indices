# ==============================================================================
# 00_run_all.R — This script runs the replication kit.
#
# Before running you should make sure all of the external data files are in place,
# if you want to re-run the empirical applications.
#
# 1. Banerjee et al. (2013) village network data (datav4.0)
#    Available at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21538
#
# 2. Amazon co-purchase graph: com-amazon.ungraph.txt (and amazon-meta.txt)
#    Available at: https://snap.stanford.edu/data/com-Amazon.html
#
# 3. Twitch Gamers dataset: large_twitch_edges.csv, large_twitch_features.csv
#    Available at: https://snap.stanford.edu/data/twitch_gamers.html
#
# All three datasets should be placed in the folder pointed to by ext_data_dir below.
# The Banerjee data should retain its original folder structure (datav4.0/datav4.0/Data/).
#
# If it helps, I have put all of the data you will need in this Dropbox folder:
# https://www.dropbox.com/scl/fo/0qliudow4iy3nsn5vnxgt/AOJyP-LzLSXCwg6BFKB6Lq0?rlkey=73480n7acexzb0kfw380dlzwh&st=5x0l3j6o&dl=0
# Download it and point ext_data_dir to wherever you store it.
#
# Output figures are written to Figures/, tables to Tables/ and derived data to DerivedData/.
# ==============================================================================

# --- External data directory --- change me!
ext_data_dir <- "C:/Users/tomru/Dropbox (Personal)/NetworkPrivacyData"

# --- Binary label simulations ---
source("binary_simulations.R")

# --- Continuous label simulations ---
source("continuous_simulations.R")

# --- Empirical applications (require external data) ---

## Karnataka villages
source("karnataka_villages_analysis.R")

## Twitch
source("twitch_analysis.R")

## Amazon
source("amazon_analysis.R")
