###################################
# Palm Functional Diversity Project
###################################
#
# This script executes all code scripts in sensible order,
# thus providing a one-click way to re-generate all output.
# Keep in mind that the total runtime of all scripts is considerable.

# SET YOUR WORKING DIRECTORY TO THE PALM_FD DIRECTORY BEFORE RUNNING THIS SCRIPT.
# File paths in all scripts are relative to the palm_FD/ directory.
# The scripts depend on the following subdirectories being present:
#   data/
#   functions/
#   output/
#   graphs/
# Only the functions/ directory is available from the Github repository.
# The data/ directory contain datasets I shouldn't publish,
# else my supervisors would kill me :/
# The output/ and graphs/ directories can be missing, in which case this
# script will create them.  
#
# See the individual scripts for comments and details of their operations.

cat("Checking directories...\n")
if (!dir.exists("data/")) {
  stop("missing 'data/' directory. The scripts cannot be run without input data.")
}
if (!dir.exists("functions/")) {
  stop("missing 'functions/' directory. The scripts use custom functions.")
}
if (!dir.exists("output/")) {
  cat("creating directory: output/\n")
  dir.create("output/")
}
if (!dir.exists("graphs/")) {
  cat("creating directory: graphs/\n")
  dir.create("graphs/")
}

cat("Gap-filling palm traits matrix:\n")
source(file = "01_trait_filling.R")

cat("Comparing genus-mean and BHPMF gap-filling:\n")
source(file = "02_trait_filling_comparison.R")

cat("Testing accuracy of gap-filling:\n")
source(file = "03_trait_filling_accuracy.R")

cat("Identifying gap-filling error:\n")
source(file = "04_trait_filling_error.R")

cat("Creating stochastic genus-level gapfilling series:\n")
source(file = "05_trait_filling_stochastic.R")

cat("Calculating functional diversity indices:\n")
source(file = "06_FD_observed.R")

cat("Standardizing FD with null model (global):\n")
source(file = "07_FD_null_model_global.R")

cat("Standardizing FD with null model (realm):\n")
source(file = "08_FD_null_model_regional.R")

cat("Standardizing FD with null model (dispersion field):\n")
source(file = "09_FD_null_model_local.R")

cat("Parsing and graphing standardized FD results:\n")
source(file = "10_FD_results_overview.R")

cat("Parsing and preparing environmental predictors for regression:\n")
source(file = "11_TDWG_environmental.R")


cat("EXTRA: Performing k-means clustering of TDWG3 units based on ADF:\n")
source(file = "adf_kmeans_clustering.R")







cat("Finished running all main scripts.\n")

