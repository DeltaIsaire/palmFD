###################################
# Palm Functional Diversity Project
###################################
#
# This script executes all main scripts in sensible order,
# thus providing a one-click way to re-generate all output.
# Keep in mind that the total runtime of all scripts is considerable (the null
# model procedure can take a few weeks on a dedicated server).

# SET YOUR WORKING DIRECTORY TO THE PALM_FD DIRECTORY BEFORE RUNNING THIS SCRIPT
# (i.e. the directory where this script is located)
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
source(file = "01_stochastic_trait_filling.R")

cat("Calculating uncorrected functional diversity indices:\n")
source(file = "02_uncorrected_functional_diversity.R")

cat("Standardizing FD with null model (Global):\n")
source(file = "03_FD_null_model_Global.R")

cat("Standardizing FD with null model (Realm):\n")
source(file = "04_FD_null_model_Realm.R")

cat("Standardizing FD with null model (ADF):\n")
source(file = "05_FD_null_model_ADF.R")

cat("Parsing and summarizing FD results:\n")
source(file = "06_FD_summary.R")

cat("Parsing and preparing environmental predictors for regression:\n")
source(file = "07_preparing_predictors.R")

cat("Performing single-predictor SAR-error modeling:\n")
source(file = "08_single_SAR_models.R")

cat("Performing multimodel averaging with OLS and SAR models:\n")
source(file = "09_multimodel_averaging.R")

cat("Creating graphs of results:\n")
source(file = "10_results_graphics.R")

cat("Finished running all main scripts.\n")

