###################################
# Palm Functional Diversity Project
###################################
#
# This script executes all TEST code scripts in sensible order,
# thus providing a one-click way to re-generate all TEST output.
# This script does NOT run the MAIN scripts, only the scripts for the TEST
# analysis.
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

cat("Generating gap-filling scenarios:\n")
source(file = "test_1_trait_filling.R")

cat("Comparing trait-filling methods:\n")
source(file = "test_2_trait_filling_comparison.R")

cat("Testing accuracy of gap-filling:\n")
source(file = "test_3_gapfilling_accuracy.R")

cat("Identifying gap-filling error:\n")
source(file = "test_4_gapfilling_error.R")

cat("Calculating functional diversity indices:\n")
source(file = "test_5_FD_indices.R")

cat("Generating FD null model \"regional\":\n")
source(file = "test_6_FD_null_model_regional.R")

cat("Finished running all test scripts.\n")

