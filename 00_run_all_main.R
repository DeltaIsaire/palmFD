###################################
# Palm Functional Diversity Project
###################################
#
# This script executes all main code scripts in sensible order,
# thus providing a one-click way to re-generate all main output.
# This script does NOT run the test scripts, only the scripts for the main
# analysis.

# SET YOUR WORKING DIRECTORY TO THE PALM_FD DIRECTORY BEFORE RUNNING THIS SCRIPT.
# File paths in all scripts are relative to the palm_FD/ directory.
# The scripts depend on the following subdirectories being present:
#   data/
#   functions/
#   output/
#   graphs/
# Only the functions/ directory is available from the repository.
# The other directories contain datasets I shouldn't publish,
# else my supervisors would kill me :/ 
#
# See the individual scripts for comments and details of their operations.

cat("Checking directories...\n")
if (!dir.exists("data/")) {
  stop("missing 'data/' directory. The scripts cannot be run without input data.")
}
if (!dir.exists("functions/")) {
  stop("missing 'functions/' directory.")
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
source(file = "1_trait_filling.R")

cat("Done.\n")

