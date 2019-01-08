###################################
# Palm Functional Diversity Project
###################################
#
# This script executes all MAIN code scripts in sensible order,
# thus providing a one-click way to re-generate all MAIN output.
# This script does NOT run the TEST scripts, only the scripts for the MAIN
# analysis.

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
source(file = "1_trait_filling.R")

cat("Finished running all main scripts.\n")

