###################################
# Palm Functional Diversity Project
###################################
#
# This main file executes all code scripts in sensible order,
# thus providing a one-click way to re-generate all output.

# SET YOUR WORKING DIRECTORY TO THE PALM_FD DIRECTORY BEFORE RUNNING THE CODE.
# File paths in all scripts are relative to the palm_FD/ directory.
# Code depends on the following subdirectories being present:
# data/, functions/, output/, graphs/
# Only the functions/ directory is available from the repository.
# The other directories contain datasets I shouldn't publish,
# 	else my supervisors would kill me :/ 
#
# See the individual scripts for comments and details of operations.

cat("Loading functions...", "\n")
source(file="functions/base.functions.R")
source(file="functions/plotting.functions.R")

cat("Gap-filling palm traits matrix:", "\n")
source(file="1_trait_filling.R")

cat("Comparing trait filling methods:", "\n")
source(file="2_trait_filling_comparison.R")

cat("Done.", "\n")
