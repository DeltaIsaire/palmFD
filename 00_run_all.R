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

print("Loading functions...",quote=FALSE)
source(file="functions/base.functions.R")

print("Gap-filling palm traits matrix:",quote=FALSE)
source(file="1_trait_filling.R")

print("Comparing trait filling methods...",quote=FALSE)
source(file="2_trait_filling_comparison.R")

print("Done.",quote=FALSE)
