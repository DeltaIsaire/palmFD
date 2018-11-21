#####################################################
# Palm FD project: trait_filling_comparison test code
#####################################################
#
# This is a companion file listing the R code used to verify the code used for
# trait filling comparisons, and to explore/verify inputs and outputs.
# The code in this file is non-essential and provided only for reference and
# bug-fixing purposes.
# Input files:
#   none.
# Generated output files:
#   none.


source(file="2_trait_filling_comparison.R")


# ----------------
# Data preparation
# ----------------
# Find out how many species are not shared between filled trait matrices:
MultiCC(list(mean  = traits.mean$species,
             BHPMF = traits.BHPMF$species),
        presence=FALSE, value=FALSE)


