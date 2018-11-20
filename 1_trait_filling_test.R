##########################################
# Palm FD project: trait_filling test code
##########################################
#
# This is a companion file listing the R code used to verify the code used for
# trait filling, and to explore/verify inputs and outputs.
# The code in this file is non-essential and provided only for reference and
# bug-fixing purposes.
# Input files:
#   none.
# Generated output files:
#   none.


source(file="1_trait_filling.R")


# ------------------------------------------------
# Compare the list of palm species of each dataset
# ------------------------------------------------
# 2. cross-check the three lists of species
missing.table <- MultiCC(species, presence=FALSE, value=FALSE)
missing.list <- MultiCC(species, presence=FALSE, value=TRUE)
# palm.dist and trait.data are the same, but palm.tree differs from the others.

# Test to check that this is true:
new.species <- list(palm.dist = present.list$palm.dist$palm.tree,
                    palm.tree = present.list$palm.tree$palm.dist)
MultiCC(new.species, presence=FALSE, value=FALSE)
# Yup, nothing missing.


# ---------------------------------
# Gap-fill trait matrix using BHPMF
# ---------------------------------
# To speed up code development or bug-testing you can run BHPMF with a test
# dataset:
test.trait.matrix <- trait.matrix[1:100, ]
test.hierarchy.matrix <- hierarchy.matrix[1:100, ]
