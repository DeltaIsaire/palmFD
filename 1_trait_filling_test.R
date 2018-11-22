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

# It is good to test how BHPMF works with different input trait matrices,
# Especially because the function doesn't like it when a species has NA
# for all provided traits.

# standard three trait matrix, including all-3 NAs, with substituted zero values;
# and associated hierarchy matrix:
trait.matrix <- as.matrix(palm.traits[, c("stem.height", "blade.length",
                                          "fruit.length")
                                      ]
                          )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
trait.matrix[, "stem.height"][which(trait.matrix[, "stem.height"] == 0)] <- 0.0001

hierarchy.matrix <- as.matrix(data.frame(species   = sub(pattern=" ",
                                                         replacement="_",
                                                         x=trait.data$SpecName),
                                         genus     = trait.data$accGenus,
                                         tribe     = trait.data$PalmTribe,
                                         subfamily = trait.data$PalmSubfamily)
                              )

# Optional subsetting to cases without all NAs:
# First extract these cases to a seperate dataframe:
to.remove <- which(rowSums(is.na(trait.matrix)) == ncol(trait.matrix))
BHPMF.missing <- palm.traits[to.remove, ]
# Then remove them from the matrices:
trait.matrix <- trait.matrix[-to.remove, ]
hierarchy.matrix <- hierarchy.matrix[-to.remove, ]
rm(to.remove)
# Result: BHPMF works.

# Optional: add a 4th dummy trait with the same value (1) for all species.
# This could be another way to deal with all NA species.
trait.matrix <- as.matrix(data.frame(trait.matrix,
                                     dummy=rep(1, times=length(trait.matrix[, 1]))
                          )          )
# Result: BHPMF works, even including species with NA for the three real traits.
# TODO: However, we have to check the results for accuracy.

# Optional: remove all but the first column from the trait matrix:
trait.matrix <- as.matrix(trait.matrix[, "stem.height"])
# Result: gives the same error "there are observations with all features missing"


# BHPMF GapFilling function call, with subset for quick processing
# and verbose=TRUE.
# Remember to always clear the preprocessing directory before use:
unlink("output/BHPMF_preprocessing", recursive=TRUE)
dir.create("output/BHPMF_preprocessing")
GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
           prediction.level=4, used.num.hierarchy.levels=3,
           mean.gap.filled.output.path="output/BHPMF.mean.gap.filled.txt",
           std.gap.filled.output.path="output/BHPMF.std.gap.filled.txt",
           tmp.dir="output/BHPMF_preprocessing", 
           rmse.plot.test.data=FALSE, verbose=TRUE)

