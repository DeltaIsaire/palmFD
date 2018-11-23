##########################################
# Palm FD project: trait_filling test code
##########################################
#
# This is a companion file listing the R code used to verify the code used for
# trait filling, and to explore/verify inputs and outputs.
# The code in this file is non-essential and provided only for reference and
# bug-fixing purposes.
# BE CAREFUL RUNNING THIS ENTIRE SCRIPT ALL AT ONCE. This is primarily a
# reference file, not necessarily a working script. Also this script
# calls BHPMF::GapFilling multiple times, so total runtime is significant.
# Input files:
#   none.
# Generated output files:
#   output/TEST.BHPMF.mean.txt
#   output/TEST.BHPMF.std.txt
#   output/TEST.BHPMF_filled_3.csv


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
# Initial tests
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
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
# Function disabled to make the script work
#GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
#           prediction.level=4, used.num.hierarchy.levels=3,
#           mean.gap.filled.output.path="output/TEST.BHPMF.mean.txt",
#           std.gap.filled.output.path="output/TEST.BHPMF.std.txt",
#           tmp.dir="output/BHPMF_preprocessing_test", 
#           rmse.plot.test.data=FALSE, verbose=TRUE)


# ---------------------------------
# Gap-fill trait matrix using BHPMF
# Advanced testing
# ---------------------------------
# We begin with the base unimputed trait matrix, which has a single 'oberved'
# trait value for each species, with gaps (NAs).
# The main script generates:
# 1. Matrix gap-filled using genus means.
# 2. Matrix gap-filled using BHPMF with our three traits of interest, extracting
#    from the BHPMF output only the estimates for the gaps in the original matrix.
# Here, we generate:
# 3. Matrix gap-filled using BHPMF with our three traits of interest,
#    using the full output of BHPMF (including estimates of 'observed' trait data)
# 4. Same as #2, but including a fourth dummy trait, with all values = 1
# 5. Same as #2, but including the following additional traits:
#    MaxStemDia_cm, MaxLeafNumber, Max_Rachis_Length_m, Max_Petiole_length_m,
#    AverageFruitWidth_cm

# In all the following, the BHPMF preprocessing directory is
# output/BHPMF_preprocessing_test
# And the output files are named according to
# TEST.BHPMF.mean.number.txt
# TEST.BHPMF.std.number.txt
# Where 'number' is the number of the test (3, 4 or 5).
# And the resulting gap-filled matrix is saved as
# TEST.BHPMF_filled_number.csv

# 3. Gap-filling using the full output of BHPMF.
# This is largely the same code as used in the main script
trait.matrix <- as.matrix(palm.traits[, c("stem.height", "blade.length",
                                          "fruit.length")
                                      ]
                          )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
hierarchy.matrix <- as.matrix(data.frame(species   = sub(pattern=" ",
                                                         replacement="_",
                                                         x=trait.data$SpecName),
                                         genus     = trait.data$accGenus,
                                         tribe     = trait.data$PalmTribe,
                                         subfamily = trait.data$PalmSubfamily)
                              )
to.remove <- which(rowSums(is.na(trait.matrix)) == ncol(trait.matrix))
BHPMF.missing <- palm.traits[to.remove, ]
trait.matrix <- trait.matrix[-to.remove, ]
hierarchy.matrix <- hierarchy.matrix[-to.remove, ]
rm(to.remove)
trait.matrix[, "stem.height"][which(trait.matrix[, "stem.height"] == 0)] <- 0.0001
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
           prediction.level=4, used.num.hierarchy.levels=3,
           mean.gap.filled.output.path="output/TEST_BHPMF.mean.3.txt",
           std.gap.filled.output.path="output/TEST_BHPMF.std.3.txt",
           tmp.dir="output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE, verbose=TRUE)
test.mean.BHPMF <- read.table(file="output/TEST_BHPMF.mean.3.txt",
                              header=TRUE, sep="	")
test.mean.BHPMF <- data.frame(species = as.character(hierarchy.matrix[, 1]),
                              genus   = as.character(hierarchy.matrix[, 2]),
                              test.mean.BHPMF)
write.csv(test.mean.BHPMF, file="output/TEST.BHPMF_filled_3.csv",
          eol="\r\n", row.names=FALSE)


# 4. Same as #2, but including a fourth dummy trait, with all values = 1
# This is largely the same code as used in the main script
trait.matrix <- as.matrix(palm.traits[, c("stem.height", "blade.length",
                                          "fruit.length")
                                      ]
                          )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
hierarchy.matrix <- as.matrix(data.frame(species   = sub(pattern=" ",
                                                         replacement="_",
                                                         x=trait.data$SpecName),
                                         genus     = trait.data$accGenus,
                                         tribe     = trait.data$PalmTribe,
                                         subfamily = trait.data$PalmSubfamily)
                              )
trait.matrix[, "stem.height"][which(trait.matrix[, "stem.height"] == 0)] <- 0.0001
trait.matrix <- as.matrix(data.frame(trait.matrix,
                                     dummy=rep(1, times=length(trait.matrix[, 1]))
                          )          )
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
           prediction.level=4, used.num.hierarchy.levels=3,
           mean.gap.filled.output.path="output/TEST_BHPMF.mean.4.txt",
           std.gap.filled.output.path="output/TEST_BHPMF.std.4.txt",
           tmp.dir="output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE, verbose=TRUE)
test.mean.BHPMF <- read.table(file="output/TEST_BHPMF.mean.4.txt",
                              header=TRUE, sep="	")
test.mean.BHPMF <- data.frame(species = as.character(hierarchy.matrix[, 1]),
                              genus   = as.character(hierarchy.matrix[, 2]),
                              test.mean.BHPMF)
test.filled.BHPMF <- GapFill(palm.traits, test.mean.BHPMF, by="species",
                               fill=c("stem.height", "blade.length",
                                      "fruit.length")
                              )
test.filled.BHPMF <- test.filled.BHPMF[complete.cases(test.filled.BHPMF), ]
write.csv(test.filled.BHPMF, file="output/TEST.BHPMF_filled_4.csv",
          eol="\r\n", row.names=FALSE)


# 5. Same as #2, but including the following additional traits:
#    MaxStemDia_cm, MaxLeafNumber, Max_Rachis_Length_m, Max_Petiole_length_m,
#    AverageFruitWidth_cm
# This is largely the same code as used in the main script
trait.matrix <- as.matrix(cbind(palm.traits[, c("stem.height", "blade.length",
                                                "fruit.length")],
                                trait.data[, c("MaxStemDia_cm", "MaxLeafNumber",
                                               "Max_Rachis_Length_m",
                                               "Max_Petiole_length_m",
                                               "AverageFruitWidth_cm")]
                          )     )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
hierarchy.matrix <- as.matrix(data.frame(species   = sub(pattern=" ",
                                                         replacement="_",
                                                         x=trait.data$SpecName),
                                         genus     = trait.data$accGenus,
                                         tribe     = trait.data$PalmTribe,
                                         subfamily = trait.data$PalmSubfamily)
                              )
to.remove <- which(rowSums(is.na(trait.matrix)) == ncol(trait.matrix))
BHPMF.missing <- palm.traits[to.remove, ]
trait.matrix <- trait.matrix[-to.remove, ]
hierarchy.matrix <- hierarchy.matrix[-to.remove, ]
rm(to.remove)
trait.matrix[, "stem.height"][which(trait.matrix[, "stem.height"] == 0)] <- 0.0001
trait.matrix[, "MaxStemDia_cm"][which(trait.matrix[, "MaxStemDia_cm"] == 0)] <-
  0.0001
trait.matrix[, "Max_Petiole_length_m"][which(trait.matrix[, "Max_Petiole_length_m"]
  == 0)] <- 0.0001
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
           prediction.level=4, used.num.hierarchy.levels=3,
           mean.gap.filled.output.path="output/TEST_BHPMF.mean.5.txt",
           std.gap.filled.output.path="output/TEST_BHPMF.std.5.txt",
           tmp.dir="output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE, verbose=TRUE)
test.mean.BHPMF <- read.table(file="output/TEST_BHPMF.mean.5.txt",
                              header=TRUE, sep="	")
test.mean.BHPMF <- data.frame(species = as.character(hierarchy.matrix[, 1]),
                              genus   = as.character(hierarchy.matrix[, 2]),
                              test.mean.BHPMF)
test.filled.BHPMF <- GapFill(palm.traits, test.mean.BHPMF, by="species",
                               fill=c("stem.height", "blade.length",
                                      "fruit.length")
                              )
test.filled.BHPMF <- test.filled.BHPMF[complete.cases(test.filled.BHPMF), ]
write.csv(test.filled.BHPMF, file="output/TEST.BHPMF_filled_5.csv",
          eol="\r\n", row.names=FALSE)

