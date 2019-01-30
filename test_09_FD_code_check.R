##########################################
# Palm FD project: Verification of FD code
##########################################
#
# In which the code for calculating FD indices and null models (mostly the custom
# functions in 'functions/functional_diversity_functions.R') is verified
# through the use of artificial data.
#
# In addition, this script assesses how often random species samples for each
# null model have to be redrawn because the S > T criterion is not satisfied.
#
# Input files:
#   output/test/tdwg3_info.csv
#   output/test/test_palm_tdwg3_pres_abs_gapfilled.csv
#   output/test/test_palm_tdwg3_pres_abs_unfilled.csv
#   output/test/test_palm_traits_transformed_gapfilled.csv
#   output/test/test_palm_traits_transformed_unfilled.csv
# Generated output files:
#   output/test/test_nullmodel_species_resampling_freq_gapfilled.csv
#   output/test/test_nullmodel_species_resampling_freq_unfilled.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)
library(reshape2)
library(sf)
library(spdep)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")


# Enable verbose reporting in the null model procedure?
verbose <- TRUE
# Number of cores to use for parallel processing (via parallel::mclapply).
# Default is 1 less than the number of available cores (so your computer
# doesn't lock up completely while running the code), and should work on most
# machines.
num.cores <- if (is.na(detectCores())) { 1 } else { max(1, detectCores() - 1) }
# Subset FD input for quick code testing?
subset <- FALSE


###################################
# Data preparation: artificial data
###################################
cat("Preparing data...\n")

# Test trait matrix
# -----------------
# Three traits, each with 4 unique values.
# Two species for every possible trait combination, so 4^3 * 2 = 128 species.
# This way not all species are unique.
higt <- c(2, 4, 6, 8)  # mean 5, sd 2.6, range 6.
blaed <- c(7, 9, 11, 13)  # mean 10, sd 2.58, range 6
fruti <- c(4, 8, 12, 16)  # mean 10, sd 5.16, range 12
# Values chosen so we can compare differing mean and differing sd.
unique <- data.frame(a = sort(rep(higt, 4 ^ 2)),
                     b = rep(sort(rep(blaed, 4)), 4),
                     c = rep(fruti, 4 ^2 )
                     )
test.traits <-
  matrix(nrow = 128,
         ncol = 3,
         byrow = FALSE,
         dimnames = list(row = paste0("species.", 1:128),
                         col = c("higt", "blaed", "fruti")
                         ),
         data = c(rep(unique[, "a"], 2),
                  rep(unique[, "b"], 2),
                  rep(unique[, "c"], 2)
                  )
         )
test.traits

# Test presence/absence matrix
# ----------------------------
# A small set of carefully chosen communities:
#   area.A  - 16 species, all with the same higt and every possible combination
#             of blaed and fruti.
#   area.B  - Like A but subsetted to fruti == c(4, 8) (only low fruti)
#   area.C  - Like A but subsetted to fruti == c(12, 16) (only high fruti)
#   area.D  - like A but subsetted to fruti == c(4, 16) (higher fruti range)
#   Area.A1 through D1: Like A-D but with two species for each trait combination.
indices <- data.frame(A  = 1:16,
                      B  = sort(c(1 + 4 * 0:3, 2 + 4 * 0:3)),
                      C  = sort(c(3 + 4 * 0:3, 4 + 4 * 0:3)),
                      D  = sort(c(1 + 4 * 0:3, 4 + 4 * 0:3)),
                      A1 = c(1:16, 1:16 + 64),
                      B1 = c(sort(c(1 + 4 * 0:3, 2 + 4 * 0:3)),
                             sort(c(1 + 4 * 0:3, 2 + 4 * 0:3)) + 64
                             ),
                      C1 = c(sort(c(3 + 4 * 0:3, 4 + 4 * 0:3)),
                             sort(c(3 + 4 * 0:3, 4 + 4 * 0:3)) + 64
                             ),
                      D1 = c(sort(c(1 + 4 * 0:3, 4 + 4 * 0:3)),
                             sort(c(1 + 4 * 0:3, 4 + 4 * 0:3)) + 64
                             )
                      )
test.pres.abs <-
  matrix(nrow = 8,
         ncol = 128,
         dimnames = list(row = colnames(indices),
                         col = rownames(test.traits)
                         )
         )
for (i in seq_along(indices)) {
  test.pres.abs[i, indices[, i]] <- 1
}
test.pres.abs[is.na(test.pres.abs)] <- 0


#########################
# Test: base FD functions
#########################
# 'RunFD' and 'SingleFD'
missing <- which(colSums(test.pres.abs) == 0)
test.pres.abs %<>% .[, -missing]
test.traits %<>% .[-missing, ]

RunFD.output <- RunFD(test.traits,
                      test.pres.abs,
                      subset = FALSE,
                      verbose = TRUE
                      )
# A-D and A1-D1 have the same FRic and FDis, so duplicate species do not affect
# FD indices. 

SingleFD.output <- SingleFD(test.traits[, -1],
                            test.pres.abs,
                            subset = FALSE,
                            verbose = TRUE
                            )
# FD:dbFD automatically removes traits with no variation. SingleFD only presents
# one trait at a time to FD:dbFD, so then if that trait has no variation,
# an error occurs. The lesson here is: be sure to only include traits with
# variation!
# Here that means excluding higt
#
# results:
# $blaed:
#   FRic is the same for all areas, as it should be. Each area has the same 4
#   unique trait values.
#   FDis is also the same for all areas, for the same reason.
# $fruti:
#   FRic reflects the trait range, being high for A, D, A1 and D1 and low for
#   the other areas. This is expected.
#   FDis reflects differing sd of the fruti distributions, again as expected.


###############
# RandomSpecies
###############
richness <- data.frame(community = rownames(test.pres.abs),
                       richness  = rowSums(test.pres.abs)
                       )

randspec.1 <- RandomSpecies(communities = list(A = c("A")),
                            richness = richness,
                            species = list(A = colnames(test.pres.abs)[1:16]),
                            trait.matrix = test.traits[, -1],
                            verbose = TRUE,
                            groups = FALSE
                            )
# I have played around with this to test it.
# IF the richness of an area is higher than the number of unique species in its
# species pool, NAs are generated, but species pools should always be large enough.
#
# Above, we learned that all included traits should have variation. Below, we
# learn that this same problem shows up when running Null Models. That means
# the variance check in RandomSpecies should not only test for <4 unique trait
# combinations, but also for trait invariance in individual traits.
# I have implemented this.
#
# Consequently, RandomSpecies only runs to completion if we remove the trait 'higt'


###########
# NullModel
###########
unlink("output/test/nullmodel_test/", recursive = TRUE)
NullModel.output <- NullModel(trait.matrix = test.traits[, -1],
                              pres.abs.matrix = test.pres.abs,
                              groups = list(group = rownames(test.pres.abs)),
                              process.dir = "output/test/nullmodel_test/",
                              iterations = 10,
                              mc.cores = num.cores,
                              subset = FALSE,
                              verbose = TRUE,
                              random.groups = TRUE,
                              single.traits = TRUE
                              )
# Getting this function to run made me change some of the argument defaults
# of this function, add one line of code, and generalize the result parsing
# code. Very useful.
# Output:
# $FRic$all.traits - very minor variation between trials, which makes sense with
#                    this data. A1-D1 all have higher FRic samples than A-D,
#                    likely because they have higher species richness. Not a higher
#                    number of unique species though. A random species sample
#                    is not likely to have duplicates for every species, which
#                    explains the higher FRic. THIS MEANS COMMUNITIES WITH HIGH
#                    TRAIT REDUNDANCY HAVE LOW FRIC COMPARED TO THEIR NULLMODEL,
#                    which actually makes sense. So no issues.
# FRic$stem.height - trait name is hardcoded, but we can live with that. All values
#                    here are the same

# Next trial: can 'NullModel' be run for a single area?
# Of course it can, although I had to tweak the code slightly to make it work.
# Here is how:
unlink("output/test/nullmodel_test/", recursive = TRUE)
NullModel.output <- NullModel(trait.matrix = test.traits[, -1],
                              pres.abs.matrix = test.pres.abs,
                              groups = list(A = c("A", "A1")),
                              process.dir = "output/test/nullmodel_test/",
                              iterations = 10,
                              mc.cores = num.cores,
                              subset = FALSE,
                              verbose = TRUE,
                              random.groups = FALSE,
                              single.traits = TRUE
                              )
# use 'groups' argument list(area = <areas>), with 'random.groups' = FALSE




###################################
# Species sampling redraw frequency
###################################
cat("Assessing frequency of random species sample redraws:\n")
cat("Preparing data...\n")
# Number of samples to test with:
samples <- 1000

# Information on tdwg3 units
tdwg3.info <- read.csv(file = "output/test/tdwg3_info.csv")

# Load datasets
pres.abs.gapfilled <- 
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
pres.abs.unfilled <- 
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_unfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
traits.gapfilled <- 
  read.csv(file = "output/test/test_palm_traits_transformed_gapfilled.csv",
           row.names = 1
           ) %>%
  as.matrix()
traits.unfilled <- 
  read.csv(file = "output/test/test_palm_traits_transformed_unfilled.csv",
           row.names = 1
           ) %>%
  as.matrix()


# Global null model
# -----------------
cat("1. For global null model\n")
global.tdwg3 <- list(global = tdwg3.info[, "tdwg3.code"])

resample.global.gapfilled <-
  ResampleTest(pres.abs.matrix = pres.abs.gapfilled,
               trait.matrix = traits.gapfilled,
               groups = global.tdwg3,
               random.groups = TRUE,
               verbose = TRUE,
               samples = samples
               )

resample.global.unfilled <-
  ResampleTest(pres.abs.matrix = pres.abs.unfilled,
               trait.matrix = traits.unfilled,
               groups = global.tdwg3,
               random.groups = TRUE,
               verbose = TRUE,
               samples = samples
               )

# Regional null model
# -------------------
cat("For Regional null model\n")
realm.tdwg3 <- list(new.world = tdwg3.info[tdwg3.info$realm == "NewWorld",
                                           "tdwg3.code"],
                    old.world.west = tdwg3.info[tdwg3.info$realm == "OWWest",
                                                "tdwg3.code"],
                    old.world.east = tdwg3.info[tdwg3.info$realm == "OWEast",
                                                "tdwg3.code"]
                    )

resample.realm.gapfilled <-
  ResampleTest(pres.abs.matrix = pres.abs.gapfilled,
               trait.matrix = traits.gapfilled,
               groups = realm.tdwg3,
               random.groups = TRUE,
               verbose = TRUE,
               samples = samples
               )

resample.realm.unfilled <-
  ResampleTest(pres.abs.matrix = pres.abs.unfilled,
               trait.matrix = traits.unfilled,
               groups = realm.tdwg3,
               random.groups = TRUE,
               verbose = TRUE,
               samples = samples
               )

# Local null model
# ----------------
# TODO: to be added


# Gather and save results
# -----------------------
resample.freq.gapfilled <-
  data.frame(community = resample.global.gapfilled[, "community"],
             global.count = resample.global.gapfilled[, 2],
             global.freq = resample.global.gapfilled[, 3],
             realm.count = resample.realm.gapfilled[, 2],
             realm.freq = resample.realm.gapfilled[, 3]
             )
write.csv(resample.freq.gapfilled,
          file = "output/test/test_nullmodel_species_resampling_freq_gapfilled.csv",
          eol = "\r\n",
          row.names = FALSE
          )

resample.freq.unfilled <-
  data.frame(community = resample.global.unfilled[, "community"],
             global.count = resample.global.unfilled[, 2],
             global.freq = resample.global.unfilled[, 3],
             realm.count = resample.realm.unfilled[, 2],
             realm.freq = resample.realm.unfilled[, 3]
             )
write.csv(resample.freq.unfilled,
          file = "output/test/test_nullmodel_species_resampling_freq_unfilled.csv",
          eol = "\r\n",
          row.names = FALSE
          )

cat("Done.\n")

