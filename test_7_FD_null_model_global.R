#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we generate FD null model two, based on 
# Random communities sampled from the same species pool.
# The definiton of species pool for this null model is all palm species
# occurring worldwide.
# Which could be the called the "global" null model.
#
# Input files:
#   output/test/tdwg3_info.csv
#   output/test/test_palm_tdwg3_pres_abs_gapfilled.csv
#   output/test/test_palm_tdwg3_pres_abs_unfilled.csv
#   output/test/test_palm_traits_transformed_gapfilled.csv
#   output/test/test_palm_traits_transformed_unfilled.csv
#   output/test/test_fd_indices_gapfilled.csv
#   output/test/test_fd_indices_unfilled.csv
# Generated output files:
#   Null model processing directory: output/test/nullmodel_global/
#   output/test/test_fd_z_scores_global_gapfilled.csv
#   output/test/test_fd_z_scores_global_unfilled.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)
library(reshape2)

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


##################
# Data preparation
##################
cat("Preparing data...\n")

# Information on tdwg3 units
# --------------------------
tdwg3.info <- read.csv(file = "output/test/tdwg3_info.csv")

# Load datasets
# -------------
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
# This will be useful:
trait.names <- colnames(traits.gapfilled)


#######################
# The Global Null Model
#######################
cat("Creating Global null model... (this will take a while)\n")

# Define grouping of tdwg3 units into realms
# ------------------------------------------
global.tdwg3 <- list(global = tdwg3.info[, "tdwg3.code"])

# Run the null model simulations
# ------------------------------
global.gapfilled <-
  NullModel(trait.matrix = traits.gapfilled,
            pres.abs.matrix = pres.abs.gapfilled,
            groups = global.tdwg3,
            process.dir = "output/test/nullmodel_global/gapfilled/",
            iterations = 100,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose
            )

global.unfilled <-
  NullModel(trait.matrix = traits.unfilled,
            pres.abs.matrix = pres.abs.unfilled,
            groups = global.tdwg3,
            process.dir = "output/test/nullmodel_global/unfilled/",
            iterations = 100,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose
            )


############################################################
# Use null model output to convert raw FD values to z-scores
############################################################
cat("Applying regional null model to raw FD indices...\n")

# Load raw FD indices
# -------------------
fd.gapfilled <- read.csv(file = "output/test/test_fd_indices_gapfilled.csv",
                      header = TRUE,
                      row.names = 1
                      )
fd.unfilled <- read.csv(file = "output/test/test_fd_indices_unfilled.csv",
                        header = TRUE,
                        row.names = 1
                        )

# Convert FD to z-cores and save output
# -------------------------------------
fd.global.gapfilled <- NullTransform(fd.gapfilled, global.gapfilled)
  write.csv(fd.global.gapfilled,
            file = "output/test/test_fd_z_scores_global_gapfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )
fd.global.unfilled <- NullTransform(fd.unfilled, global.unfilled)
  write.csv(fd.global.unfilled,
            file = "output/test/test_fd_z_scores_global_unfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )

cat("Done.\n")

