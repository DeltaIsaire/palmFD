##################################################
# Palm FD project: Local null model for FD indices
##################################################
#
# In which we generate FD null model three, based on 
# Random communities sampled from the same species pool.
# The definiton of species pool for this null model is the Assemblage
# Dispersion Field (ADF) weighted by proportion of shared species.
# This could be the called the "local" null model.
#
# Input files:
#   dir 'output/observed_FD/' with 100 observed FD files
#   output/observed_FD/community_trait_means_genus_mean.csv
#   output/palm_tdwg3_pres_abs_gapfilled.csv
# Generated output files:
#   < output dir specified below, with 100 FD z-score files >
#   < nullmodel processing dir specified below, containing for each iteration
#     1 file with null FD values >


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)
library(reshape2)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")
source(file = "functions/weighted_ADF_null_model_functions.R")


# Subset FD input for quick code testing?
subset <- FALSE
# Enable verbose reporting in the FD calculation?
verbose <- TRUE
# Null model processing directory (with trailing slash!)
nm.dir <- "output/null_model_local/iterations/"
# Null model output directory (with trailing slash!)
output.dir <- "output/null_model_local/"
# Common first part of filename for all output files:
header <- "FD_z_scores_local_"
# Number of cores to use for parallel processing. Default is 95% of available cores.
num.cores <- 
  if (!is.na(detectCores())) {
    floor(detectCores() * 0.95)
  } else {
    1
    warning("unable to detect cores. Parallel processing is NOT used!")
    cat("unable to detect cores. Parallel processing is NOT used!")
  }


##################
# Data preparation
##################
cat("Preparing data...\n")

# Information on tdwg3 units
# --------------------------
tdwg3.info <- read.csv(file = "output/tdwg3_info.csv", header = TRUE)
# data frame with 368 observations of 9 variables:
# $tdwg3.code   - 3-letter code for each botanical country (tdwg3 unit).
# $tdwg3.name   - Name of the botanical country.
# $realm        - Factor with 3 levels "NewWorld", "OWEast", "OWWest".
# $flora.region - Floristic region to which the botanical country belongs
# $lat          - Latitude.
# $lon          - #Longitude
# is.island     - binary factor indicating whether botanical country is an island
# has.palms     - binary factor indicating whether botanical country has palms
# palm.richness - Palm species richness in the botanical country, including
#                 known absence (richness 0).


# Load transformed gapfilled trait matrices
# -----------------------------------------
traits.mean <-
  read.csv(file = "output/trait_matrices/palm_trait_matrix_genus_mean.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
# This will be useful:
trait.names <- colnames(traits.mean)
traits.gapfilled <- vector("list", length = 100)
for (i in seq_along(traits.gapfilled)) {
  traits.gapfilled[[i]] <-
    read.csv(file = paste0("output/trait_matrices/palm_trait_matrix_filled_",
                           i,
                           ".csv"
                           ),
             row.names = 1,
             check.names = FALSE
             ) %>%
    as.matrix()
}

# Load presence/absence matrix
# ----------------------------
pres.abs.matrix <- 
  read.csv(file = "output/palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()


#######################
# The Local Null Model
#######################
cat("Creating Local null model... (this will take a moment)\n")

if (!dir.exists(output.dir)) {
  cat("creating directory:", output.dir, "\n")
  dir.create(output.dir)
}

# Define species pools using weighted Assemblage Dispersion Fields (ADF)
# ----------------------------------------------------------------------
# That entails a unique species pool for each botanical country.
#
# Determine weighted species pools:
communities <- rownames(pres.abs.matrix)
adf <- ADF(communities, pres.abs.matrix)
species.pools <- WeighADF(adf, pres.abs.matrix)

# Checksum: are all the species pools large enough?
test <- data.frame(richness = rowSums(pres.abs.matrix),
                   pool.size = unlist(llply(species.pools, nrow)),
                   row.names = rownames(pres.abs.matrix)
                   )
if (any(test[, "richness"] > test[, "pool.size"])) {
  stop("Some ADF species pools are too small")
}


# Function to run the local null model
# ------------------------------------
RunLocal <- function(trait.matrix, pres.abs.matrix, id, header) {
#   id: a unique identifier, used for naming the produced files.
#       should match the id used for the observed FD file
#   header: character string giving common first part of the name for all
#           output files

  # Run null model for the given dataset
  local.gapfilled <-
    NullADF(trait.matrix = trait.matrix,
            pres.abs.matrix = pres.abs.matrix,
            species.pools = species.pools,
            process.dir = nm.dir,
            iterations = 1000,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            single.traits = TRUE,
            fast = TRUE
            )

  # Using the null model output, find the z-scores (SES) for the observed FD
  # first, load observed FD indices
  fd.observed <- read.csv(file = paste0("output/observed_FD/",
                                        "FD_observed_",
                                        id,
                                        ".csv"
                                        ),
                          header = TRUE,
                          row.names = 1
                          )
  # Second, calculate the z-scores
  fd.local <- NullTransform(fd.observed, local.gapfilled)
  # Finally, save results
  # remember fd.local is a list
  cat("Saving output...\n")
  for (i in seq_along(fd.local)) {
    write.csv(fd.local[[i]],
              file = paste0(output.dir,
                            header,
                            id,
                            "_",
                            names(fd.local)[i],
                            ".csv"
                            ),
              eol = "\r\n",
              row.names = TRUE
              )
  }

  # When all is completed, delete the process.dir
  # test for existence of the LAST file saved (the sd output)
  if (file.exists(paste0(output.dir,
                         header,
                         id,
                         "_",
                         names(fd.local)[length(fd.local)],
                         ".csv"
                         )
                   )
      ) {
    if (identical(DetectOS(), "windows")) {
      remove.dir <- substr(nm.dir, 1, (nchar(nm.dir) - 1))
    } else {
      remove.dir <- nm.dir
    }
    unlink(remove.dir, recursive = TRUE)
  }

  return (0)
}


# Local null model for genus-mean filled data
# -------------------------------------------
# This is time-consuming: run only if the output does not yet exist
cat("For genus-mean filled dataset:\n")
id <- "genus_mean"
if (!file.exists(paste0(output.dir,
                        header,
                        id,
                        "_",
                        "sds",
                        ".csv"
                        )
                 )
    ) {
  RunLocal(trait.matrix = traits.mean,
           pres.abs.matrix = pres.abs.matrix,
           id = id
           )
}
cat("Done.\n")

# Local null model for stochastic genus-level filled data
# -------------------------------------------------------
# Run how many samples? (max 100)
samples <- 30
# This is time-consuming: run only if the output does not yet exist
cat("For stochastic genus-level filled data:\n")
for (i in seq_len(samples)) {
  cat("Sample", i, "\n")
  id <- i
  if (!file.exists(paste0(output.dir,
                          header,
                          id,
                          "_",
                          "sds",
                          ".csv"
                          )
                   )
      ) {
    RunLocal(trait.matrix = traits.gapfilled[[i]],
             pres.abs.matrix = pres.abs.matrix,
             id = id
             )
  }
}

# Average null model outputs for the stochastic genus-level filled data
# ---------------------------------------------------------------------
cat("Averaging results over", samples, "stochastic runs...\n")
# To get the mean means, mean sds and mean z-scores over the (up to) 100
# stochastic datasets.

for (stat in c("means", "sds", "z.scores")) {
  StochasticMeans(stat = stat,
                  samples = samples,
                  output.dir = output.dir,
                  header = header
                  )
}


cat("Done.\n")

