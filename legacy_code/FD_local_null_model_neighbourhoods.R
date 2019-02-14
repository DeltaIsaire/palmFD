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


# Define neighbourhood for each TDWG3 unit
# ----------------------------------------
# Load spatial map:
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")
tdwg.spatial <- as(tdwg.map, "Spatial")
# Find neighbours:
local <- poly2nb(tdwg.spatial, queen = TRUE)
# we don't want a fancy neighbour object, just a list
class(local) <- "list"
# This is a list of TDWG3 ID numbers. We can find the corresponding TDWG3 codes:
tdwg3.code <- as.data.frame(tdwg.map)[, "LEVEL_3_CO"]
# Then do substitution magic
local.tdwg3 <- llply(local, function(x) { tdwg3.code[x] } )
names(local.tdwg3) <- tdwg3.code
# Species pool includes not only neighbours, but also the area itself,
# so add these to the groups
for (i in seq_along(local.tdwg3)) {
  local.tdwg3[[i]] <- c(local.tdwg3[[i]], names(local.tdwg3)[i])
  local.tdwg3[[i]] %<>% .[order(.)]
}
# Not done yet: species pool should include at least one neighbour, not solely
# the tdwg3 unit itself. This happens for islands, which have no direct
# borders and therefore no neighbours.
# Solution: for each tdwg3 unit with length (species pool) == 1, include
# the nearest neighbour.
near.neigh <- 
  coordinates(tdwg.spatial) %>%
  knearneigh(., k = 1) %>%
  knn2nb()
class(near.neigh) <- "list"
# And merge:
for (i in seq_along(local.tdwg3)) {
  if (identical(length(local.tdwg3[[i]]), as.integer(1))) {
    local.tdwg3[[i]] %<>%
      c(., tdwg3.code[near.neigh[[i]]]) %>%
      unique()
  }
}

local.tdwg3
# List of 368:
#   For each of the 368 TDWG3 units, a character vector giving the names of the
#   TDWG3 unit plus the surrounding TDWG3 units OR the nearest neighbouring
#   TDWG3 unit.

# Run the null model simulations
# ------------------------------
local.gapfilled <-
  NullModel(trait.matrix = traits.gapfilled,
            pres.abs.matrix = pres.abs.gapfilled,
            groups = local.tdwg3,
            process.dir = "output/test/nullmodel_local/gapfilled/",
            iterations = 100,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = FALSE
            )

local.unfilled <-
  NullModel(trait.matrix = traits.unfilled,
            pres.abs.matrix = pres.abs.unfilled,
            groups = local.tdwg3,
            process.dir = "output/test/nullmodel_local/unfilled/",
            iterations = 100,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = FALSE
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
fd.local.gapfilled <- NullTransform(fd.gapfilled, local.gapfilled)
  write.csv(fd.local.gapfilled,
            file = "output/test/test_fd_z_scores_local_gapfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )
fd.local.unfilled <- NullTransform(fd.unfilled, local.unfilled)
  write.csv(fd.local.unfilled,
            file = "output/test/test_fd_z_scores_local_unfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )

cat("Done.\n")
