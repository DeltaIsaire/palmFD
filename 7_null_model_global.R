###################################################
# Palm FD project: Functional Diversity calculation
###################################################
#
# In which we generate FD null model two, based on 
# Random communities sampled from the same species pool.
# The definiton of species pool for this null model is all palm species
# occurring worldwide.
# Which could be the called the "global" null model.
#
# Input files:
#   output/tdwg3_info.csv
#   output/FD_palm_tdwg3_pres_abs_gapfilled.csv
#   output/FD_palm_tdwg3_pres_abs_unfilled.csv
#   output/FD_palm_traits_transformed_gapfilled.csv
#   output/FD_palm_traits_transformed_unfilled.csv
#   output/FD_indices_gapfilled.csv
#   output/FD_indices_unfilled.csv
# Generated output files:
#   Null model processing directory: output/nullmodel_global/
#   output/FD_z_scores_global_gapfilled.csv
#   output/FD_z_scores_global_unfilled.csv
#   graphs/nullmodel_global_z_convergence_FRic.svg
#   graphs/nullmodel_global_z_convergence_FDis.svg


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)
library(reshape2)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")
source(file = "functions/plotting_functions.R")


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
tdwg3.info <- read.csv(file = "output/tdwg3_info.csv")

# Load datasets
# -------------
pres.abs.gapfilled <- 
  read.csv(file = "output/FD_palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
pres.abs.unfilled <- 
  read.csv(file = "output/FD_palm_tdwg3_pres_abs_unfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
traits.gapfilled <- 
  read.csv(file = "output/FD_palm_traits_transformed_gapfilled.csv",
           row.names = 1
           ) %>%
  as.matrix()
traits.unfilled <- 
  read.csv(file = "output/FD_palm_traits_transformed_unfilled.csv",
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
            process.dir = "output/nullmodel_global/gapfilled/",
            iterations = 400,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = TRUE
            )

global.unfilled <-
  NullModel(trait.matrix = traits.unfilled,
            pres.abs.matrix = pres.abs.unfilled,
            groups = global.tdwg3,
            process.dir = "output/nullmodel_global/unfilled/",
            iterations = 400,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = TRUE
            )


############################################################
# Use null model output to convert raw FD values to z-scores
############################################################
cat("Applying regional null model to raw FD indices...\n")

# Load raw FD indices
# -------------------
fd.gapfilled <- read.csv(file = "output/FD_indices_gapfilled.csv",
                      header = TRUE,
                      row.names = 1
                      )
fd.unfilled <- read.csv(file = "output/FD_indices_unfilled.csv",
                        header = TRUE,
                        row.names = 1
                        )

# Convert FD to z-cores and save output
# -------------------------------------
fd.global.gapfilled <- NullTransform(fd.gapfilled, global.gapfilled)
  write.csv(fd.global.gapfilled,
            file = "output/FD_z_scores_global_gapfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )
fd.global.unfilled <- NullTransform(fd.unfilled, global.unfilled)
  write.csv(fd.global.unfilled,
            file = "output/FD_z_scores_global_unfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )

###############################
# Visualize z-score convergence
###############################
cat("Visualizing z-score convergence...\n")
# For a few carefully chosen communities, we can visualize the convergence
# of their z-scores as the number of iterations increases.
# How? By computing the z-score incrementally.
# Which communities? The ones with highest/lowest raw FRic and FDis.
# That's MLY, WAU, SCZ and TCI, respectively.
areas <- sort(c("MLY", "WAU", "SCZ", "TCI"))
#
# For practical purposes, we can focus only on the all-traits FD for the gapfilled
# dataset. The unfilled and single-trait data is assumed to converge in a similar
# fashion.

# Subset data to a form usable by ZScore()
# FRic:
temp <- fd.gapfilled[rownames(fd.gapfilled) %in% areas, 1]
raw.FRic.subset <- structure(temp, names = areas)
null.FRic.subset <-
  global.gapfilled[[1]] [[1]] [, colnames(global.gapfilled[[1]] [[1]]) %in%
                                   areas]
# FDis:
temp <- fd.gapfilled[rownames(fd.gapfilled) %in% areas, 5]
raw.FDis.subset <- structure(temp, names = areas)
null.FDis.subset <-
  global.gapfilled[[2]] [[1]] [, colnames(global.gapfilled[[2]] [[1]]) %in%
                                   areas]

# Initialize output
# z.scores require at least 2 iterations, so the first one doesn't count
z.scores.FRic <- matrix(ncol = length(areas),
                        nrow = nrow(null.FRic.subset) - 1,
                        dimnames = list(row = seq_len(nrow(null.FRic.subset) - 1),
                                        col = areas
                                        )
                        )
z.scores.FDis <- matrix(ncol = length(areas),
                        nrow = nrow(null.FDis.subset) - 1,
                        dimnames = list(row = seq_len(nrow(null.FDis.subset) - 1),
                                        col = areas
                                        )
                        )

# calculate incremental means
# Instead of NullTransform() we use ZScore() directly.
# z.scores require at least 2 iterations, so the first one doesn't count
for (i in seq_len(nrow(null.FRic.subset))[-1]) {
  z.scores.FRic[(i - 1), ] <-
    ZScore(raw.FRic.subset, null.FRic.subset[(1:i), ,drop = FALSE])
}
for (i in seq_len(nrow(null.FDis.subset))[-1]) {
  z.scores.FDis[(i - 1), ] <-
    ZScore(raw.FDis.subset, null.FDis.subset[(1:i), ,drop = FALSE])
}

# Visualize results
MyPlot <- function(z.scores, index) {
  par(mfrow = c(2, 2))
  for (i in seq_along(areas)) {
    Scatterplot(x = seq_len(nrow(z.scores)),
                y = z.scores[, i],
                xlab = "iterations",
                ylab = "z.score",
                title = paste(index,
                              "z.score convergence in",
                              colnames(z.scores)[i]
                              )
                )
  }
}
GraphSVG(MyPlot(z.scores.FRic, "FRic"),
         file = "graphs/nullmodel_global_z_convergence_FRic.svg",
         width = 8,
         height = 8
         )
GraphSVG(MyPlot(z.scores.FDis, "FDis"),
         file = "graphs/nullmodel_global_z_convergence_FDis.svg",
         width = 8,
         height = 8
         )

cat("Done.\n")

