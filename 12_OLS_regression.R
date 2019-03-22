########################################
# Palm FD project: OLS regression models
########################################
#
# In which the functional diversity indices from the null model outut (response
# variables) are correlated to the environmental data, as well as to species
# richness and proportion of species endemism, using OLS regression.
#
# Because canopy height has some missing data, the analysis is performed with and
# without canopy height; and the spatial distribution of missing values in canopy
# height is investigated.
#
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(corrplot)
library(leaps)
library(car)
library(sf)
library(ggplot2)
library(scales)

source(file = "functions/base_functions.R")
source(file = "functions/OLS_regression_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")

# Model output directory (with trailing slash!)
output.dir <- "output/OLS_models/"


###############
# Load datasets
###############
cat("loading data...\n")

# Functional diversity indices
# ----------------------------
# A list with all the FD index data.
# Note these dataframes include richness and endemism.
fd.indices <- vector("list", length = 0)
fd.names <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
              "FRic.fruit.length", "FDis.all.traits", "FDis.stem.height",
              "FDis.blade.length", "FDis.fruit.length"
              )
for (index in fd.names) {
  fd.indices [[index]] <-
    read.csv(file = paste0("output/FD_summary_", index, ".csv"), row.names = 1)
}


# Environmental data
# ------------------
env <- read.csv(file = "output/tdwg3_environmental_predictors.csv",
                row.names = 1
                )

# Additionally, merge richness and endemism into env
indices <- match(rownames(fd.indices[[1]]), rownames(env))
env[, "palm.richness"] <- fd.indices[[1]] [indices, "palm.richness"]
env[, "endemism"] <- fd.indices[[1]] [indices, "percent.endemism"]
# Richness and endemism may need to be transformed to be more normal
# Using Histogram(), the following transformation was deemed necessary:
env[, "palm.richness"] %<>% log()
# Furthermore, endemism is zero-inflated, which could make it unsuitable for
# regular OLS regression. Keep this in mind.

# Additionally, merge contemporary climate PCA axes into env
clim.pca <- read.csv(file = "output/tdwg3_environmental_predictors_climate_pca_axes.csv",
                     row.names = 1
                     )
indices <- match(rownames(fd.indices[[1]]), rownames(env))
env[, c(paste0("clim.", names(clim.pca)))] <-
  clim.pca[indices, ]

# Be aware that canopy height (CH_Mean and CH_Range) has 14 missing values,
# and the LGM data has one missing value.
# One option is running the analysis with and without canopy height, to see if
# removing these 14 values affects the results.
env.noCH <- env[, -which(names(env) %in% c("CH_Mean", "CH_Range"))]


# TDWG3 data
# ----------
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")

tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")


#########################
# Single-predictor models
#########################
cat("Generating single-predictor models:\n")

if (!dir.exists(output.dir)) {
  cat("creating directory:", output.dir, "\n")
  dir.create(output.dir)
}

# Functions to do the hard work
# -----------------------------
# Wrapper function to run SingleModel() and process and save results
RunMM <- function(responses, predictors, name) {
# Name: character string as unique data identifier. Used as suffix for file names.
  output <- SingleModels(responses, predictors)
  write.csv(output[[1]],
            file = paste0(output.dir, "OLS_single_Rsq_", name, ".csv"),
            eol = "\r\n"
            )
  write.csv(output[[2]],
            file = paste0(output.dir, "OLS_single_slope_", name, ".csv"),
            eol = "\r\n"
            )
  write.csv(output[[3]],
            file = paste0(output.dir, "OLS_single_p_value_", name, ".csv"),
            eol = "\r\n"
            )
  return (0)
}

# Function to extract data of each null model from fd.indices
GetFD <- function(null.model) {
  fd.list <- llply(fd.indices, function(x) { x[, paste0(null.model, ".SES")] } )
  df <- as.data.frame(fd.list)
  rownames(df) <- rownames(fd.indices[[1]])
  df
}


# Call the functions for each null model dataset
# ----------------------------------------------
cat("(1) For null model Global...\n")
fd.global <- GetFD("global")
RunMM(fd.global, env, "global")
RunMM(fd.global, env.noCH, "global_noCH")

cat("(2) For null model Realm...\n")
fd.realm <- GetFD("realm")
RunMM(fd.realm, env, "realm")
RunMM(fd.realm, env.noCH, "realm_noCH")

cat("(3) For null model ADF...\n")
fd.adf <- GetFD("adf")
RunMM(fd.adf, env, "adf")
RunMM(fd.adf, env.noCH, "adf_noCH")


# Summarize results
# -----------------
cat("Summarizing single-predictor model results...\n")
SingleSummary <- function(ch) {
# 1. For each null model, for each of the 8 FD indices, get the R-squared of
# statistically significant predictors, combining results into a list.
#
# ch: logical indicating whether to use the data INcluding canopy height
  canopy <- ifelse(ch, "", "_noCH")
  single.predictors <- vector("list", length = 3)
  names(single.predictors) <- c("global", "realm", "adf")
  for (i in seq_along(single.predictors)) {
    # Load the single-model data
    p.vals <- read.csv(file = paste0(output.dir,
                                     "OLS_single_p_value_",
                                     names(single.predictors)[i],
                                     canopy,
                                     ".csv"
                                     ),
                       row.names = 1
                       )
    rsq.vals <- read.csv(file = paste0(output.dir,
                                       "OLS_single_Rsq_",
                                       names(single.predictors)[i],
                                       canopy,
                                       ".csv"
                                       ),
                         row.names = 1
                         )
    # For each FD index, filter data
    fd.list <- vector("list", length = 8)
    names(fd.list) <- rownames(p.vals)
    for (j in seq_len(nrow(p.vals))) {
      indices <- which(p.vals[j, ] < 0.05)
      best.preds <- rsq.vals[j, indices, drop = FALSE]
      rownames(best.preds) <- "R-squared"
      if (ncol(best.preds) > 0) {
        fd.list[[j]] <- sort(best.preds, decreasing = TRUE)
      } else {
        fd.list[[j]] <- "no significant predictors"
      }
    }
  single.predictors[[i]] <- fd.list
  }
  single.predictors
}

results.single <- SingleSummary(ch = TRUE)
results.single.noCH <- SingleSummary(ch = FALSE)



#############################################
# Investigate missing values in canopy height
#############################################
cat("Investigating missing values in canopy height...\n")
# Same missing values for both CH_Mean and CH_Range.
# Which TDWG3 units are affected?
tdwg.observed <-
  fd.indices[[1]] [!is.na(fd.indices[[1]] [, "observed"]),
                   "observed",
                   drop = FALSE] %>%
  rownames()
tdwg.ch <- env[match(tdwg.observed, rownames(env)), c("CH_Mean", "CH_Range")]
affected <- rownames(tdwg.ch)[!complete.cases(tdwg.ch)]

# What do we know about these places?
tdwg3.info[tdwg3.info$tdwg3.code %in% affected, ]
# Most of these are islands
# Most of these are in Old World East (which has the most islands anyway)
# Palm richness is 4-25, so not huge but not insignificant.

# Visualize the locations
# -----------------------
factor <- rep(TRUE, nrow(tdwg3.info))
indices <- CrossCheck(x = tdwg3.info[, "tdwg3.code"],
                      y = affected,
                      presence = TRUE,
                      value = FALSE
                      )
factor[indices] <- FALSE
factor %<>% as.factor()
factor[!tdwg3.info[, "tdwg3.code"] %in% tdwg.observed] <- NA
ggsave(SpatialPlotFactor(tdwg.map,
                         factor = factor,
                         factor.name =  "Data for CH?",
                         title = "TDWG3 units with data for canopy height"
                         ),
       filename = paste0(output.dir, "canopy_height_missing_values.png"),
       width = 8,
       height = 4
       )
# Not a completely random distribution, but not too clumped either.



#################################
# Multi-predictor model selection
#################################
cat("Selecting best multi-predictor models:\n")

# Predictor selection
# -------------------
# Predictor preprocessing showed that some predictors cannot be used together
# in the same model, due to extreme collinearity.
# Here, these predictors are removed. Endemism is also removed.
env.subset <- env[, !colnames(env) %in% c("srtm_alt_mean", "lgm_ens_Tmean",
                                          "lgm_ens_Pmean", "bio5_mean",
                                          "bio6_mean", "endemism"
                                          )
                  ]
# Now, this needs to be split into two versions: one where climate is described
# by the bioclim variables, and one where PCA axes are used instead.
pca.axes <- paste0("clim.PC", 1:4)
env.subset.clim <- env.subset[, !colnames(env.subset) %in% pca.axes]
bioclim.vars <- paste0("bio", 1:19, "_mean")
env.subset.pca <- env.subset[, !colnames(env.subset) %in% bioclim.vars]


# Function to parse the output of MultiSelect() and select the best model
# ----------------------------------------------------------------------
# On theoretical grounds, the best best model will be the one found by
# cross-validation with the exhaustive search method; because exhaustive search
# does not get stuck in a local optimum in the search space but searches the
# complete space, and cross-validation is the best method to prevent over-fitting.
# The downside is that cross-validation divides the observations randomly
# into k folds, which introduces stochasticity into the model selection process.
#
# The additional 15 selections can be seen as reference models with which the
# single best model can be compared. Taking the means of those 15 models will
# provide a reference.
EvalMM <- function(multimod) {
  # extract best model performance data
  performance <- array(data = NA,
                       dim = c(dim(multimod[[1]] [["performance"]]), 4),
                       dimnames = list(rownames(multimod[[1]] [["performance"]]),
                                       colnames(multimod[[1]] [["performance"]]),
                                       names(multimod)
                                       )
                       )
  for (i in seq_along(multimod)) {
    performance[, , i] <- multimod[[i]] [["performance"]]
  }
  # summarize performance
  mat <- matrix(data = NA,
                nrow = 4,
                ncol = 5,
                dimnames= list(c("best.mod", "mean", "min", "max"),
                               colnames(performance)
                               )
                )
  mat[1, ] <- performance[1, , 1]
  mat[2, ] <- apply(performance, 2, mean)
  mat[3, ] <- apply(performance, 2, min)
  mat[4, ] <- apply(performance, 2, max)
  mat
}


# Apply model selection for all response variables
# ------------------------------------------------
cat("Applying model selection for all response variables:\n")
# Using either (1) bioclim data, or (2) climate pca axes.

# Init output lists
fd.list <- vector("list", length = length(fd.names))
names(fd.list) <- fd.names
null.models <- c("global", "realm", "adf")
formulae.pca <- rep(list(fd.list), 3)
names(formulae) <- null.models
mods.pca <- formulae
performances.pca <- formulae

formulae.clim <- formulae.pca
mods.clim <- mods.pca
performances.clim <- performances.pca

# Run selection process:
cat("(1) for climate as pca axes:\n")
for (model in null.models) {
  for (index in fd.names) {
    cat(model, index, "...\n")
    multimod <-
      MultiSelect(response.var = fd.indices[[index]] [, paste0(model, ".SES")],
                  response.name = "index",
                  predictors = env.subset.pca
                  )
    formulae.pca[[model]] [[index]] <-
      multimod[["exhaustive"]] [["formulae"]] [["cross.validation"]]
    mods.pca[[model]] [[index]] <-
      multimod[["exhaustive"]] [["mods"]] [["cross.validation"]]
    performances.pca[[model]] [[index]] <- EvalMM(multimod)
  }
}
cat("(2) for climate as seperate bioclim predictors:\n")
for (model in null.models) {
  for (index in fd.names) {
    cat(model, index, "...\n")
    multimod <-
      MultiSelect(response.var = fd.indices[[index]] [, paste0(model, ".SES")],
                  response.name = "index",
                  predictors = env.subset
                  )
    formulae.clim[[model]] [[index]] <-
      multimod[["exhaustive"]] [["formulae"]] [["cross.validation"]]
    mods.clim[[model]] [[index]] <-
      multimod[["exhaustive"]] [["mods"]] [["cross.validation"]]
    performances.clim[[model]] [[index]] <- EvalMM(multimod)
  }
}


cat("Done.\n")

