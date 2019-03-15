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
source(file = "functions/plotting_functions_ggplot2.R")

# Model output directory (with trailing slash!)
output.dir <- "output/OLS_models/"


###############
# Load datasets
###############
cat("loading data...\n")

# Functional diversity indices
# ----------------------------
fd.global <- 
  read.csv(file = "output/null_model_global/FD_z_scores_global_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  .[complete.cases(.), ]

fd.realm <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  .[complete.cases(.), ]

fd.adf <-
  read.csv(file = "output/null_model_local/FD_z_scores_local_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  .[complete.cases(.), ]



# Environmental data
# ------------------
env <- read.csv(file = "output/tdwg3_environmental_predictors.csv",
                row.names = 1
                )
# Additionally, merge richness and endemism into the data
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")
env$tdwg3.code <- rownames(env)
env <- merge(x = env,
             y = tdwg3.info[, c("tdwg3.code", "experimental.richness", "endemism")],
             by = "tdwg3.code",
             all = TRUE
             )
rownames(env) <- env[, "tdwg3.code"]
env <- env[, -1]
# Be aware that canopy height (CH_Mean and CH_Range) have 14 missing values,
# and the LGM variables have one missing value.


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
# Function to generate single-predictor model data
MultiModel <- function(fd, predictors) {
  # init output list
  mat <- matrix(data = NA,
                nrow = ncol(fd),
                ncol = ncol(predictors),
                dimnames = list(colnames(fd), colnames(predictors))
                )
  output <- list(mat, mat, mat)
  names(output) <- c("Rsq", "slope", "p-value")

  # Generate single-predictor models for each FD metric
  for (i in seq_along(fd)) {
    response <- colnames(fd)[i]
    model.data <- predictors
    model.data[, response] <- fd[, response]
    mat <- SingleOLS(model.data, response, standardize = TRUE)
    output[[1]] [i, ] <- mat[, 1]
    output[[2]] [i, ] <- mat[, 2]
    output[[3]] [i, ] <- mat[, 3]
  }
  output
}

# Wrapper function to run MultiModel and process and save results
RunSingle <- function(fd, predictors, name) {
# Name: character string as unique data identifier. Used as suffix for file names.
  single <- MultiModel(fd, predictors)
  write.csv(single[[1]],
            file = paste0(output.dir, "OLS_single_Rsq_", name, ".csv"),
            eol = "\r\n"
            )
  write.csv(single[[2]],
            file = paste0(output.dir, "OLS_single_slope_", name, ".csv"),
            eol = "\r\n"
            )
  write.csv(single[[3]],
            file = paste0(output.dir, "OLS_single_p_value_", name, ".csv"),
            eol = "\r\n"
            )
  return (0)
}


# Global null model
# -----------------
cat("(1) For null model Global...\n")
# Subset environmental data to the extent of the FD data
indices <- CrossCheck(x = rownames(env),
                      y = rownames(fd.global),
                      presence = TRUE,
                      value = FALSE
                      )
env.global <- env[indices, ]
env.global.noCH <- env[indices, -which(names(env) %in% c("CH_Mean", "CH_Range"))]

RunSingle(fd.global, env.global, "global")
RunSingle(fd.global, env.global.noCH, "global_noCH")


# Realm null model
# ----------------
cat("(2) For null model Realm...\n")
# Subset environmental data to the extent of the FD data
indices <- CrossCheck(x = rownames(env),
                      y = rownames(fd.realm),
                      presence = TRUE,
                      value = FALSE
                      )
env.realm <- env[indices, ]
env.realm.noCH <- env[indices, -which(names(env) %in% c("CH_Mean", "CH_Range"))]

RunSingle(fd.realm, env.realm, "realm")
RunSingle(fd.realm, env.realm.noCH, "realm_noCH")


# adf null model
# --------------
cat("(3) For null model ADF...\n")
# Subset environmental data to the extent of the FD data
indices <- CrossCheck(x = rownames(env),
                      y = rownames(fd.adf),
                      presence = TRUE,
                      value = FALSE
                      )
env.adf <- env[indices, ]
env.adf.noCH <- env[indices, -which(names(env) %in% c("CH_Mean", "CH_Range"))]

RunSingle(fd.adf, env.adf, "adf")
RunSingle(fd.adf, env.adf.noCH, "adf_noCH")


################################################
# Filtering predictors by collinearity using VIF
################################################
cat("Filtering predictors by collinearity using VIF:\n")

# Function to do the hard work
# ----------------------------
FilterPred <- function(fd, predictors, name) {
# Name: character string as unique data identifier. Used as suffix for file names.

  # init output list
  output <- vector("list", length = ncol(fd))
  names(output) <- colnames(fd)
  # Select non-collinear predictors for each FD metric
  for (i in seq_along(fd)) {
    response <- colnames(fd)[i]
    model.data <- predictors
    model.data[, response] <- fd[, response]
    output[[i]] <- AutoVIF(model.data, response, standardize = TRUE, threshold = 3)
  }
  # parse and save output
  if (length(unique(output)) == 1) {
    cat("Selected predictors are identical for each response variable.\n")
    write.csv(output[[1]],
              file = paste0(output.dir,
                            "OLS_noncollinear_predictors_", 
                            name,
                            ".csv"
                            ),
              eol = "\r\n"
              )
  } else {
    cat("Selected predictors differ between response variables.",
        "Output is NOT saved. We apologize for the inconvenience :(\n"
        )
  }
  output
}

# Apply the function
# ------------------
# The input data is already generated above, but we have to exclude predictors
# with linear dependencies. Hence,
# we use bio7 instead of bio5 + bio6,
# and LGM climate anomalies instead of absolute LGM climate.
exclude <- c("bio5_mean", "bio6_mean", "lgm_ens_Tmean", "lgm_ens_Pmean")

cat("(1) for null model Global...\n")
env.global %<>% .[, -which(names(.) %in% exclude)]
env.global.noCH %<>% .[, -which(names(.) %in% exclude)]

a <- FilterPred(fd.global, env.global, "global")
FilterPred(fd.global, env.global.noCH, "global_noCH")

cat("(2) for null model Realm...\n")
env.realm %<>% .[, -which(names(.) %in% exclude)]
env.realm.noCH %<>% .[, -which(names(.) %in% exclude)]

a <- FilterPred(fd.realm, env.realm, "realm")
FilterPred(fd.realm, env.realm.noCH, "realm_noCH")

cat("(2) for null model ADF...\n")
env.adf %<>% .[, -which(names(.) %in% exclude)]
env.adf.noCH %<>% .[, -which(names(.) %in% exclude)]

a <- FilterPred(fd.adf, env.adf, "adf")
FilterPred(fd.adf, env.adf.noCH, "adf_noCH")


#############################################
# Investigate missing values in canopy height
#############################################
cat("Investigating missing values in canopy height...\n")
# Same missing values for both CH_Mean and CH_Range.
# Which TDWG3 units are affected:
affected <- CrossCheck(x = rownames(fd.global),
                       y = rownames(env)[is.na(env[, "CH_Mean"])],
                       presence = TRUE,
                       value = TRUE
                       )

# What do we know about these places:
tdwg3.info[tdwg3.info$tdwg3.code %in% affected, ]
# Most of these are islands
# Most of these are in Old World East
# Palm richness is 4-25, so not huge but not insignificant.

# Visualize the locations
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")
factor <- rep(TRUE, nrow(tdwg3.info))
indices <- CrossCheck(x = tdwg3.info[, "tdwg3.code"],
                      y = affected,
                      presence = TRUE,
                      value = FALSE
                      )
factor[indices] <- FALSE
factor %<>% as.factor()
factor[!tdwg3.info[, "tdwg3.code"] %in% rownames(fd.global)] <- NA
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




cat("Done.\n")

