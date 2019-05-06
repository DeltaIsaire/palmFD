########################################
# Palm FD project: OLS regression models
########################################
#
# In which the functional diversity indices from the null model output (response
# variables) are correlated to the environmental data, as well as to species
# richness and proportion of species endemism, using OLS regression.
#
# Canopy height has 14 missing values. Here, we always include canopy height.
# In other words, the full dataset is subsetted to complete cases.
#
# This new script focuses on a selection of predictor variables relevant to
# our specific hypothesis (i.e. environmental heterogeneity variables). Furthermore,
# For functional dispersion (FDis) we consider only the observed FDis, not the
# nullmodel-corrected values (SES).
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

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/OLS_regression_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")

# Model output directory (with trailing slash!)
output.dir <- "output/OLS_models_new/"


if (!dir.exists(output.dir)) {
  cat("creating directory:", output.dir, "\n")
  dir.create(output.dir)
}

set.seed(999)


###########################
# Load and prepare datasets
###########################
cat("Preparing data...\n")

# TDWG3 data
# ----------
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")

tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")

# List with the tdwg3 units in each realm
realm.tdwg3 <- vector("list", length = 3)
names(realm.tdwg3) <- levels(tdwg3.info[, "realm"])
for (i in 1:3) {
  realm.tdwg3[[i]] <-
    tdwg3.info[tdwg3.info[, "realm"] == names(realm.tdwg3)[i], "tdwg3.code"]
}


# Functional diversity indices
# ----------------------------
# A list with all the FD index data.
# Note these dataframes include richness and endemism variables.
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
env <- read.csv(file = "output/tdwg3_environmental_predictors_new.csv",
                row.names = 1
                )


# Merge meta-predictors realm, richness and endemism into env
# -----------------------------------------------------------
# Richness and endemism
indices <- match(rownames(fd.indices[[1]]), rownames(env))
env[, "palm.richness"] <- fd.indices[[1]] [indices, "palm.richness"]
env[, "endemism"] <- fd.indices[[1]] [indices, "percent.endemism"]
# Richness and endemism may need to be transformed to be more normal
# Using Histogram(), the following transformation was deemed necessary:
env[, "palm.richness"] %<>% log10()
# Furthermore, endemism is zero-inflated, although the zeros are true zeros.

# realm
indices <- match(rownames(env), tdwg3.info[, "tdwg3.code"])
env[, "realm"] <- tdwg3.info[indices, "realm"]

# Account for NAs in predictors AND responses.
# Subsetting to complete cases would reduce the number of rows (tdwg3.units), which
# hinders data merging. Instead, propagate NAs across rows.
# Check for NAs in predictors:
env.complete <- env
env.complete[!complete.cases(env.complete), ] <- NA
# 14 NAs from canopy height + 1 NA from LGM anomalies.
# We also want a version without canopy height:
env.complete.noch <- env[, !colnames(env) %in% "CH_SD"]
env.complete.noch[!complete.cases(env.complete.noch), ] <- NA
# Only one NA, from LGM anomalies
#
# Additionally, check for NAs in responses: the exclusion of 100% endemic
# communities for the ADF null model
indices <- !complete.cases(fd.indices[["FRic.all.traits"]] [, "adf.SES"])
env.complete[indices, ] <- NA
# 4 NAs from 100% endemism, but 3 of those overlap with canopy height NAs.
# Final sample size is 116  # sum(complete.cases(env.complete))
env.complete.noch[indices, ] <- NA
# 4 NAs from 100% endemism. Final sample size is 126.

# Save result for downstream use
write.csv(env.complete,
          file = "output/tdwg3_predictors_complete.csv",
          eol = "\r\n"
          )
write.csv(env.complete.noch,
          file = "output/tdwg3_predictors_complete_noch.csv",
          eol = "\r\n"
          )



################################################
# Functions for applying Single-predictor models
################################################
cat("Loading functions for single-predictor analysis...\n")

# Wrapper function to run SingleModel() and process and save results
# ------------------------------------------------------------------
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


# Function to summarize single model outputs (1)
# ----------------------------------------------
CombineRealm <- function(statistic, fd.index, null.model, realm, filter.p = TRUE) {
# Filenames generated by RunMM have names with the following structure:
#   ("OLS_single_", statistic, "_", fd.index, "_", "null.model", "_", realm, ".csv")
#
# This is a base-level function to combine realm data and filter by p-value.
#
# Statistic: model statistic to summarize, one of c("Rsq", "slope", "p_value")
# fd.index: FD index to summarize, one of c("FRic", "FDis")
# null.model: nullmodel name to summarize
# realm: character vector giving realm names to summarize
# filter.p: whether to replace statistics from models where p > 0.05 with NA
  output <- vector("list", length = length(realm))
  names(output) <- realm
  for (i in seq_along(realm)) {
    # load requested statistic
    stat <- read.csv(file = paste0(output.dir, "OLS_single_", statistic, "_",
                                   fd.index, "_", null.model, "_", realm[i],
                                   ".csv")
                     )
     if (filter.p) {
       # load corresponding p values
       p.vals <- read.csv(file = paste0(output.dir, "OLS_single_p_value_",
                                        fd.index, "_", null.model, "_", realm[i],
                                        ".csv")
                          )
       # filter data for p < 0.05
       fd.list <- vector("list", length = nrow(stat))
       names(fd.list) <- rownames(stat)
       for (j in seq_len(nrow(stat))) {
         nonsignificant <- suppressWarnings(which(p.vals[j, ] > 0.05))
         # warnings suppressed because the 1st column is not numeric (originally
         # rownames) which is not in fact a problem
         stat[j, nonsignificant] <- "n.s."
       }
     }
    # Asssign (filtered) stat data to output list
    output[[i]] <- stat
  }
  # collapse output into a dataframe. Row and Col names should match.
  ldply(output, function(x) { x } , .id = "realm.subset")
}


# Function to summarize single model outputs (2)
# ----------------------------------------------
SingleSummary <- function(statistic, fd.index, null.models, realm,
                          filter.p = TRUE) {
# Filenames generated by RunMM have names with the following structure:
#   ("OLS_single_", statistic, "_", fd.index, "_", "null.model", "_", realm, ".csv")
#
# This function is meant to summarize files for different 'null.model' and 'realm'.
#
# Statistic: model statistic to summarize, one of c("Rsq", "slope", "p_value")
# fd.index: FD index to summarize, one of c("FRic", "FDis")
# null.models: character vector giving nullmodel names to summarize
# realm: character vector giving realm names to summarize
# filter.p: whether to replace statistics from models where p > 0.05 with NA

  # First, for each null.model, summarize over realm, combining results into a list
  null.model.list <- vector("list", length = length(null.models))
  names(null.model.list) <- null.models
  for (i in seq_along(null.model.list)) {
    null.model.list[[i]] <-
      CombineRealm(statistic, fd.index, null.models[i], realm, filter.p)
  }
  # Second, collapse resulting list into a dataframe
  ldply(null.model.list, function(x) { x } , .id = "null.model")
}


# Wrapper function to apply RunMM() to all data + to each realm subset
#---------------------------------------------------------------------
RunSingles <- function(fd.indices, fd.index, name, env.complete) {
# fd.indices: the list object 'fd.indices' or a subset thereof
# fd.index: FD index to analyze, one of c("FRic", "FDis")
# name: name of column to extract from each df in 'fd.indices' via GetFD(),
#       also passed to RunMM() as suffix for filenames
  # full dataset:
  fd.global <- GetFD(fd.indices, name)
  RunMM(fd.global,
        env.complete,
        paste0(fd.index, "_", name, "_all")
        )
  # realm subsets:
  fd.global.realms <- RealmSubset(fd.global)
  env.complete.realms <-
    RealmSubset(env.complete[, !colnames(env.complete) %in% "realm"])
  for (i in seq_along(fd.global.realms)) {
    RunMM(fd.global.realms[[i]],
          env.complete.realms[[i]],
          paste0(fd.index, "_", name, "_", names(env.complete.realms)[[i]])
          )
  }
  return (0)
}



##################################
# Generate single-predictor models
##################################
cat("Generating single-predictor OLS models:\n")
# Let's put those functions to work.

# Prepare names for easy subsetting
fric <- names(fd.indices)[1:4]
fdis <- names(fd.indices)[5:8]

# NOTE: standardization to mean 0 and unit variance is applied to predictors,
# making slopes comparable.
# For the predictor dataset we will use the 'noch' version, i.e. no canopy height.


# For functional richness (FRic)
# ------------------------------
cat("(1) FRic for null model Global...\n")
RunSingles(fd.indices[fric], "FRic", "global.SES", env.complete.noch)

cat("(2) FRic for null model Realm...\n")
RunSingles(fd.indices[fric], "FRic", "realm.SES", env.complete.noch)

cat("(3) FRic for null model Realm without MDG...\n")
RunSingles(fd.indices[fric], "FRic", "realm.SES.noMDG", env.complete.noch)

cat("(4) FRic for null model ADF...\n")
RunSingles(fd.indices[fric], "FRic", "adf.SES", env.complete.noch)

cat("(5) FRic for observed data...\n")
RunSingles(fd.indices[fric], "FRic", "observed", env.complete.noch)


# For functional dispersion (FDis)
# -------------------------------
cat("(6) FDis observed...\n")
RunSingles(fd.indices[fdis], "FDis", "observed", env.complete.noch)

cat("(7) FDis for null model Global...\n")
RunSingles(fd.indices[fdis], "FDis", "global.SES", env.complete.noch)

cat("(8) FDis for null model Realm...\n")
RunSingles(fd.indices[fdis], "FDis", "realm.SES", env.complete.noch)

cat("(9) FDis for null model Realm without MDG...\n")
RunSingles(fd.indices[fdis], "FDis", "realm.SES.noMDG", env.complete.noch)

cat("(10) FDis for null model ADF...\n")
RunSingles(fd.indices[fdis], "FDis", "adf.SES", env.complete.noch)


# Summarize results
# -----------------
cat("Summarizing single-predictor model results...\n")

null.models <- c("global.SES", "realm.SES", "realm.SES.noMDG", "adf.SES", "observed"
)
# For FRic
fric.rsq <- SingleSummary("Rsq",
                          "FRic",
                          null.models,
                          c("all", "NewWorld", "OWWest", "OWEast"),
                          filter.p = TRUE
                          )
write.csv(fric.rsq,
          file = paste0(output.dir, "00_OLS_single_Rsq_FRic_combined.csv"),
          row.names = FALSE,
          quote = FALSE
          )

# For FDis
fdis.rsq <- SingleSummary("Rsq",
                          "FDis",
                          null.models,
                          c("all", "NewWorld", "OWWest", "OWEast"),
                          filter.p = TRUE
                          )
write.csv(fdis.rsq,
          file = paste0(output.dir, "00_OLS_single_Rsq_FDis_combined.csv"),
          row.names = FALSE,
          quote = FALSE
          )



############################
# Investigate missing values
############################
cat("Investigating missing values...\n")
# NOTE: this is legacy code relating to the case where canopy height is included.
# Most NAs are due to NAs in canopy height.

# Which TDWG3 units are affected?
tdwg.observed <-
  fd.indices[[1]] [!is.na(fd.indices[[1]] [, "observed"]),
                   "observed",
                   drop = FALSE] %>%
  rownames()
tdwg.complete <- rownames(env.complete)[complete.cases(env.complete)]
affected <- tdwg.observed[!tdwg.observed %in% tdwg.complete]

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
                         factor.name =  "Complete data?",
                         title = "TDWG3 units with complete environmental data"
                         ),
       filename = paste0(output.dir, "environmental_data_missing_values.png"),
       width = 8,
       height = 4
       )
# Not a completely random distribution, but not too clumped either.



###############################################
# Functions for multi-predictor model selection
###############################################
cat("Loading functions for multi-predictor analysis...\n")

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
# provide a frame of reference.
#
# UPDATE: using the below performance metrics, you can verify that the exhaustive
# cross-validation model choice is indeed a good overall choice. 
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


# Function to apply mod selection automatically to all fd.indices
# ---------------------------------------------------------------
ApplyModSelect <- function(fd.indices, fd.names, predictors, null.models, k) {
# fd.indices: the object 'fd.indices' created above.
# fd.names: character vector giving names of the response variables, e.g. the
#           'fd.names' object created above, or a subset thereof.
# predictors: dataframe with predictor variables
# null.models: names of null models; or, technically, the column names in
#              fd.indices for which to run the analysis, e.g. "realm.SES".
# k: number of folds to use for cross-validation; passed to MultiSelect() and
#    from there to SelectOLS()

  # Init 3 output lists: for formulae, fitted models, and performance metrics
  fd.list <- vector("list", length = length(fd.names))
  names(fd.list) <- fd.names
  formulae <- rep(list(fd.list), length(null.models))
  names(formulae) <- null.models
  mods <- formulae
  performances <- formulae

  # Run selection process:
  for (model in null.models) {
    for (index in fd.names) {
      cat(model, index, "...\n")
      multimod <-
        MultiSelect(response.var = fd.indices[[index]] [, model],
                    response.name = index,
                    predictors = predictors,
                    k = k
                    )
      formulae[[model]] [[index]] <-
        multimod[["exhaustive"]] [["formulae"]] [["cross.validation"]]
      mods[[model]] [[index]] <-
        multimod[["exhaustive"]] [["mods"]] [["cross.validation"]]
      performances[[model]] [[index]] <- EvalMM(multimod)
    }
  }

  # Combine and return output
  list(formulae = formulae, mods = mods, performances = performances)
}


# Function to extract all Rsq from an ApplyModSelect() output object into a vector
# --------------------------------------------------------------------------------
GetBestRsq <- function(select) {
  num.ind <- length(names(select[["performances"]] [[1]]))
  num.mod <- length(names(select[["performances"]]))
  rsq <- numeric(length = num.ind * num.mod)
  i <- 1
  for (model in names(select[["performances"]])) {
    for (index in names(select[["performances"]] [[1]])) {
#      cat("model", model, "index", index, "product", i, "\n")
      rsq[i] <- select[["performances"]] [[model]] [[index]] [1, 2]
      # This extracts the data of row 1 ("best.mod") column 2 ("Rsq")
      names(rsq)[i] <- paste0(model, ".", index)
      i <- i + 1
    }
  }
  rsq
}


# Function to parse best model formulae for saving to disk
# --------------------------------------------------------
# Saving .csv files implies/requires coercion to dataframe format, so let's make
# a dataframe with the response variables, a single character column
# giving the formula of the best model, and a factor indicating the null model.
# While we're at it, throw in the best Rsq as well.
GetBestMods <- function(select) {
  resp.names <- names(select[["formulae"]] [[1]])
  num.resp <- length(resp.names)
  num.mods <- length(names(select[["formulae"]]))
  mods <- data.frame(null.model   = rep(names(select[["formulae"]]),
                                        each = num.resp
                                        ),
                     fd.index     = rep(resp.names, num.mods),
                     Rsq          = GetBestRsq(select),
                     best.formula = rep(NA, length = num.mods * num.resp)
                     )
  formulae <- character(length = num.mods * num.resp)
  for (i in seq_along(select[["formulae"]])) {
    formulae[(1 + num.resp * (i - 1)):(num.resp + num.resp * (i - 1))] <-
      ldply(select[["formulae"]] [[i]], as.character) [, 4]
  }
  mods[, "best.formula"] <- formulae
  mods 
}


# Wrapper function to run multi-model selection
# ---------------------------------------------
# Four predictor sets: all, eh, stability, meta (see below)
# Four cases: total data + 3 realm subsets
# That is sixteen cases in total. Generating 16 files is acceptable, but copying
# code 16 times is not. We need a wrapper function.
ApplyMMS <- function(fd.indices, fd.names, null.models, env.complete, predictors,
                     k, identifier) {
  # Subset env.complete to selected predictors:
  env.subset <- env.complete[, predictors]
  # Run for worldwide data:
  cat("For complete dataset:\n")
  write.csv(GetBestMods(ApplyModSelect(fd.indices,
                                       fd.names,
                                       env.subset,
                                       null.models,
                                       k = k
                                       )
                        ),
            file = paste0(output.dir, "OLS_best_models_", identifier, "_all.csv"),
            eol = "\r\n",
            row.names = FALSE
            )

  # Generate realm subset:
  env.subset.realms <- RealmSubset(env.subset[, !colnames(env.subset) %in% "realm"])
  # Subset fd.indices. This requires some magic because ApplyModSelect() needs
  # a clean subsetted list.
  # First, apply RealmSubset() to fd.indices:
  fd.indices.realms <- llply(fd.indices, RealmSubset)
  # This is a list of length x, with sublists of length y.
  # This must be transformed to a list of length y, with sublists of length x.
  fd.indices.transformed <- vector("list", length = length(fd.indices.realms[[1]]))
  names(fd.indices.transformed) <- names(fd.indices.realms[[1]])
  for (i in seq_along(fd.indices.transformed)) {
    fd.indices.transformed[[i]] <- llply(fd.indices.realms, function(x) { x[[i]] } )
  }

  # run for subsets:
  for (i in seq_along(env.subset.realms)) {
    cat("For data subset", names(env.subset.realms)[i], ":\n")
    write.csv(GetBestMods(ApplyModSelect(fd.indices.transformed[[i]],
                                         fd.names,
                                         env.subset.realms[[i]],
                                         null.models,
                                         k = k
                                         )
                          ),
              file = paste0(output.dir, "OLS_best_models_", identifier, "_",
                            names(env.subset.realms)[i], ".csv"),
              eol = "\r\n",
              row.names = FALSE
              )
  }
  return (0)
}



########################################
# Generating best multi-predictor models
########################################
cat("Generating best multi-predictor OLS models:\n")
# NOTE: excluding canopy height

# Preparation of input parameters for all runs:
null.models <- c("global.SES", "realm.SES", "realm.SES.noMDG", "adf.SES", "observed")
fd.names
k = 4

cat("(1) Predictors related to heterogeneity...\n")
ApplyMMS(fd.indices,
         fd.names,
         null.models,
         env.complete.noch,
         predictors = c("alt_range", "soilcount", "bio1_sd", "bio12_sd",
                        "bio1_mean", "bio12_mean"),
         k,
         identifier = "EH"
         )

cat("(2) Predictors related to stability...\n")
ApplyMMS(fd.indices,
         fd.names,
         null.models,
         env.complete.noch,
         predictors = c("bio4_mean", "bio15_mean", "lgm_Tano", "lgm_Pano"),
         k,
         identifier = "stability"
         )

cat("(3) Meta-predictors: palm richness, endemism, and realm...\n")
ApplyMMS(fd.indices,
         fd.names,
         null.models,
         env.complete.noch,
         predictors = c("palm.richness", "endemism", "realm"),
         k,
         identifier = "meta"
         )

cat("(4) All environmental predictors (EH + stability)...\n")
ApplyMMS(fd.indices,
         fd.names,
         null.models,
         env.complete.noch,
         predictors = c("alt_range", "soilcount", "bio1_sd", "bio12_sd",
                        "bio4_mean", "bio15_mean", "lgm_Tano", "lgm_Pano",
                        "bio1_mean", "bio12_mean"),
         k,
         identifier = "envir"
         )



cat("Done.\n")

