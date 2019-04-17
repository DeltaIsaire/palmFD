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
# ---------------------------------------
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

# Account for NAs in predictors. Subsetting to complete cases would exclude
# tdwg3.units, which hinders data merging. Instead, propagate NAs across rows.
env.complete <- env
env.complete[!complete.cases(env.complete), ] <- NA



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


# Function to extract data of each null model from fd.indices
# -----------------------------------------------------------
GetFD <- function(fd.indices, name) {
#  fd.indices: list of FD index dataframes (generated above)
#  name: name of column to extract from each df in 'fd.indices'
  fd.list <- llply(fd.indices, function(x) { x[, name] } )
  df <- as.data.frame(fd.list)
  rownames(df) <- rownames(fd.indices[[1]])
  df
}


# Function to subset the input data to realm, combining results into a list
# -------------------------------------------------------------------------
RealmSubset <- function(x) {
# Where x is a dataframe, such as the output of GetFD(), with tdwg3 codes as
# rownames
  output <- vector("list", length = length(realm.tdwg3))
  names(output) <- names(realm.tdwg3)
  for (i in seq_along(output)) {
    output[[i]] <- x[rownames(x) %in% realm.tdwg3[[i]], ]
  }
  output
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
         stat[j, nonsignificant] <- NA
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
cat("Generating single-predictor models:\n")
# Let's put those functions to work

# Prepare names for easy subsetting
fric <- names(fd.indices)[1:4]
fdis <- names(fd.indices)[5:8]

# NOTE: standardization to mean 0 and unit variance is applied to predictors,
# making slopes comparable.


# For functional richness (FRic)
# ------------------------------
cat("(1) FRic for null model Global...\n")
RunSingles(fd.indices[fric], "FRic", "global.SES", env.complete)

cat("(2) FRic for null model Realm...\n")
RunSingles(fd.indices[fric], "FRic", "realm.SES", env.complete)

cat("(3) FRic for null model ADF...\n")
RunSingles(fd.indices[fric], "FRic", "adf.SES", env.complete)


# For functional dispersion (FDis)
# -------------------------------
cat("(4) FDis observed...\n")
RunSingles(fd.indices[fdis], "FDis", "observed", env.complete)


# Summarize results
# -----------------
cat("Summarizing single-predictor model results...\n")

# For FRic
fric.rsq <- SingleSummary("Rsq",
                          "FRic",
                          c("global.SES", "realm.SES", "adf.SES"),
                          c("all", "NewWorld", "OWWest", "OWEast"),
                          filter.p = TRUE
                          )
write.csv(fric.rsq,
          file = paste0(output.dir, "OLS_single_Rsq_FRic_combined.csv"),
          row.names = FALSE
          )

# For FDis
fdis.rsq <- SingleSummary("Rsq",
                          "FDis",
                          "observed",
                          c("NewWorld", "OWWest", "OWEast"),
                          filter.p = TRUE
                          )
write.csv(fdis.rsq,
          file = paste0(output.dir, "OLS_single_Rsq_FDis_combined.csv"),
          row.names = FALSE
          )



############################
# Investigate missing values
############################
cat("Investigating missing values...\n")

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



#################################
# Multi-predictor model selection
#################################
cat("Selecting best multi-predictor models:\n")

# Predictor selection
# -------------------
# Predictor preprocessing showed that some predictors cannot be used together
# in the same model, due to extreme collinearity.
# Here, these predictors are removed. Endemism is also removed because it is
# zero-inflated.
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
ApplyModSelect <- function(fd.indices, fd.names, predictors) {
# fd.indices: the object 'fd.indices' created above.
# fd.names: character vector giving names of the response variables
# predictors: dataframe with predictor variables

  # Init 3 output lists: for formulae, fitted models, and performance metrics
  fd.list <- vector("list", length = length(fd.names))
  names(fd.list) <- fd.names
  null.models <- c("global", "realm", "adf")
  formulae <- rep(list(fd.list), 3)
  names(formulae) <- null.models
  mods <- formulae
  performances <- formulae

  # Run selection process:
  for (model in null.models) {
    for (index in fd.names) {
      cat(model, index, "...\n")
      multimod <-
        MultiSelect(response.var = fd.indices[[index]] [, paste0(model, ".SES")],
                    response.name = index,
                    predictors = predictors
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


# Apply model selection for all response variables
# ------------------------------------------------
cat("Applying model selection for all response variables:\n")
# Using either (1) bioclim data, or (2) climate pca axes.
# Additionally, either including or excluding canopy height (which has 14 NAs).
# Additionally, either including or excluding c("realm", "palm.richness")
# So that's 8 different runs on 3 * 8 = 24 response variables...
#
cat("(1) climatic variables + realm + palm.richness + CH...\n")
select.one <- ApplyModSelect(fd.indices, fd.names, env.subset.clim)

cat("(2) climatic variables + realm + palm.richness, excluding CH...\n")
env.subset.clim.noCH <-
  env.subset.clim[, !names(env.subset.clim) %in% c("CH_Mean", "CH_Range")]
select.two <- ApplyModSelect(fd.indices, fd.names, env.subset.clim.noCH)

cat("(3) climatic variables + CH, excluding realm + palm.richness...\n")
env.subset.clim.noRealm <-
  env.subset.clim[, !names(env.subset.clim) %in% c("realm", "palm.richness")]
select.three <- ApplyModSelect(fd.indices, fd.names, env.subset.clim.noRealm)

cat("(4) climatic variables, excluding realm + palm.richness + CH...\n")
env.subset.clim.none <-
  env.subset.clim.noCH[, !names(env.subset.clim.noCH) %in% c("realm",
                                                             "palm.richness"
                                                             )
                       ]
select.four <- ApplyModSelect(fd.indices, fd.names, env.subset.clim.none)


cat("(5) climate PCA axes + realm + palm.richness + CH...\n")
select.five <- ApplyModSelect(fd.indices, fd.names, env.subset.pca)

cat("(6) climate PCA axes + realm + palm.richness, excluding CH...\n")
env.subset.pca.noCH <-
  env.subset.pca[, !names(env.subset.pca) %in% c("CH_Mean", "CH_Range")]
select.six <- ApplyModSelect(fd.indices, fd.names, env.subset.pca.noCH)

cat("(7) climate PCA axes + CH, excluding realm + palm.richness...\n")
env.subset.pca.noRealm <-
  env.subset.pca[, !names(env.subset.pca) %in% c("realm", "palm.richness")]
select.seven <- ApplyModSelect(fd.indices, fd.names, env.subset.pca.noRealm)

cat("(8) climate PCA axes, excluding realm + palm.richness + CH...\n")
env.subset.pca.none <-
  env.subset.pca.noCH[, !names(env.subset.pca.noCH) %in% c("realm",
                                                           "palm.richness"
                                                             )
                       ]
select.eight <- ApplyModSelect(fd.indices, fd.names, env.subset.pca.none)


# Parse and save Rsq data
# -----------------------
# What we want is a dataframe with Rsq for all response variables (rows) and all
# 8 predictor datasets (columns).

# init matrix
Rsq <-
  matrix(data = NA,
         nrow = 8 * 3,
         ncol = 8,
         dimnames = list(rep(fd.names, 3),
                         c("clim.all", "clim.noCH", "clim.norealm", "clim.none",
                           "pca.all", "pca.noCH", "pca.norealm", "pca.none"
                           )
                         )
         )
# Function to extract all Rsq from an output object into a vector
GetRsq <- function(select) {
  rsq <- numeric(length = 8 * 3)
  i <- 1
  for (model in 1:3) {
    for (index in 1:8) {
#      cat("model", model, "index", index, "product", i, "\n")
      rsq[i] <- select[["performances"]] [[model]] [[index]] [1, 2]
      i <- i + 1
    }
  }
  rsq
}

# Gather data and save
rsq <- data.frame(null.model       = rep(c("global", "realm", "adf"), each = 8),
                  fd.index         = rep(fd.names, 3),
                  clim.all         = GetRsq(select.one),
                  clim.noCH        = GetRsq(select.two),
                  clim.norealm     = GetRsq(select.three),
                  clim.none        = GetRsq(select.four),
                  clim.pca.all     = GetRsq(select.five),
                  clim.pca.noCH    = GetRsq(select.six),
                  clim.pca.norealm = GetRsq(select.seven),
                  clim.pca.none    = GetRsq(select.eight)
                  )
write.csv(rsq,
          file = paste0(output.dir, "OLS_best_models_rsq.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

# Parse and save best model formulae
# ----------------------------------
# With focus on 'select.one' because Rsq analysis indicates this run contains the
# best combination of predictors.

# We'll write a function to extract the formulae.
# Saving .csv files implies/requires coercion to dataframe format, so let's make
# a dataframe with the response variables, a single character column
# giving the formula of the best model, and a factor indicating the null model.
GetBestMods <- function(select) {
  pred.names <- names(select[["formulae"]] [[1]])
  mods <- data.frame(null.model = rep(names(select[["formulae"]]),
                                      each = length(pred.names)
                                      ),
                     fd.index = rep(pred.names, 3),
                     best.formula = rep(NA, length = 3 * length(pred.names))
                     )
  formulae <- character(length = 3 * length(pred.names))
  for (i in seq_along(select[["formulae"]])) {
    formulae[(1 + 8 * (i - 1)):(8 + 8 * (i - 1))] <-
      ldply(select[["formulae"]] [[i]], as.character) [, 4]
  }
  mods[, "best.formula"] <- formulae
  mods 
}

# apply function and save:
best.mods <- GetBestMods(select.one)

write.csv(best.mods,
          file = paste0(output.dir, "OLS_best_models_predictors.csv"),
          eol = "\r\n",
          row.names = FALSE
          )


#############################################
# Multi-Models related to specific hypotheses
#############################################
cat("Fitting and evaluating models for specific hypotheses:\n")
# Essentially, here we focus on narrow subsets of predictor variables, to test
# how each group of predictors relates to FD.

cat("(1) Non-climatic environmental variation...\n")
# Select only non-clim environmental predictors:
env.nonclim <- env.subset[, c("alt_range", "soilcount", "CH_Mean", "CH_Range")]
# Apply best model selection:
select.nonclim <- ApplyModSelect(fd.indices, fd.names, env.nonclim)
# Extract and save best model formulae
write.csv(GetBestMods(select.nonclim),
          file = paste0(output.dir, "OLS_best_models_predictors_nonclim.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

cat("(2) Climate change anomaly since LGM...\n")
# Select only LGM climate anomaly predictor:
env.lgm <- env.subset[, c("lgm_Tano", "lgm_Pano")]
# Apply best model selection:
select.lgm <- ApplyModSelect(fd.indices, fd.names, env.lgm)
# Extract and save best model formulae
write.csv(GetBestMods(select.lgm),
          file = paste0(output.dir, "OLS_best_models_predictors_lgm.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

cat("(3) Contemporary climate:\n")
cat("(3a) Bioclim variables...\n")
env.clim <- env.subset[, c("bio1_mean", "bio2_mean", "bio3_mean", "bio4_mean",
                           "bio7_mean", "bio8_mean", "bio9_mean", "bio10_mean",
                           "bio11_mean", "bio12_mean", "bio13_mean", "bio14_mean",
                           "bio15_mean", "bio16_mean", "bio17_mean", "bio18_mean",
                           "bio19_mean"
                           )
                       ]
# Apply best model selection:
select.clim <- ApplyModSelect(fd.indices, fd.names, env.clim)
# Extract and save best model formulae
write.csv(GetBestMods(select.clim),
          file = paste0(output.dir, "OLS_best_models_predictors_clim.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

cat("(3b) Temperature and Precipitation seasonality...\n")
env.seas <- env.subset[, c("bio4_mean", "bio15_mean")]
# Apply best model selection:
select.seas <- ApplyModSelect(fd.indices, fd.names, env.seas)
# Extract and save best model formulae
write.csv(GetBestMods(select.seas),
          file = paste0(output.dir, "OLS_best_models_predictors_seas.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

cat("(3c) Contemporary climate PCA axes...\n")
env.clim.pca <- env.subset[, c("clim.PC1", "clim.PC2", "clim.PC3", "clim.PC4")]
# Apply best model selection:
select.clim.pca <- ApplyModSelect(fd.indices, fd.names, env.clim.pca)
# Extract and save best model formulae
write.csv(GetBestMods(select.clim.pca),
          file = paste0(output.dir, "OLS_best_models_predictors_clim_pca.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

cat("(4) Biogeographical parameters...\n")
env.geo <- env.subset[, c("realm", "palm.richness")]
# Apply best model selection:
select.geo <- ApplyModSelect(fd.indices, fd.names, env.geo)
# Extract and save best model formulae
write.csv(GetBestMods(select.geo),
          file = paste0(output.dir, "OLS_best_models_predictors_geo.csv"),
          eol = "\r\n",
          row.names = FALSE
          )

# Extract and save best model Rsq for all hypothesis cases
# --------------------------------------------------------
rsq2 <- data.frame(null.model = rep(c("global", "realm", "adf"), each = 8),
                   fd.index   = rep(fd.names, 3),
                   nonclim    = GetRsq(select.nonclim),
                   lgm        = GetRsq(select.lgm),
                   clim       = GetRsq(select.clim),
                   clim.pca   = GetRsq(select.clim.pca),
                   clim.seas  = GetRsq(select.seas),
                   geo        = GetRsq(select.geo)
                   )
write.csv(rsq2,
          file = paste0(output.dir, "OLS_best_models_rsq_hypotheses.csv"),
          eol = "\r\n",
          row.names = FALSE
          )





cat("Done.\n")

