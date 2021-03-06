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

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/OLS_regression_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")

# Model output directory (with trailing slash!)
output.dir <- "output/OLS_models/"

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

# Additionally, merge richness, endemism and realm into env
indices <- match(rownames(fd.indices[[1]]), rownames(env))
env[, "palm.richness"] <- fd.indices[[1]] [indices, "palm.richness"]
env[, "endemism"] <- fd.indices[[1]] [indices, "percent.endemism"]
indices <- match(rownames(env), tdwg3.info[, "tdwg3.code"])
env[, "realm"] <- tdwg3.info[indices, "realm"]
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

# Function to subset the input data to realm, combining results into a list
# -------------------------------------------------------------------------
RealmSubset <- function(x) {
  output <- vector("list", length = length(realm.tdwg3))
  names(output) <- names(realm.tdwg3)
  for (i in seq_along(output)) {
    output[[i]] <- x[rownames(x) %in% realm.tdwg3[[i]], ]
  }
  output
}


# Call the functions for each null model dataset
# ----------------------------------------------
# In addition to the full dataset, it is good to look at subsets for each of the
# three realms. To avoid combinatorial explosion, for realm subsets we only look
# at the case where Canopy Height is INcluded. Results on the full data indicate
# that including CH does not strongly affect the results.
cat("(1) For null model Global...\n")
# full dataset:
fd.global <- GetFD("global")
RunMM(fd.global, env, "global")
RunMM(fd.global, env.noCH, "global_noCH")
# realm subsets:
fd.global.realms <- RealmSubset(fd.global)
env.realms <- RealmSubset(env[, !colnames(env) %in% "realm"])
for (i in 1:3) {
  RunMM(fd.global.realms[[i]],
        env.realms[[i]],
        paste0("global_", names(env.realms)[[i]])
        )
}

cat("(2) For null model Realm...\n")
fd.realm <- GetFD("realm")
RunMM(fd.realm, env, "realm")
RunMM(fd.realm, env.noCH, "realm_noCH")
fd.realm.realms <- RealmSubset(fd.realm)
for (i in 1:3) {
  RunMM(fd.realm.realms[[i]],
        env.realms[[i]],
        paste0("realm_", names(env.realms)[[i]])
        )
}

cat("(3) For null model ADF...\n")
fd.adf <- GetFD("adf")
RunMM(fd.adf, env, "adf")
RunMM(fd.adf, env.noCH, "adf_noCH")
fd.adf.realms <- RealmSubset(fd.adf)
for (i in 1:3) {
  RunMM(fd.adf.realms[[i]],
        env.realms[[i]],
        paste0("adf_", names(env.realms)[[i]])
        )
}


# Summarize results
# -----------------
cat("Summarizing single-predictor model results...\n")
SingleSummary <- function(ch, name.ext = "") {
# 1. For each null model, for each of the 8 FD indices, get the R-squared of
# statistically significant predictors, combining results into a list.
#
# ch: logical indicating whether to use the data INcluding canopy height
# name.ext: name extension of the datafiles (e.g. in 'global_NewWorld' the name
# extension is '_NewWorld')
  canopy <- ifelse(ch, "", "_noCH")
  single.predictors <- vector("list", length = 3)
  names(single.predictors) <-
    c("global", "realm", "adf") %>%
    paste0(., name.ext)
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

# Apply for full dataset:
results.single <- SingleSummary(ch = TRUE)
results.single.noCH <- SingleSummary(ch = FALSE)
# Apply for realm subsets:
results.single.NewWorld <- SingleSummary(ch = TRUE, name.ext = "_NewWorld")
results.single.OWWest <- SingleSummary(ch = TRUE, name.ext = "_OWWest")
results.single.OWEast <- SingleSummary(ch = TRUE, name.ext = "_OWEast")


# Function to parse and save the Rsq of significant predictors
# ------------------------------------------------------------
SaveResults <- function(results, file) {
  GetNames <- function(x) {
    unique(unlist(llply(x, names)))
  }

  cols <- unique(unlist(llply(results, GetNames)))

  Parse <- function(x) {
    rows <- names(x)
    mat <- matrix(data = NA,
                  nrow = length(rows),
                  ncol = length(cols),
                  dimnames = list(rows, cols)
                  )
    for (i in seq_along(x)) {
      indices <- match(names(x[[i]]), colnames(mat))
      mat[i, indices] <- suppressWarnings(as.numeric(x[[i]]))
      # warning suppression because x[[i]] could be of class chr instead of df.
      # That is not an error, so the warning can be ignored
    }
    mat
  }
  output.list <- llply(results, Parse)
  output <- output.list[[1]]
  if (length(output.list) > 1) {
    for (i in 2:length(output.list)) {
      output <- rbind(output, output.list[[i]])
    }
  }
  output %<>% as.data.frame()
  output <- data.frame(null.model = rep(names(results),
                                        each = nrow(output) / length(results)
                                        ),
                       fd.index   = rownames(output),
                       output,
                       row.names = seq_len(nrow(output))
                       )
  write.csv(output, file = file, eol = "\r\n", row.names = FALSE)
  output
}

# Call the function
# -----------------
SaveResults(results.single,
            file = paste0(output.dir, "OLS_single_Rsq_all.csv")
            )
SaveResults(results.single.noCH,
            file = paste0(output.dir, "OLS_single_Rsq_all_noCH.csv")
            )
SaveResults(results.single.NewWorld,
            file = paste0(output.dir, "OLS_single_Rsq_all_NewWorld.csv")
            )
SaveResults(results.single.OWWest,
            file = paste0(output.dir, "OLS_single_Rsq_all_OWWest.csv")
            )
SaveResults(results.single.OWEast,
            file = paste0(output.dir, "OLS_single_Rsq_all_OWEast.csv")
            )



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

