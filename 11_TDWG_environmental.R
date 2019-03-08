###########################################
# Palm FD project: TDWG3 environmental data
###########################################
#
# In which environmental data at the resolution of botanical countries (TDWG3
# units) is prepared for regression analysis. This includes assessment of extreme
# values, testing for approximate normality and transforming if required, and
# standardizing the (transformed) predictors to mean 0 and unit variance.
#
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(vegan)
library(corrplot)
library(car)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


#####################################
# Load and process environmental data
#####################################
cat("Loading and processing environmental data...\n")

# First, load included TDWG3 units, for subsetting
# ------------------------------------------------
tdwg3.included <-
  read.csv(file = "output/null_model_global/FD_z_scores_global_genus_mean_z.scores.csv",
           row.names = 1
           ) %>%
  rownames()

# Basic info on tdwg3 units
# -------------------------
tdwg3.info <- read.csv(file = "output/tdwg3_info.csv")
# The column 'palm.richness' reflects the full known palm richness.
# However, our data analysis uses a subset of all palm species, so the experimental
# richness values are lower.
# Add those as a new column:
pres.abs.matrix <- 
  read.csv(file = "output/palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
experimental.richness <-
  merge(x = tdwg3.info[, c("tdwg3.code", "palm.richness")],
        y = data.frame(tdwg3.code = rownames(pres.abs.matrix),
                       richness   = rowSums(pres.abs.matrix)
                       ),
        by = "tdwg3.code",
        all.x = TRUE
        ) %>%
  .[, "richness"]
tdwg3.info$experimental.richness <- experimental.richness

# Environmental data
# ------------------
env.old <- read.csv(file = "data/TDWG_Environment_AllData_2019Jan.csv")
env.new <- read.csv(file = "data/TDWG_Environment_AllData_2019Feb.csv")


# Non-climatic environmental
# --------------------------
env.nonclim <- 
  env.new[, c("LEVEL_3_CO", "srtm_alt_mean", "alt_range", "soilcount",
                           "CH_Mean", "CH_Range")] %>%
  .[.$LEVEL_3_CO %in% tdwg3.included, ]
# This data is the same as the old data.
#
# srtm_alt_mean is highly left-skewed.
# alt_range is also left-skewed.
# soilcount is slightly left-skewed, and too uniform to be normal.
# mean canopy height is decently normal.
# canopy height range is right-skewed.
#
# None of these variables have clear outliers or faulty values.
#
# Transform the skewed variables:
env.nonclim[, "srtm_alt_mean"] %<>% log()
env.nonclim[, "alt_range"] %<>% sqrt()


# Contemporary climate
# --------------------
env.clim <-
  env.new[, c(1, 75:93)] %>%
  .[.$LEVEL_3_CO %in% tdwg3.included, ]
# Not the same data as the 'old' variables, but each of the old variables matches
# pretty well with the new bioclim data, so we can safely use the new data.
#
# bio1 has 2 extreme values on the low end, and is right-skewed.
#   These extremes are CHT (Tibet) and WHM (West Himalaya).
# bio2 is somewhat bimodal.
# bio3 is slightly right-skewed.
# bio4 is left-skewed and three-modal. This is temperature seasonality, which is
#   likely correlated with latitude. Latitudinal variation in seasonality does
#   make sense though.
# bio5 has 1-2 extreme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West Himalaya).
# bio6 has 2 extreme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West Himalaya).
# bio7 is left-skewed.
# bio8 has 2-3 extreme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West Himalaya), and ALA (Alabama, USA).
# bio9 has 1 extreme value on the low end, and is right-skewed.
#   This extreme value is CHT (Tibet)
# bio10 has 2 extreme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West Himalaya).
# bio11 has 2 exterme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West Himalaya).
# bio12 has 1 extreme value on the high end.
#   This extreme value is SCZ (Santa Cruz Is.).
# bio13 is slightly skewed.
# bio14 is highly right-skewed, and has 1 extreme value on the high end.
#   this extreme value is SCZ (santa Cruz Is.).
# bio15 is fairly normal, but with ~2 rather high values. This is precipitation
#   seasonality.
#   The extreme values are CHA (Chad) and NGR (Niger).
# bio16 tends towards uniformity.
# bio17 is highly left-skewed, with 1 extreme value on the high end.
#   This extreme value is SCZ (Santa Cruz Is.)
# bio18 is left-skewed.
# bio19 is strongly left-skewed, with 2-3 extreme values on the high end.
#   These extreme values are SIE (Sierra Leone), BIS (Bismarck Archipelago),
#   and SCZ (Santa Cruz Is.).
#
# There is likely a strong spatial structure in all these variables.
# Most extreme values are linked to CHT (Tibet), WHM (West Himalaya),
# or SCZ (Santa Cruz Is.)
# Likely all extremes are real values.
#
# Transform the skewed variables:
env.clim[, "bio1_mean"] %<>% . ^ 2
env.clim[, "bio4_mean"] %<>% log()
env.clim[, "bio9_mean"] %<>% . ^ 2
env.clim[, "bio14_mean"] %<>% sqrt()
env.clim[, "bio17_mean"] %<>% sqrt()
env.clim[, "bio18_mean"] %<>% sqrt()
env.clim[, "bio19_mean"] %<>% sqrt()


# LGM climate (anomaly)
# ---------------------
env.lgm <-
  data.frame(env.new[, c("LEVEL_3_CO", "lgm_ens_Tmean", "lgm_ens_Pmean")],
             lgm_Tano = (env.new[, "lgm_ens_Tmean"] / 10) - env.new[, "bio1_mean"],
             lgm_Pano = (env.new[, "lgm_ens_Pmean"] / 10) - env.new[, "bio12_mean"]
             ) %>%
  .[.$LEVEL_3_CO %in% tdwg3.included, ]
# The old data was based only on CCSM and MIROC, while the new data
# includes more sources. So new should be 'better', and some deviations are
# to be expected, which is OK.
# The new data is in K * 10 while the old data is in K.
# They do however have the same SD, which is mildly weird.
# NOTE: The new LGM data for FIJ (Fiji) is faulty due to calculation glitches.
env.lgm[match("FIJ", env.lgm$LEVEL_3_CO), -1] <- NA
#
# LGM mean temperature is right-skewed, with 2 extreme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West himalaya), which makes sense.
# LGM Mean precipitation is slightly bimodal but otherwise fine.
# LGM temperature anomaly has two extremely high values. These are CHT (Tibet)
# and WHM (West Himalaya).
# LGM Precipitation anomaly has three extremely low values. These are
# SCZ (Santa Cruz Is.), BIS (Bismarck Archipelago), CRL (Caroline Is.).
#
# No transformations necessary. These outliers are a little troubling however.


##############################################################
# Gather and standardize the (transformed) predictor variables
##############################################################
cat("Standardizing and saving predictor variables...\n")
# Note that standardized values here are mostly for reference.
# In downstream analyses, certain TDWG3 units may be excluded, at which point
# standardized values would have to be recalculated.

# Gather predictors and pad rows to include all TDWG3 units
predictors <-
  merge(env.nonclim, env.clim, by = "LEVEL_3_CO") %>%
  merge(., env.lgm, by = "LEVEL_3_CO") %>%
  merge(.,
        tdwg3.info[, "tdwg3.code", drop = FALSE],
        by.x = "LEVEL_3_CO",
        by.y = "tdwg3.code",
        all.y = TRUE,
        sort = TRUE
        )
rownames(predictors) <- predictors[, "LEVEL_3_CO"]
predictors %<>% .[, -1]

# Standardize all predictors
predictors.std <- apply(predictors, 2, scale, center = TRUE, scale = TRUE)
rownames(predictors.std) <- rownames(predictors)

# save both versions
write.csv(predictors,
          file = "output/tdwg3_environmental_predictors.csv",
          eol = "\r\n",
          row.names = TRUE
          )
write.csv(predictors.std,
          file = "output/tdwg3_environmental_predictors_standardized.csv",
          eol = "\r\n",
          row.names = TRUE
          )


#################################
# Exploratory PCA and Correlogram
#################################
cat("Performing exploratory PCA...\n")

# Principal Components Analysis
# -----------------------------
# We will use PCA based on a correlation matrix, i.e. with standardization of 
# variables. This is suitable for environmental variables.

pca.env <- rda(predictors[complete.cases(predictors), ], scale = TRUE)
summary(pca.env)
# >95% of variation is explained with 9 axes, which is a rough indication of how many
# predictors could be important, at most, in a full model.
# OrdinationPlot(pca.env)


# Correlations between predictors
# -------------------------------
cat("Creating correlation plot of predictors...\n")

corr.matrix <- 
  predictors[complete.cases(predictors), ] %>%
  apply(., 2, scale, center = TRUE, scale = TRUE) %>%
  cor()

round(corr.matrix, digits = 2)
# Predictably, this is a mess. Going through it is entirely possible, but it is
# easier to create a correlogram.
GraphSVG(corrplot(corr.matrix, method = "square", order = "FPC"),
         file = "graphs/environmental_predictors_corrplot.svg",
         width = 8,
         height = 8
         )


cat("Done.\n")

