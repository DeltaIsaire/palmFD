###########################################
# Palm FD project: TDWG3 environmental data
###########################################
#
# In which environmental data at the resolution of botanical countries (TDWG3
# units) is prepared for regression analysis. This includes assessment of extreme
# values, testing for approximate normality and transforming if required.
#
# This 'new' script version focuses on selected key predictor values, thought to
# be relevant in light of our main hypothesis (that FD is positively associated 
# with greater environmental heterogeneity). As such we focus on key predictors
# measuring variability, rather then absolute values (e.g. variability in mean
# annual temperature, rather than said mean).
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


#########################
# Load environmental data
#########################
cat("Loading environmental data...\n")

# Basic info on tdwg3 units
# -------------------------
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")

# List of included tdwg3 units
tdwg3.included <-
  tdwg3.info[, "tdwg3.code"] %>%
  .[tdwg3.info[, "experimental.richness"] > 0] %>%
  .[complete.cases(.)] %>%
  as.character()



# Environmental data
# ------------------
env.new <- read.csv(file = "data/TDWG_Environment_AllData_2019Feb.csv")


###############################
# Subset and process predictors
###############################
cat("Processing evironmental data...\n")

# Non-climatic environmental
# --------------------------
# Desired biogeographical heterogeneity predictors from 'env.new' include:
#   alt_range  (altitude range)
#   soilcount  (number of different soil types present)
#   CH_SD      (standard deviation of canopy height)
env.nonclim <- env.new[, c("LEVEL_3_CO", "alt_range", "soilcount", "CH_SD")]
# For appropriate evaluation, subset the data to datapoints we are actually
# investigating, i.e. to those tdwg3 units included in our study.
env.nonclim %<>% .[.[, "LEVEL_3_CO"] %in% tdwg3.included, ]
#
# CH_SD has 14 missing values. This was investigated in the old script, and deemed
# not problematic.
#
# Distributions visually examined with Histogram().
# The following transformations were deemed necessary to improve normality:
env.nonclim[, "alt_range"] %<>% sqrt()


# Climatic heterogeneity
# ----------------------
# Desired climatic heterogeneity predictors from 'env.new' include:
#   bio1_sd   (standard deviation of annual temperature)
#   bio12_sd  (standard deviation of annual precipitation)
env.clim.h <-
  env.new[, c("LEVEL_3_CO", "bio1_sd", "bio12_sd")] %>%
  .[.[, "LEVEL_3_CO"] %in% tdwg3.included, ]
#
# No missing values.
#
# Distributions visually examined with Histogram().
# The following transformations were deemed necessary to improve normality:
env.clim.h[, "bio1_sd"] %<>% log10()
env.clim.h[, "bio12_sd"] %<>% sqrt()
#
# Transformed bio1_sd has 1 low outlier, for TCI (Turks-Caicos Island)


# Climatic stability
# ------------------
# Desired climatic stability predictors from 'env.new' include:
#   "bio4_mean"   (Temperature seasonality)
#   'bio15_mean"  (Precipitation seasonality)
#   LGM temperature anomaly
#   LGM precipitation anomaly
# The LGM anomalies need to be calculated from LGM mean and contemporary mean values
env.clim.stable <-
  data.frame(env.new[, c("LEVEL_3_CO", "bio4_mean", "bio15_mean")],
             lgm_Tano = (env.new[, "lgm_ens_Tmean"] / 10) - env.new[, "bio1_mean"],
             lgm_Pano = (env.new[, "lgm_ens_Pmean"] / 10) - env.new[, "bio12_mean"]
             ) %>%
  .[.$LEVEL_3_CO %in% tdwg3.included, ]
# NOTE: The new LGM data for FIJ (Fiji) is faulty due to calculation glitches.
env.clim.stable[env.clim.stable[, "LEVEL_3_CO"] == "FIJ",
                c("lgm_Tano", "lgm_Pano")
                ] <- NA
#
# 1 missing value for the LGM anomalies, per the above.
#
# Distributions visually examined with Histogram().
# The following transformations were deemed necessary to improve normality:
env.clim.stable[ "bio4_mean"] %<>% log10()
#
# LGM temperature anomaly has two extremely high values. These are CHT (Tibet)
# and WHM (West Himalaya).
# LGM Precipitation anomaly has three extremely low values. These are
# SCZ (Santa Cruz Is.), BIS (Bismarck Archipelago), CRL (Caroline Is.).



#####################################
# Gather and save predictor variables
#####################################
cat("Merging and saving predictor dataframes...\n")

# Gather predictors,
# and pad rows (with NA values) to include all TDWG3 units, for easier
# downstream coding.
predictors <-
  merge(env.nonclim, env.clim.h, by = "LEVEL_3_CO") %>%
  merge(., env.clim.stable, by = "LEVEL_3_CO") %>%
  merge(.,
        tdwg3.info[, "tdwg3.code", drop = FALSE],
        by.x = "LEVEL_3_CO",
        by.y = "tdwg3.code",
        all.y = TRUE,
        sort = FALSE
        )
rownames(predictors) <- predictors[, "LEVEL_3_CO"]
predictors %<>% .[, -1]


# Save data
write.csv(predictors,
          file = "output/tdwg3_environmental_predictors_new.csv",
          eol = "\r\n",
          row.names = TRUE
          )


cat("Done.\n")

