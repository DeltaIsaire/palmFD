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
env.new <- read.csv(file = "data/TDWG_Environment_AllData_20May2019.csv")



###############################
# Subset and process predictors
###############################
cat("Processing evironmental data...\n")

# Non-climatic environmental
# --------------------------
# Desired biogeographical heterogeneity predictors from 'env.new' include:
#   alt_range  (altitude range) (meters)
#   soilcount  (number of different soil types present) (count)
#   CH_SD      (standard deviation of canopy height) (meters)
env.nonclim <- env.new[, c("LEVEL_3_CO", "alt_range", "soilcount", "CH_SD")]
# For appropriate evaluation, subset the data to datapoints we are actually
# investigating, i.e. to those tdwg3 units included in our study.
env.nonclim %<>% .[.[, "LEVEL_3_CO"] %in% tdwg3.included, ]
#
# CH_SD has 14 missing values.
#
# Distributions visually examined with Histogram().


# Climatic heterogeneity
# ----------------------
# Desired climatic heterogeneity predictors from 'env.new' include:
#   bio1_sd   (standard deviation of annual temperature) (degrees * 10)
#   bio12_sd  (standard deviation of annual precipitation) (mm)
env.clim.h <-
  env.new[, c("LEVEL_3_CO", "bio1_sd", "bio12_sd")] %>%
  .[.[, "LEVEL_3_CO"] %in% tdwg3.included, ]
#
# No missing values.
#
# Distributions visually examined with Histogram().
# bio1_sd should be converted to degrees
env.clim.h[, "bio1_sd"] %<>% { . / 10 }
#


# Update: we should also include mean temperature and precipitation
#   bio1_mean (degrees celsius *10)
#   bio12_mean (mm)
env.clim <-
  env.new[, c("LEVEL_3_CO", "bio1_mean", "bio12_mean")] %>%
  .[.[, "LEVEL_3_CO"] %in% tdwg3.included, ]
# Convert bio1_mean to Celsius:
env.clim$bio1_mean %<>% { . / 10 }
#
# No missing values.
#
# Distributions visually examined with Histogram().
# Transforming is not going to do us much good.


# Climatic stability
# ------------------
# Desired climatic stability predictors from 'env.new' include:
#   "bio4_mean"   (Temperature seasonality) (degrees * 100)
#   'bio15_mean"  (Precipitation seasonality) (unitless, coefficient of variation)
#   LGM temperature anomaly (degrees C)
#   LGM precipitation anomaly (mm/year)
#   'PLIO_T_ANOM' (pliocene temp anomaly) (degrees C * 10)
#   'PLIO_P_ANOM' (pliocene prec anomaly) (mm/year)
#   Miocene temperature anomaly (degrees C)
#   Miocene precipitation anomaly (mm/year)
# The LGM anomalies need to be calculated from LGM mean and contemporary mean values
# PLIO_T_ANOM and miocene_temp_anom_mean should be divided by 10 to convert to C.
env.clim.stable <-
  data.frame(env.new[, c("LEVEL_3_CO", "bio4_mean", "bio15_mean")],
             lgm_Tano  = (env.new[, "bio1_mean"] / 10) -
                         (env.new[, "lgm_ens_Tmean"] / 10 - 273.15),
             lgm_Pano  = env.new[, "bio12_mean"] -
                         (env.new[, "lgm_ens_Pmean"] / 10),
             plio_Tano = env.new[, "PLIO_T_ANOM"] / 10,
             plio_Pano = env.new[, "PLIO_P_ANOM"],
             mio_Tano  = env.new[, "miocene_temp_anom_mean"] / 10,
             mio_Pano  = env.new[, "miocene_prec_anom_mean"]
             ) %>%
  .[.$LEVEL_3_CO %in% tdwg3.included, ]
# NOTE: The new LGM data for FIJ (Fiji) is faulty due to calculation glitches.
env.clim.stable[env.clim.stable[, "LEVEL_3_CO"] == "FIJ",
                c("lgm_Tano", "lgm_Pano")
                ] <- NA
#
# 1 missing value for the LGM anomalies, per the above.
# We should transform bio4_mean to degrees
env.clim.stable$bio4_mean %<>% { . / 100 }

#
# Distributions visually examined with Histogram().
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
  merge(., env.clim, by = "LEVEL_3_CO") %>%
  merge(., env.clim.stable, by = "LEVEL_3_CO") %>%
  merge(.,
        tdwg3.info[, "tdwg3.code", drop = FALSE],
        by.x = "LEVEL_3_CO",
        by.y = "tdwg3.code",
        all.y = TRUE,
        sort = TRUE
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

