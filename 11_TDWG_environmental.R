###########################################
# Palm FD project: TDWG3 environmental data
###########################################
#
# In which environmental data at the resolution of botanical countries (TDWG3
# units) is prepared for regression analysis.
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

# Directory for saving plots (with a trailing slash):
plot.dir = "graphs/environmental/"


###################################
# Load and parse environmental data
###################################
cat("Loading and parsing environmental data...\n")

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


# LGM climate (anomaly)
# ---------------------
env.lgm <-
  env.new[, c("LEVEL_3_CO", "lgm_ens_Tmean", "lgm_ens_Pmean")] %>%
  .[.$LEVEL_3_CO %in% tdwg3.included, ]
# The old data was based only on CCSM and MIROC, while the new data
# includes more sources. So new should be 'better', and some deviations are
# to be expected, which is OK.
# The new data is in K * 10 while the old data is in K.
# They do however have the same SD, which is mildly weird.
#
# NOTE: The new LGM data for FIJ (Fiji) is faulty due to calculation glitches.
env.lgm[match("FIJ", env.lgm$LEVEL_3_CO), -1] <- NA
#
# LGM mean temperature is right-skewed, with 2 extreme values on the low end.
#   These extremes are CHT (Tibet) and WHM (West himalaya), which makes sense.
# LGM Mean precipitation is slightly bimodal but otherwise fine.


#######################################
# Exploratory PCA on environmental data
#######################################
# Just to get some idea of collinearities.
# We will use PCA based on a correlation matrix, i.e. with standardization of 
# variables. This is suitable for environmental variables.

# PCA on non-climatic environmental variables
# -------------------------------------------
pca.nonclim <- rda(env.nonclim[complete.cases(env.nonclim), -1], scale = TRUE)
summary(pca.nonclim)
OrdinationPlot(pca.nonclim)
# Looks like mean altitude and altitude range are strongly related.
# Canopy Height mean and range are also related.
# For both altitude and canopy height, pick either range or mean. For our study,
# the range might make the most sense.

# PCA on contemporary climate variables
# -------------------------------------
pca.clim <- rda(env.clim[complete.cases(env.clim), -1], scale = TRUE)
summary(pca.clim)
# From this, it looks like only 4-6 axes define most (88-97%) variation in
# contemporary climate. So for future regressions, we may have to focus on
# a rather limited subset of the 19 bioclim variables.
OrdinationPlot(pca.clim)
# This plot suggests a more nuanced view. Additional collinearity analysis will
# be needed.
# The climatic variables cluster roughly around two axes, related to
# (1) temperature (mean and extremes)
# (2) everything else: precipitation, and temperature seasonality.

# PCA on LGM climate
# ------------------
# These are only two variables, but still interesting to see a PCA of it.
pca.lgm <- rda(env.lgm[complete.cases(env.lgm), -1], scale = TRUE)
summary(pca.lgm)
OrdinationPlot(pca.lgm)
# That does not look like a correlation at all. Which is good.

# PCA on combined environmental data
# ----------------------------------
env.all <- merge(env.nonclim, merge(env.clim, env.lgm))
pca.all <- rda(env.all[complete.cases(env.all), -1], scale = TRUE)
summary(pca.all)
# >95% of variation is explained with 8 variables, which is a rough initial
# indication of how many predictors could be important, at most.
OrdinationPlot(pca.all)
# The only new insight here is that the LGM climate is likely strongly correlated
# with contemporary climate, especially precipitation.


###################################################
# Assessing collinearity in environmental variables
###################################################
# Correlation matrix
# ------------------
corr.matrix <- 
  cor(env.all[complete.cases(env.all), -1]) %>%
  round(., digits = 2)
# Predictably, this is a mess. Going through it is entirely possible, but it might
# be easier to create a correlogram.
corrplot(corr.matrix, method = "square", order = "FPC")
# Pay attention to the dark blue and dark red.
# Clear collinearities are the following:
# Mean altitude is strongly negatively related to temperature, LGM temperature
#   and quarter temperatures (bio6, bio 8-11).
# Altitude range shows similar correlations, but less strong.
# Mean annual temperature is strongly related to mean altitude, LGM temperature
#   and quarter temperatures (bio6, bio 8-11)
# mean diurnal range (bio2) is related to annual temperature range (bio7),
#   annual precipitation and LGM precipitation.
# Temperature seasonality is related to bio6-7 and bio11
# Temperature of coldest month is related to temperature of warmest quarter.
# annual temperature range is strongly related to temperature of coldest month
# bio8-11 are moderately correlated with each other.
# bio 12-19 are weakly to moderately correlated with each other, especially
#  bio15 (precipitation seasonality) with bio14 and bio17.
# LGM temperature is strongly related to mean altitude, mean annual temperature
#   and bio8-11.
# LGM precipitation is strongly related to annual precipitation, and to
#  bio13-19 (all precipitation-related).
#
# In summary, keeping in mind the PCA results, we can tentatively identify
# The following variables as containing the main variation:
# - Mean annual temperature (related to lgm Temperature, bio6, 8-11, mean
#   altitude, and to a lesser extent temperature range and temperature seasonality)
# - annual precipitation (related to LGM precipitation and other prec variables)
# - Precipitation seasonality
# - Temperature seasonality (related to bio2, bio 6-7)
# - canopy height mean
# - canopy height range
# - soil count
# - bio3 (isothermality)
#
# So altitude is not in there (makes sense),
# LGM climate is not in there (sort of understandable),
# And the variation in temperature and precipitation is captured pretty well
# with the primary variables mean temp/precip and seasonality of temp/precip.




# In addition to simple collinearity, there might be some multicollinearity, which
# is hard to find from a correlogram or a correlation matrix. Using VIF will be
# a good complement.


# Variance Inflation Factor (VIF)
# -------------------------------
# TODO: go selection based on VIF, with a threshold of 2-3 (zuur et al 2010), or
# 5 if you're feeling generous.
# NOTE: VIF works on regression model output, so it will have to wait until we
#       start running regressions.



cat("Done.\n")

