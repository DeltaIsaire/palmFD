###########################################
# Palm FD project: TDWG3 environmental data
###########################################
#
# In which environmental data at the resolution of botanical countries (TDWG3
# units) is prepared for regression analysis. This includes assessment of extreme
# values, testing for approximate normality and transforming if required, and
# standardizing the (transformed) predictors to mean 0 and unit variance.
# Additionally, PCA on predictor variables is performed, and PCA coordinates
# of contemporary climate data are saved for use as alternative predictors.
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
source(file = "functions/OLS_regression_functions.R")


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


###########################################
# Exploratory PCA and Collinearity analysis
###########################################
cat("Performing exploratory PCA and producing corrplots...\n")

# We will use PCA based on a correlation matrix, i.e. with standardization of 
# variables. This is suitable for environmental variables.
# Corrplots will help identify how predictors are correlated with each other,
# and which predictors (if any) are redundant. Additionally, VIF can be used
# to identify the most collinear predictors.

# All data
# --------
cat("(1) for all predictors together...\n")
pca.env <- rda(predictors[complete.cases(predictors), ], scale = TRUE)
if (FALSE) {
  summary(pca.env)
  OrdinationPlot(pca.env)
}
GraphSVG(Correlogram(predictors),
         file = "graphs/environmental_predictors_corrplot.svg",
         width = 16,
         height = 16
         )
x <- predictors[!names(predictors) %in% c("bio5_mean", "bio6_mean",
                                          "lgm_ens_Tmean", "lgm_ens_Pmean"
                                          )
                ]
x$response <- rnorm(n = nrow(x), mean = 10, sd = 1)
AutoVIF(x, response = "response", threshold = 3)

# There are 26 predictors in total.
# >95% of variation is explained with 9 axes, which is a rough indication of how
# many predictors could be important, at most, in a full model.
# There is no clear clustering in the plot. Arguably, these are too many predictors
# to draw any clear conclusions. Best look at some subsets.
# Based on VIF, non-collinear predictors are limited to a set of 10:
# alt_range, soilcount, CH_mean, CH_Range, bio2, bio3, bio8, bio13, bio19,
# and lgm precipitation anomaly.


# Nonclimatic environmental data
# ------------------------------
cat("(2) for Nonclimatic environmental predictors...\n")
pca.nonclim <- rda(env.nonclim[complete.cases(env.nonclim), ] [, -1], scale = TRUE)
if (FALSE) {
  summary(pca.nonclim)
  OrdinationPlot(pca.nonclim)
}
GraphSVG(Correlogram(env.nonclim[, -1]),
         file = "graphs/environmental_predictors_corrplot_nonclim.svg",
         width = 6,
         height = 6
         )
x <- env.nonclim
x$response <- rnorm(n = nrow(x), mean = 10, sd = 1)
AutoVIF(x, response = "response", threshold = 3)
# There are 5 nonclimatic environmental predictors.
# The cumulative proportion explained is 0.5 / 0.76 / 0.9 / 0.97 / 1,
# indicating that 2-3 predictors will capture most variation.
# Mean altitude and altitude range are strongly correlated and shouldn't be used
# in the same model. Pick one: for our purposes, range is the best choice because
# it measures variation; although VIF is highest for alt_range.


# Contemporary Climate data
# -------------------------
cat("(3) for contemporary climate predictors...\n")
pca.clim <- rda(env.clim[complete.cases(env.clim), ] [, -1], scale = TRUE)
if (FALSE) {
  summary(pca.clim)
  OrdinationPlot(pca.clim)
}
GraphSVG(Correlogram(env.clim[, -1]),
         file = "graphs/environmental_predictors_corrplot_clim.svg",
         width = 12,
         height = 12
         )
x <- env.clim[, !names(env.clim) %in% c("bio5_mean", "bio6_mean")]
x$response <- rnorm(n = nrow(x), mean = 10, sd = 1)
AutoVIF(x, response = "response", threshold = 3)
# There are 18 included climate predictors.
# >95% of variation is explained with 6 PCA axes, so there is definitely a lot
# of redundancy. The first two axes explain 73%. 
# An 18x18 correlation matrix is quite large, but there are a number of strong
# correlations.
# Based on VIF, the unique variables are bio3, bio9, bio10, bio13, bio15, bio19,
# which is 6 predictors. A curiously coincidental number.
# 3/10/13/15 have the lowest VIF scores, so are most unique.
# Those four are Isothermality (T diurnal range / annual T range), mean T of warmest
# quarter, Precipitation of wettest month, and precipitation seasonality.
# The additional variables selected by VIF are mean T of driest quarter, and
# Precipitation of coldest quarter.
#
# Save the first 4 axes as alternative predictors. Together these explain 89.3%
# of the variation in climate.
# Remember these PCA scores are based on all bioclim variables except bio5 and bio6.
env.clim.pca <- as.data.frame(summary(pca.clim) [["sites"]] [, 1:4])
# Fix rownames and pad to all TDWG3 units:
env.clim.pca[, "tdwg3.code"] <- env.new[rownames(env.clim.pca), "LEVEL_3_CO"]
env.clim.pca %<>%
  merge(.,
        tdwg3.info[, "tdwg3.code", drop = FALSE],
        by = "tdwg3.code",
        all.y = TRUE,
        sort = TRUE
        )
rownames(env.clim.pca) <- env.clim.pca[, "tdwg3.code"]
env.clim.pca %<>% .[, -1]
if (!identical(rownames(env.clim.pca), rownames(predictors))) {
  stop("climatic PCA scores aren't parsed properly")
}
# Save final result:
write.csv(env.clim.pca,
          file = "output/tdwg3_environmental_predictors_climate_pca_axes.csv",
          eol = "\r\n",
          row.names = TRUE
          )


# LGM climate data
# ----------------
cat("(4) for LGM climate predictors...\n")
pca.lgm <- rda(env.lgm[complete.cases(env.lgm), ] [, -1], scale = TRUE)
if (FALSE) {
  summary(pca.lgm)
  OrdinationPlot(pca.lgm)
}
GraphSVG(Correlogram(env.lgm[, -1]),
         file = "graphs/environmental_predictors_corrplot_lgm.svg",
         width = 6,
         height = 5
         )
x <- env.lgm
x$response <- rnorm(n = nrow(x), mean = 10, sd = 1)
AutoVIF(x, response = "response", threshold = 3)
# There are four predictors here.
# 3 PCA axes explain 99% of variation. The ordination plot and corrplot reveal
# this is due to a very strong relation between LGM T and LGM anomT.
# Based on VIF, LGM T is the predictor to exclude.
# There is very little collinearity between the other 3 variables. 


# Synthesis: predictor choice
# ---------------------------
# Nonclim:
#   Altitude is often related to climate, but altitude range could capture environ-
#   mental heterogeneity not covered by climatic variables. So Mean altitude can be
#   discarded in favor of alt_range. 
#   A similar argument can NOT be made for canopy height, because CH_Mean is not
#   strongly correlated with CH_Range (r 0.47). Soil count is similarly non-
#   collinear.
#  >>> choose alt_range, soilcount, CH_Range, CH_Mean
#
# LGM climate:
#   LGM temperature is too strongly related to Temperature anomaly.
#   Both T and P anomalies are linearily related to contemporary + LGM T and P.
#   This means we can only have 2 out of 3 of the following: contemporary climate,
#   LGM climate, and LGM climate anomaly. Contemporary climate should not be
#   eliminated, thus forcing a choice between LGM climate or LGM anomalies.
#   Filtering effects due to environmental change are best captured by the
#   anomalies, so let's go with those.
#  >>> choose LGM_Tano, LGM_Pano
#
# Contemporary climate:
#   Theoretically, measures of change and variability, such as T and P seasonality,
#   should be most strongly correlated with functional diversity. That still leaves
#   lots of options. Here, using PCA axes could be a good choice, particularly
#   considering that the first 2 axes already explain 73% of variation.
#   The alternative is careful selection.
#   Theoretically, the most interesting variables are (bio1), bio3, bio4, (bio12),
#   and bio15.
x <- env.clim[, names(env.clim) %in% c("bio1_mean", "bio3_mean", "bio4_mean",
                                       "bio12_mean", "bio15_mean"
                                       )
              ]
x$response <- rnorm(n = nrow(x), mean = 10, sd = 1)
AutoVIF(x, response = "response", threshold = 3)
#   VIF says bio4 (T seasonality) is collinear with the other four predictors.
#
# This script has saved the first four axes of PCA on all climate data (except bio5
# and bio6). These axes are by definition non-collinear, and will capture the
# complete influence of climate. The cumulative proportion explained is
# 43.9% / 72.9% / 81.7% / 89.3%.


cat("Done.\n")

