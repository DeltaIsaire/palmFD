
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


