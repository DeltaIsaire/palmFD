###################################################
# Palm FD project: Preparation for SAR error models
###################################################
#
# In which data and custom functions are prepared for correlating FD indices with
# environmental predictors, using simultaneous autoregressive error (SAR error)
# models.
# This includes assessing spatial autocorrelaton with Moran's I, calculating and
# choosing neighbourhoods for the spatial weights matrix, and developing model
# selection procedures.
#
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(ncf)
library(sf)
library(spdep)
library(ggplot2)
library(scales)
library(viridis)

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/SAR_regression_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")

# Graphics output directory (with trailing slash!)
plot.dir <- "graphs/SAR_correlograms/"

# Single-predictor model output directory (with trailing slash!)
output.dir <- "output/SAR_models/"


if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
}
if (!dir.exists(output.dir)) {
  cat("creating directory:", output.dir, "\n")
  dir.create(output.dir)
}

set.seed(125)



############################################
# Load functional diversity and spatial data
############################################
cat("loading data...\n")

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
fric <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
          "FRic.fruit.length")
fdis <- c("FDis.all.traits", "FDis.stem.height", "FDis.blade.length",
          "FDis.fruit.length")
for (index in fd.names) {
  fd.indices [[index]] <-
    read.csv(file = paste0("output/FD_summary_", index, ".csv"), row.names = 1)
}


# Environmental data
# ------------------
env.complete <- read.csv(file = "output/tdwg3_predictors_complete.csv",
                         row.names = 1
                         )



############################
# Correlograms and Moran's I
############################
cat("Constructing correlograms...\n")

# First, a wrapper function for ncf::correlog() appropriate for our situation
# ---------------------------------------------------------------------------
MoranPlot <- function(corr, title = NULL, xlab = "Distance class (km)",
                      ylab = "Moran's I") {
# corr: object of class 'correlog', returned by ncf::correlog()
  plot(corr$correlation ~ corr$mean.of.class,
       type = "b",
       main = title,
       pch = 21,
       bg = "black",
       ylab = ylab,
       xlab = xlab,
       ylim = c(-0.5, 1)
       )
  abline(h = 0, lty = 1, col = "black")
}


# Second, a wrapper function to run ncf::correlog()
# -------------------------------------------------
Correlogram <- function(obs, increment = 500, resamp = 999, plot = TRUE, ...) {
# obs: a vector of observed data at each spatial point (TDWG3 unit), for
#         all 368 TDWG3 units.
# title: title for Moran's I plot (e.g. a name for 'obs'). The defaults are
#        approximately appropriate for our data.
# other arguments: see correlog(). Increment is in kilometers in our case.
#
# Returns: creates a plot of Moran's I, and returns the output of correlog()

  df <- data.frame(x = tdwg3.info[, "lon"],
                   y = tdwg3.info[, "lat"],
                   z = obs
                   )
  df %<>% .[complete.cases(.), ]
  corr <- correlog(x = df$x,
                   y = df$y,
                   z = df$z,
                   increment = increment,
                   resamp = resamp,
                   latlon = TRUE
                   )
  if (plot) { MoranPlot(corr = corr, ...) }
  invisible(corr)
}


# Asses spacial autocorrelation (SAC) in FD indices based on all traits
# ---------------------------------------------------------------------
# Quality correlogs take a long time, so execute only if output doesn't exist

cat("(1) For FRic...\n")
# Another small function for automating the FRic plots
MakePlots <- function(x, index) {
  par(mfrow = c(1, 4), mar = c(4.1, 4.1, 4.1, 1.0))
  for (model in c("global.SES", "realm.SES", "realm.SES.noMDG", "adf.SES")) {
    Correlogram(x[, model], plot = TRUE, title = paste0(index, " (", model, ")"))
  }
}

for (i in seq_along(fric)) {
  filename <- paste0(plot.dir, "Morans_I_", fric[i], ".svg")
  if (!file.exists(filename)) {
    GraphSVG(MakePlots(fd.indices[[fric[i]]], index = fric[i]),
             file = filename,
             height = 4,
             width = 16
             )
  }
}

cat("(2) For FDis...\n")
for (i in seq_along(fdis)) {
  filename <- paste0(plot.dir, "Morans_I_", fdis[i], ".svg")
  if (!file.exists(filename)) {
    GraphSVG(Correlogram(fd.indices[[fdis[i]]] [, "observed"],
                         plot = TRUE,
                         title = paste0(fdis[i], " (observed)")
                         ),
             file = filename,
             height = 4,
             width = 5
             )
  }
}



##############################################
# Define spatial weights matrix for SAR models
##############################################
cat("Defining spatial weights matrix...\n")

# Neighbourhood: Sphere of Influence
# ----------------------------------
# It has been determined that the best neighbourhood model to use for our palms
# data in TDWG3 units is the "sphere of influence". It is defined as follows:
#   1. For each polygon, draw a circle around the centroid with radius = distance
#      to nearest neighbour
#   2. Polygons whose "influence" circles overlap are neighbours.
# This is implemented via function SoiNB()
tdwg.map.subset <- tdwg.map[complete.cases(env.complete), ]
nb.soi <- SoiNB(tdwg.map.subset)


# Spatial weights matrix
# ----------------------
# Is created from a neighbourhood object.
# We should use a row standardized weights matrix.
# In order to decrease type II error, it may be worth to weight by distance.
# This is done as follows:
#   1. For each polygon, the (great circle) distances to neighbours:
#   2. The greatest (great circle) distance between any two neighbours:
#   3. For each polygon, weight neighbour distances by the greatest distance
#   4. Create the distance-weighted weights matrix
# spatial weights matrix creation is implemented via function SWMat()
nb.soi.swmat <- SWMat(nb.soi, style = "W")

nb.soi.swmat.dw <- SWMat(nb.soi, tdwg.map.subset, dist.weight = TRUE, style = "W")


# Visualize the default neighbourhood structure
# ---------------------------------------------
# Default meaning the case of all 117 TDWG3 units for which we have data, i.e.
# the 'nb.soi.swmat' neighbourhood defined above.
# The default plot method is plot.nb() from package spdep:
# plot(x = nb.soi.swmat, coords = tdwg.coords)
# Combining it directly with ggplot is difficult, but we can recreate the aesthetic.
# This is implemented with SpatialPlotNB() and Nb2Segments(), which we can call
# with a wrapper function:
NbPlot <- function(nb, tdwg.map, tdwg.map.subset, title = NULL, subtitle = NULL,
                   filename) {
  segments <- Nb2Segments(nb, tdwg.map.subset)
  presence <- tdwg.map$LEVEL_3_CO %in% tdwg.map.subset$LEVEL_3_CO
  presence[!presence] <- NA
  soi.plot <- SpatialPlotNB(tdwg.map[!tdwg.map$LEVEL_3_CO == "ANT", ],
                            presence[!tdwg.map$LEVEL_3_CO == "ANT"],
                            segments,
                            title = title,
                            subtitle = subtitle
                          )
  ggsave(soi.plot, filename = filename, width = 8, height = 4)
  invisible(soi.plot)
}

# Call it for the default SOI neighbourhood:
NbPlot(nb.soi,
       tdwg.map,
       tdwg.map.subset,
       title = "Sphere of Influence Neighbourhood",
       subtitle = "Based on TDWG3 units for which we have complete data",
       filename = paste0(plot.dir, "SOI_neighbourhood_default.png")
       )



#############################
# Single-predictor SAR models
#############################
cat("Generating single-predictor SAR error models:\n")

# - spatial weights matrix with and without distance-weighting
# - all single predictors
# - FRic + FDis
# - null models global, realm, realmnoMDG, adf, and observed
# - global dataset and subset for each three realms
# That is a lot of combinations!
# Fortunately, we can borrow code from the OLS regressions, so developing functions
# for single-predictor SAR models will is relatively easy.


cat("test 1...\n")
model.data <- env.complete
model.data[, "null"] <- runif(n = nrow(model.data))
model.data[, "response"] <- fd.indices[[1]] [, "global.SES"]

result <- SingleSAR(model.data,
                    listw = nb.soi.swmat,
                    response = "response"
                    )


cat("test 2...\n")
responses <- fd.indices[[1]] [, c("global.SES", "realm.SES", "adf.SES")]
# TODO: if you use 'realm.SES.noMDG' you need a special neighbourhood, where MDG
# is removed!

result2 <- MultiSingleSAR(responses = responses,
                          predictors = env.complete,
                          listw = nb.soi.swmat
                          )

cat("test 3...\n")
result3 <- RunMSSAR(name = paste0(output.dir, "SAR_single_test"),
                    responses = responses,
                    predictors = env.complete,
                    listw = nb.soi.swmat
                    )

cat("test 4...\n")

# Function to run all models via RunMSSAR
# ---------------------------------------
AllSARSingles <- function(fd.indices, colname, predictors, tdwg.map,
                          dist.weight = FALSE, name.all, ...) {
# Args:
#   fd.indices: the list object 'fd.indices' or a subset thereof
#   colname: name of column to extract from each df in 'fd.indices' via GetFD()
#   predictors: dataframe with predictor variables
#   tdwg.map: the map object as created by read_sf(). Is automatically
#             subsetted as required. Used for creating the spatial weights matrix.
#   dist.weight: whether to apply neighbour distance weighting in the spatial
#                weights matrix. See SWMat()
#   name.all: 'name' argument passed to RunMSSAR(). Should be of length 1; suffixes
#         will be added for the cases 'full' and each realm subset
#   ...: Additional arguments to pass to RunMSSAR()

  fd.all <- GetFD(fd.indices, colname)

  # full dataset:
  cat(name.all, "full...\n")
  ind.complete <- complete.cases(fd.all) & complete.cases(predictors)
  tdwg.map.subset <- tdwg.map[ind.complete, ]
  nb <- SoiNB(tdwg.map.subset)
  swmat <- SWMat(nb, tdwg.map.subset, dist.weight = dist.weight, style = "W")
  RunMSSAR(name = paste0(name.all, "_full"),
           responses = fd.all,
           predictors = predictors,
           listw = swmat,
           ...
           )

  # realm subsets:
  fd.all.realms <- RealmSubset(fd.all)
  predictors.realms <- RealmSubset(predictors[, !colnames(predictors) %in% "realm"])
  for (i in seq_along(fd.all.realms)) {
    cat(name.all, names(fd.all.realms)[i], "...\n")
    ind.complete <-
      complete.cases(fd.all.realms[[i]]) & complete.cases(predictors.realms[[i]])
    tdwg.ind <- tdwg.map$LEVEL_3_CO %in% rownames(fd.all.realms[[i]])[ind.complete]
    tdwg.map.subset <- tdwg.map[tdwg.ind, ]
    nb <- SoiNB(tdwg.map.subset)
    swmat.realm <- SWMat(nb, tdwg.map.subset, dist.weight = dist.weight, style = "W")
    RunMSSAR(name = paste0(name.all, "_", names(fd.all.realms)[i]),
             responses = fd.all.realms[[i]],
             predictors = predictors.realms[[i]],
             listw = swmat.realm,
             ...
             )
  }
  return (0)
}

# call the function
# -----------------
AllSARSingles(fd.indices[fric],
              colname = "global.SES",
              predictors = env.complete,
              tdwg.map = tdwg.map,
              dist.weight = FALSE,
              name.all = paste0(output.dir, "SAR_single_test_FRic_global"),
              standardize = TRUE,
              numeric.only = FALSE
              )






cat("Done.\n")

