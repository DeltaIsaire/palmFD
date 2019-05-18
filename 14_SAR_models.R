###################################################
# Palm FD project: Preparation for SAR error models
###################################################
#
# In which data and custom functions are prepared for correlating FD indices with
# environmental predictors, using simultaneous autoregressive error (SAR error)
# models.
# This includes assessing spatial autocorrelaton with Moran's I, calculating and
# choosing neighbourhoods for the spatial weights matrix, and implementing
# single-predictor models.
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
library(spatialreg)

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
env.complete.noch <- read.csv(file = "output/tdwg3_predictors_complete_noch.csv",
                              row.names = 1
                              )
# The 'noch' (no canopy height) version has more data points, because canopy
# height has some missing values. Since canopy height is a questionable predictor,
# We don't want to unnecessarily throw away data.



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

null.models <- c("global.SES", "realm.SES", "realm.SES.noMDG", "adf.SES", "observed")

# Another small function for automating the FRic plots
MakePlots <- function(x, index) {
  par(mfrow = c(5, 1), mar = c(4.1, 4.1, 4.1, 1.0))
  for (model in null.models) {
    Correlogram(x[, model], plot = TRUE, title = paste0(index, " (", model, ")"))
  }
}

cat("(1) For FRic...\n")
for (i in seq_along(fric)) {
  filename <- paste0(plot.dir, "Morans_I_", fric[i], ".svg")
  if (!file.exists(filename)) {
    GraphSVG(MakePlots(fd.indices[[fric[i]]], index = fric[i]),
             file = filename,
             height = 20,
             width = 5
             )
  }
}

cat("(2) For FDis...\n")
for (i in seq_along(fdis)) {
  filename <- paste0(plot.dir, "Morans_I_", fdis[i], ".svg")
  if (!file.exists(filename)) {
    GraphSVG(MakePlots(fd.indices[[fdis[i]]], index = fdis[i]),
             file = filename,
             height = 20,
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
tdwg.map.subset <- tdwg.map[complete.cases(env.complete.noch), ]
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
# Note that this function applies a correct for Samoa, because SOI neighbourhood
# is graph-based, which doesn't wrap around the globe (as it doesn't use great
# circle distance)
nb.soi.swmat <- SWMat(nb.soi, style = "W")

nb.soi.swmat.dw <- SWMat(nb.soi, tdwg.map.subset, dist.weight = TRUE, style = "W")


# Visualize the default neighbourhood structure
# ---------------------------------------------
# Default meaning the case of all 126 TDWG3 units for which we have data, i.e.
# the 'nb.soi.swmat' neighbourhood defined above.
# The default plot method is plot.nb() from package spdep:
# plot(x = nb.soi.swmat, coords = tdwg.coords)
# Combining it directly with ggplot is difficult, but we can recreate the aesthetic.
# This is implemented with SpatialPlotNB() and Nb2Segments(), which we can call
# with the wrapper function NbPlot()


# Call it for the default SOI neighbourhood:
NbPlot(nb.soi,
       tdwg.map,
       tdwg.map.subset,
       title = "Sphere of Influence Neighbourhood",
       subtitle = "Based on TDWG3 units with complete data",
       filename = paste0(plot.dir, "SOI_neighbourhood_default.png")
       )



#############################
# Single-predictor SAR models
#############################
cat("Generating single-predictor SAR error models:\n")

# Wrapper function for awesomeness
# --------------------------------
RunSARSingles <- function(fd.indices, colname, name.all, dist.weight = FALSE,
                          double.std = FALSE) {
# Generate the data, then immediately parse it
  AllSARSingles(fd.indices,
                colname = colname,
                predictors = env.complete.noch,
                tdwg.map = tdwg.map,
                dist.weight = dist.weight,
                name.all = name.all,
                standardize = TRUE,
                numeric.only = FALSE,
                double.std = double.std
                )
  ParseSARSingle(name.all = name.all,
                 cases = c("full", "NewWorld", "OWWest", "OWEast"),
                 statistics = c("slope", "p.value", "moran.p", "pseudo.Rsq"),
                 filter.p = TRUE
                 )
}


# Run the single-predictor models
# -------------------------------
# For FRic:
null.models
SAR.FRic <- vector("list", length = length(null.models))
names(SAR.FRic) <- null.models
for (i in seq_along(SAR.FRic)) {
  cat("For FRic", null.models[i], "\n")
  SAR.FRic[[i]] <-
    RunSARSingles(fd.indices[fric],
                  colname = null.models[i],
                  name.all = paste0(output.dir, "SAR_single_FRic_", null.models[i]),
                  dist.weight = FALSE,
                  double.std = TRUE
                  )
}
# For FRic With distance weighted spatial weights matrix:
SAR.FRic.dw <- SAR.FRic
for (i in seq_along(SAR.FRic.dw)) {
  cat("For FRic", null.models[i], "with distance-weighted swmat\n")
  SAR.FRic.dw[[i]] <-
    RunSARSingles(fd.indices[fric],
                  colname = null.models[i],
                  name.all = paste0(output.dir,
                                    "SAR_single_FRic_dw_",
                                    null.models[i]
                                    ),
                  dist.weight = TRUE,
                  double.std = TRUE
                  )
}

# For FDis:
null.models
SAR.FDis <- vector("list", length = length(null.models))
names(SAR.FDis) <- null.models
for (i in seq_along(SAR.FDis)) {
  cat("For FDis", null.models[i], "\n")
  SAR.FDis[[i]] <-
    RunSARSingles(fd.indices[fdis],
                  colname = null.models[i],
                  name.all = paste0(output.dir, "SAR_single_FDis_", null.models[i]),
                  dist.weight = FALSE,
                  double.std = TRUE
                  )
}
# For FDis with distance weighted spatial weights matrix:
SAR.FDis.dw <- SAR.FDis
for (i in seq_along(SAR.FDis.dw)) {
  cat("For FDis", null.models[i], "with distance-weighted swmat\n")
  SAR.FDis.dw[[i]] <-
    RunSARSingles(fd.indices[fdis],
                  colname = null.models[i],
                  name.all = paste0(output.dir,
                                    "SAR_single_FDis_dw_",
                                    null.models[i]
                                    ),
                  dist.weight = TRUE,
                  double.std = TRUE
                  )
}


# Parse and save results
# ----------------------
cat("Parsing and saving results...\n")

DoParse <- function(output, filename) {
  for (i in seq_along(output[[1]])) {
    write.csv(ldply(output, function(x) { x[[i]] } , .id = "null.model" ),
              file = paste0(output.dir,
                            filename,
                            "_",
                            names(output[[1]])[i],
                            ".csv"
                            ),
              eol = "\r\n",
              row.names = FALSE,
              quote = FALSE
              )
  }
}

DoParse(SAR.FRic, filename = "00_SAR_single_output_FRic")
DoParse(SAR.FRic.dw, filename = "00_SAR_single_output_FRic_dw")
DoParse(SAR.FDis, filename = "00_SAR_single_output_FDis")
DoParse(SAR.FDis.dw, filename = "00_SAR_single_output_FDis_dw")



cat("Done.\n")

