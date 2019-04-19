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
tdwg.map.subset <- tdwg.map[complete.cases(env.complete), ]
tdwg.spatial <- as(tdwg.map.subset, "Spatial")
tdwg.coords <- coordinates(tdwg.spatial)
nb.dlny <- tri2nb(tdwg.coords)
nb.soi <-
  soi.graph(nb.dlny, tdwg.coords) %>%
  graph2nb()


# Spatial weights matrix
# ----------------------
# Using row standardization (style = "W" in nb2listw())
nb.soi.swmat <- nb2listw(nb.soi, style = "W")
# Neat statistics available via summary(nb.soi.wtm)

# In order to decrease type II error, it may be worth to weight by distance.
# Thus we can create an alternative spatial weights matrix:
#   1. For each polygon, the (great circle) distances to neighbours:
nb.soi.dist <- nbdists(nb.soi, tdwg.coords, longlat = TRUE)
#   2. The greatest (great circle) distance between any two neighbours:
nb.soi.maxdist <- max(unlist(nb.soi.dist))
#   3. For each polygon, weight neighbour distances by the greatest distance
nb.soi.wdist <- lapply(nb.soi.dist, function(x) { 1 - x / nb.soi.maxdist } )
#   4. Create distance-weighted weights matrix:
nb.soi.swmat.dw <- nb2listw(nb.soi, glist = nb.soi.wdist, style = "W")


# Visualize the neighbourhood structure
# -------------------------------------
# The default plot method is plot.nb() from package spdep:
# plot(x = nb.soi.swmat, coords = tdwg.coords)
# Combining it directly with ggplot is difficult, but we can recreate the aesthetic.
#
# Function to plot the graph
SpatialPlotNB <- function(tdwg.map, presence, segments, title = NULL,
                          subtitle = NULL) {
# Segments: data.frame with coordinates for segments. See geom_segment().
#           The order of columns in the dataframe should be
#           c("x", "y", "xend", "yend")
  tdwg.plot <- ggplot(data = tdwg.map) +
                 geom_sf(size = 0.25,
                         color = "white",
                         aes(fill = presence),
                         show.legend = FALSE
                         ) +
                 # This magically only adds axes:
                 geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
                 labs(x = NULL,
                      y = NULL,
                      title = title,
                      subtitle = subtitle
                      ) +
                 geom_point(data = tdwg.map[which(presence), ],
                            aes(x = Long, y = Lat),
                            pch = 21,
                            size = 2
                            ) +
                 geom_segment(data = segments,
                              aes(x = segments[, 1],
                                  y = segments[, 2],
                                  xend = segments[, 3],
                                  yend = segments[, 4]
                                  ),
                              size = 0.35,
                              show.legend = FALSE
                              ) +
                 scale_fill_viridis(na.value = "#C0C0C0",
                                    discrete = is.discrete(presence),
                                    option = "inferno",
                                    begin = 0.75
                                    )
  tdwg.plot
}

# Function to create segments dataframe from a neighbourhood object
Nb2Segments <- function(nb, tdwg.map, presence) {
  nb.num <- laply(nb, length)
  present.names <- tdwg.map$LEVEL_3_CO[which(presence)]
  segments <-
    matrix(data = NA,
           nrow = sum(nb.num),
           ncol = 4,
           dimnames = list(rep(present.names, times = nb.num),
                           c("x", "y", "xend", "yend")
                           )
                     )
  segments[, 1] <- tdwg.map$Long[match(rownames(segments), tdwg.map$LEVEL_3_CO)]
  segments[, 2] <- tdwg.map$Lat[match(rownames(segments), tdwg.map$LEVEL_3_CO)]
  indices <- match(present.names[unlist(nb)], tdwg.map$LEVEL_3_CO)
  segments[, 3] <- tdwg.map$Long[indices]
  segments[, 4] <- tdwg.map$Lat[indices]
  as.data.frame(segments)
}

# Make the plot
presence <- complete.cases(env.complete)
presence[!presence] <- NA
segments <- Nb2Segments(nb.soi, tdwg.map, presence)
soi.plot <- SpatialPlotNB(tdwg.map[!tdwg.map$LEVEL_3_CO == "ANT", ],
                          presence[!tdwg.map$LEVEL_3_CO == "ANT"],
                          segments,
                          title = "Sphere of Influence Neighbourhood",
                          subtitle = "Based on TDWG3 units for which we have complete data"
                          )
ggsave(soi.plot,
       filename = paste0(plot.dir, "SOI_neighbourhood.png"),
       width = 8,
       height = 4
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




cat("Done.\n")

