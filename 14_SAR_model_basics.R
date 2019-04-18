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

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")

# Graphics output directory (with trailing slash!)
plot.dir <- "graphs/SAR_correlograms/"


if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
}


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
MakePlots <- function(x, index) {
  par(mfrow = c(1, 4), mar = c(4.1, 4.1, 4.1, 1.0))
  for (model in c("global.SES", "realm.SES", "realm.SES.noMDG", "adf.SES")) {
    Correlogram(x[, model], plot = TRUE, title = paste0(index, " (", model, ")"))
  }
}

# For FRic:
for (i in seq_along(fric)) {
  GraphSVG(MakePlots(fd.indices[[fric[i]]], index = fric[i]),
           file = paste0(plot.dir, "Morans_I_", fric[i], ".svg"),
           height = 4,
           width = 16
           )
}

# For FDis:
for (i in seq_along(fdis)) {
  GraphSVG(Correlogram(fd.indices[[fdis[i]]] [, "observed"],
                       plot = TRUE,
                       title = paste0(fdis[i], " (observed)")
                       ),
           file = paste0(plot.dir, "Morans_I_", fdis[i], ".svg"),
           height = 4,
           width = 5
           )
}




cat("Done.\n")

