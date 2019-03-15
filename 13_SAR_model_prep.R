###################################################
# Palm FD project: Preparation for SAR error models
###################################################
#
# In which data and custom functions are prepared for correlating FD indices with
# environmental predictors, using simultaneous autoregressive error (SAR error) models.
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

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


############################################
# Load functional diversity and spatial data
############################################
cat("loading data...\n")

# For now, we will focus on the all-traits data
fric <- read.csv(file = "output/FD_summary_FRic.csv", header = TRUE, row.names = 1)
fdis <- read.csv(file = "output/FD_summary_FDis.csv", header = TRUE, row.names = 1)

tdwg3.info <- read.csv(file = "output/tdwg3_info.csv")


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

Correlogram <- function(obs, increment = 500, resamp = 999, title = NULL,
                        plot = TRUE) {
# obs: a vector of observed data at each spatial point (TDWG3 unit), for
#         all 368 TDWG3 units.
# title: title for Moran's I plot (e.g. a name for 'obs'). The defaults are
#        approximately appropriate for our data.
# other arguments: see correlog()
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
  if (plot) { MoranPlot(corr, title) }
  invisible(corr)
}

# Asses spacial autocorrelation (SAC) in SES FD indices
# -----------------------------------------------------
sac.fric.global <- Correlogram(fric$FRic.global.SES, plot = FALSE)
sac.fric.realm <- Correlogram(fric$FRic.realm.SES, plot = FALSE)
sac.fric.adf <- Correlogram(fric$FRic.adf.SES, plot = FALSE)

sac.fdis.global <- Correlogram(fdis$FDis.global.SES, plot = FALSE)
sac.fdis.realm <- Correlogram(fdis$FDis.realm.SES, plot = FALSE)
sac.fdis.adf <- Correlogram(fdis$FDis.adf.SES, plot = FALSE)

PlotAll <- function() {
  par(mfrow = c(2, 3), mar = c(4.1, 4.1, 4.1, 1.0))
  MoranPlot(sac.fric.global, title = "FRic (global)")
  MoranPlot(sac.fric.realm, title = "FRic (realm)")
  MoranPlot(sac.fric.adf, title = "FRic (adf)")
  MoranPlot(sac.fdis.global, title = "FDis (global)")
  MoranPlot(sac.fdis.realm, title = "FDis (realm)")
  MoranPlot(sac.fdis.adf, title = "FDis (adf)")
}

GraphSVG(PlotAll(),
         file = "graphs/FD_all_traits_Morans_I.svg",
         height = 6,
         width = 12
         )




cat("Done.\n")

