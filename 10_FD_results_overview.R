##########################################################
# Palm FD project: Summary and visualization of FD results
##########################################################
#
# In which the calculated FD indices (both observed and null model z-scores)
# are summarized and visualized.
#
# Input files:
#   output/observed_FD/FD_observed_genus_mean.csv
#   output/observed_FD/FD_observed_stochastic_mean_of.csv
#   output/null_model_global/FD_z_scores_global_genus_mean_z.scores.csv
#   output/null_model_global/FD_z_scores_global_stochastic_mean_of_z.scores.csv
#   output/null_model_regional/FD_z_scores_regional_genus_mean_z.scores.csv
#   output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_z.scores.csv
#   output/tdwg3_info.csv
#   output/palm_tdwg3_pres_abs_gapfilled.csv
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(sf)
library(ggplot2)
library(scales)

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")

# Directory for saving plots (with a trailing slash):
plot.dir = "graphs/FD_distributions/"


###################
# Load spatial data
###################

# Info on tdwg3 units
# -------------------
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

# Spatial map of tdwg3 units
# --------------------------
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")
# Note I am using the 'sf' spatial format instead of 'SpatialPolygonsDataFrame'.
# The sf package with function sf::read_sf is more modern than using rgdal::readOGR.
# That modernity translates into improved performance.


####################
# Load FD index data
####################
cat("Preparing data...\n")

# After loading, immediately pad the data to include all tdwg3 units.
# This is helpful for plotting the data.
Temp <- function(x) {
  x[, "tdwg3.code"] <- rownames(x)
  merge(x = tdwg3.info[, "tdwg3.code", drop = FALSE],
        y = x,
        by = "tdwg3.code",
        all.x = TRUE,
        sort = TRUE
        )
}

# Observed FD
# -----------
fd.observed.mean <-
  read.csv(file = "output/observed_FD/FD_observed_genus_mean.csv",
           row.names = 1
           ) %>%
  Temp()
fd.observed.stochastic <-
  read.csv(file = "output/observed_FD/FD_observed_stochastic_mean_of.csv",
           row.names = 1
           ) %>%
  Temp()

# Global null model z-scores
# --------------------------
fd.global.mean <-
  read.csv(file = "output/null_model_global/FD_z_scores_global_genus_mean_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()
fd.global.stochastic <-
  read.csv(file = "output/null_model_global/FD_z_scores_global_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()

# Regional null model z-scores
# ----------------------------
fd.regional.mean <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_genus_mean_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()
fd.regional.stochastic <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()

# Local null model z-scores
# -------------------------
fd.local.mean <-
  read.csv(file = "output/null_model_local/FD_z_scores_local_genus_mean_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()
fd.local.stochastic <-
  read.csv(file = "output/null_model_local/FD_z_scores_local_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()

# Community mean trait values
# ---------------------------
# Keep in mind these are log10-transformed
cwm.observed.mean <- read.csv(file = "output/test/observed_FD/test_community_trait_means_genus_mean.csv",
                              row.names = 1
                              ) %>%
  Temp()
cwm.observed.stochastic <- read.csv(file = "output/test/observed_FD/test_community_trait_means_genus_mean.csv",
                                    row.names = 1
                                    ) %>%
  Temp()


##############################################
# Generalized tdwg3 spatial plotting functions
##############################################
# CODING NOTE:
# With sf and ggplot2, instead of using geom_polygon(), we use geom_sf() with the
# 'sf' object. And instead of ggplot2::fortify(), we pretend the 'sf' object
# is a dataframe and simply add our data to it as a column.
# See the function SpatialPlot() below for an example.
# The geom_sf() function is included in the github version of ggplot2, but not
# in the CRAN version. 
# For info and tutorial see https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
#

# Main plotting function
# ----------------------
SpatialPlot <- function(tdwg.map, vector, vector.name, title = NULL,
                        subtitle = NULL) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units).
#         Data must be a continuous numeric variable.
# vector.name: character string giving name for the vector (used in legend)
  tdwg.map$vector <- vector
  subset <- tdwg.map[!is.na(tdwg.map$vector), ]
  min <- ifelse(min(vector, na.rm = TRUE) < 0,
                min(vector, na.rm = TRUE),
                0 - 2 * .Machine$double.eps
                )
  max <- ifelse(max(vector, na.rm = TRUE) > 0,
                max(vector, na.rm = TRUE),
                0 + 2 * .Machine$double.eps
                )
  ggplot(data = tdwg.map) + 
         geom_sf(size = 0.15, color = "black") +
         # This magically only adds axes:
         geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
         # While this does NOT add axes but does add points:
         geom_point(data = subset,
                    mapping = aes(x = Long, y = Lat, fill = vector),
                    size = 3,
                    shape = 21) +
         labs(fill = vector.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              ) +
         scale_fill_gradientn(colours = c("blue", "cyan", "white", "yellow", "red"),
                              values = rescale(c(min, 
                                                 0 - .Machine$double.eps,
                                                 0,
                                                 0 + .Machine$double.eps,
                                                 max
                                                 )
                                               )
                              )
}

# Alternative plotting function using filled polygons instead of points
# ---------------------------------------------------------------------
SpatialPlotFill <- function(tdwg.map, vector, vector.name, title = NULL) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units)
# vector.name: character string giving name for the vector (used in legend)
  tdwg.map$vector <- vector
  tdwg.plot <- ggplot(data = tdwg.map) + 
                 geom_sf(size = 0.15, color = "black", aes(fill = vector)) +
                 labs(fill = vector.name) +
                 ggtitle(title) +
                 # This magically only adds axes:
                 geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
                 xlab("Longitude") +
                 ylab("Latitude")
  tdwg.plot
}


##########################################################
# Plots of FRic and FDis: observed and null model z-scores
##########################################################
cat("Plotting FRic distributions...\n")

if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
}

# Sample plotting code
# --------------------
if (FALSE) {
ggsave(plot = SpatialPlot(tdwg.map,
                          vector = fd.observed.mean[, "FRic.all.traits"],
                          vector.name = "FRic",
                          title = "Observed Functional Richness",
                          subtitle = "Genus-mean filled data"
                          ),
       filename = paste0(plot.dir,
                         "FRic_observed_mean",
                         ".png"
                         ),
       width = 8,
       height = 4
       )
}

# Generate all plots
# ------------------
# This we do with one master function that generates 16 graphs for the
# FD indices based on all traits together.
MakePlots <- function() {
  for (model in c("observed", "global", "regional", "local")) {
    for (index in c("FRic", "FDis")) {
      for (source in c("genus_mean", "stochastic_mean")) {
        if (index == "FRic") {
          index.long <- "Functional richness"
        } else {
          index.long <- "Functional dispersion"
        }
        if (source == "genus_mean") {
          objname <- paste0("fd.", model, ".mean")
        } else {
          objname <- paste0("fd.", model, ".stochastic")
        }
        ggsave(plot = SpatialPlot(tdwg.map,
                                  vector = get(objname)[, paste0(index,
                                                                 ".all.traits"
                                                                 )
                                                        ],
                                  vector.name = index,
                                  title = paste("Observed", index.long),
                                  subtitle = paste(source, "filled data")
                                  ),
               filename = paste0(plot.dir,
                                 index,
                                 "_", model, "_",
                                 source,
                                 ".png"
                                 ),
               width = 8,
               height = 4
               )
      }
    }
  }
  return (0)
}

MakePlots()


##################
# Additional plots
##################
cat("Creating additional plots...\n")

# Three Realms
# ------------
ggsave(plot = SpatialPlotFill(tdwg.map,
                              vector = tdwg3.info[, "realm"],
                              vector.name = "Realm",
                              title = "Assignment of TDWG3 units to realms"
                              ),
       filename = paste0(plot.dir,
                         "TDWG3_realm_distribution",
                         ".png"
                         ),
       width = 8,
       height = 4
       )

# Community average trait values
# ------------------------------
MakePlots <- function() {
  for (trait in c("stem.height", "blade.length", "fruit.length")) {
    for (source in c("genus_mean", "stochastic_mean")) {
      if (source == "genus_mean") {
        objname <- paste0("cwm.observed.mean")
      } else {
        objname <- paste0("cwm.observed.stochastic")
      }

  ggsave(plot = SpatialPlot(tdwg.map,
                            vector = get(objname)[, trait],
                            vector.name = paste0("log10(", trait, ")"),
                            title = paste("Community mean", trait)
                            ),
         filename = paste0(plot.dir,
                           "TDWG3_cwm_",
                           source,
                           "_",
                           trait,
                           ".png"
                           ),
         width = 8,
         height = 4
         )
    }
  }
}

MakePlots()


cat("Done.\n")

