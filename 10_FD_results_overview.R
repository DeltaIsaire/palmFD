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


cat("Preparing data...\n")
#####################
# Info on tdwg3 units
#####################
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


# For comparison: Palm community percent endemism
# -----------------------------------------------
cat("Calculating endemism percentage of palm communities...\n")

endemism <- numeric(length = nrow(pres.abs.matrix))
for (row in seq_len(nrow(pres.abs.matrix))) {
  sp.indices <- pres.abs.matrix[row, ] > 0
  total.species = sum(sp.indices)

  subset <- pres.abs.matrix[-row, sp.indices]
  endemics <- sum(colSums(subset) == 0)
  endemism[row] <- endemics / total.species
}
tdwg3.info <-
  merge(x = tdwg3.info,
        y = data.frame(tdwg3.code = rownames(pres.abs.matrix),
                       endemism   = endemism
                       ),
        by = "tdwg3.code",
        all.x = TRUE
        )

# Save new tdwg3.info
# -------------------
write.csv(tdwg3.info,
          file = "output/tdwg3_info_v2.csv",
          eol = "\r\n",
          row.names = FALSE
          )


###################
# Load spatial data
###################
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")
# Note I am using the 'sf' spatial format instead of 'SpatialPolygonsDataFrame'.
# The sf package with function sf::read_sf is more modern than using rgdal::readOGR.
# That modernity translates into improved performance.


####################
# Load FD index data
####################
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
# Concurrently, I have assessed the shape of the FD distributions, and presence
# of outliers, using Histogram() and Multihist().


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
# FRic: no exceptionally extreme values. Single-trait FRic is roughly normally
#       distributed, although stem height and blade length are somewhat bimodal.
#       All-traits FRic is strongly one-tailed.
# FDis: no exceptionally extreme values. Both single-trait and all-traits FRic
#       is roughly normally distributed.
# Mean and Stochastic data are not qualitatively different.


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
# FRic: no exceptionally extreme values. Both single-trait and all-traits FRic
#       are roughly normal, although blade length is bimodal.
# FDis: no exceptionally extreme values. Both single-triat and all-traits FDis are
#       roughly normal.
# Mean and stochastic data are not qualitatively different.


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
# FRic: no exceptionally extreme values. Both all-traits and single-traits FRic
#       are roughly normally distributed, but stem height and blade length are
#       somewhat left-tailed.
# FDis: Both all-traits and single-traits FRic are roughly normally distributed.
#       there are a few especially low values for all-traits FDis, for CUB and
#       MDG (Cuba and Madagascar). Values not problematically extreme, just notable.
# Mean and stochastic data are not qualitatively different.


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
# Clean infinities from the data (there can be infinities for HAW, NFK, NWC and SEY
# because these have 100% endemic palm communities)
fd.local.mean <-
  lapply(fd.local.mean, function(x) { replace(x, is.infinite(x), NA) } ) %>%
  as.data.frame()
fd.local.stochastic <-
  lapply(fd.local.stochastic, function(x) { replace(x, is.infinite(x), NA) } ) %>%
  as.data.frame()
# Fric: both single-trait and all-traits FRic are roughly normally distributed.
#       There are a few very low values for stem height: CUB and MDG (Cuba  and
#       Madagscar) and all.traits: MDG (Madagascar).
#       In the stochastic data, there is also an extreme low value for blade length,
#       for MDG (Madagascar).
# FDis: both single-trait and all-traits FDis are roughly normally distributed.
#       There are a few very low values for fruit length and all.traits, for both
#       mean and stochastic datasets. These values belong to MDG (Madagascar).


# Community mean trait values
# ---------------------------
# Keep in mind these are log10-transformed
cwm.observed.mean <- read.csv(file = "output/observed_FD/community_trait_means_genus_mean.csv",
                              row.names = 1
                              ) %>%
  Temp()
cwm.observed.stochastic <- read.csv(file = "output/observed_FD/community_trait_means_stochastic_mean_of.csv",
                                    row.names = 1
                                    ) %>%
  Temp()
# These are all roughly normally distributed, without extreme values, except for
# a rather large mean fruit size in MLI (Mali).


# Summary of normality and outlier assessment
# -------------------------------------------
# Overall, it looks like no transformations are necessary to use these data in
# regression models. 
# CUB and MDG (cuba and Madagascar) have unusually low functional diversity,
# according to some metrics. These low values are real, not erroneous.


#############################
# Summary files of FD indices
#############################
cat("Saving summary file of FD indices...\n")
# Create and save a dataframe with all all-traits FRic values for each TDWG3 unit,
# including null model means and SDs. Using data from stochastic gapfilling.

# Load null model means and SDs:
global.means <-
  read.csv(file = "output/null_model_global/FD_z_scores_global_stochastic_mean_of_means.csv",
           row.names = 1
           ) %>%
  Temp()
global.sds <-
  read.csv(file = "output/null_model_global/FD_z_scores_global_stochastic_mean_of_sds.csv",
           row.names = 1
           ) %>%
  Temp()
regional.means <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_means.csv",
           row.names = 1
           ) %>%
  Temp()
regional.sds <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_sds.csv",
           row.names = 1
           ) %>%
  Temp()
local.means <-
  read.csv(file = "output/null_model_local/FD_z_scores_local_stochastic_mean_of_means.csv",
           row.names = 1
           ) %>%
  Temp()
local.sds <-
  read.csv(file = "output/null_model_local/FD_z_scores_local_stochastic_mean_of_sds.csv",
           row.names = 1
           ) %>%
  Temp()

summary.FRic <-
  data.frame(tdwg3.name       = tdwg3.info[, "tdwg3.name"],
             tdwg3.realm      = tdwg3.info[, "realm"],
             palm.richness    = tdwg3.info[, "experimental.richness"],
             percent.endemism = tdwg3.info[, "endemism"],
             FRic.observed    = fd.observed.stochastic[, "FRic.all.traits"],
             FRic.global.SES  = fd.global.stochastic[, "FRic.all.traits"],
             FRic.global.mean = global.means[, "FRic.all.traits"],
             FRic.global.sd   = global.sds[, "FRic.all.traits"],
             FRic.realm.SES   = fd.regional.stochastic[, "FRic.all.traits"],
             FRic.realm.mean  = regional.means[, "FRic.all.traits"],
             FRic.realm.sd    = regional.sds[, "FRic.all.traits"],
             FRic.adf.SES     = fd.local.stochastic[, "FRic.all.traits"],
             FRic.adf.mean    = local.means[, "FRic.all.traits"],
             FRic.adf.sd      = local.sds[, "FRic.all.traits"],
             row.names = fd.observed.stochastic[, "tdwg3.code"]
             )
write.csv(summary.FRic,
          file = "output/FD_summary_FRic.csv",
          eol = "\r\n",
          row.names = TRUE
          )

# And the same thing for FDis:
summary.FDis <-
  data.frame(tdwg3.name       = tdwg3.info[, "tdwg3.name"],
             tdwg3.realm      = tdwg3.info[, "realm"],
             palm.richness    = tdwg3.info[, "experimental.richness"],
             percent.endemism = tdwg3.info[, "endemism"],
             FDis.observed    = fd.observed.stochastic[, "FDis.all.traits"],
             FDis.global.SES  = fd.global.stochastic[, "FDis.all.traits"],
             FDis.global.mean = global.means[, "FDis.all.traits"],
             FDis.global.sd   = global.sds[, "FDis.all.traits"],
             FDis.realm.SES   = fd.regional.stochastic[, "FDis.all.traits"],
             FDis.realm.mean  = regional.means[, "FDis.all.traits"],
             FDis.realm.sd    = regional.sds[, "FDis.all.traits"],
             FDis.adf.SES     = fd.local.stochastic[, "FDis.all.traits"],
             FDis.adf.mean    = local.means[, "FDis.all.traits"],
             FDis.adf.sd      = local.sds[, "FDis.all.traits"],
             row.names = fd.observed.stochastic[, "tdwg3.code"]
             )
write.csv(summary.FDis,
          file = "output/FD_summary_FDis.csv",
          eol = "\r\n",
          row.names = TRUE
          )


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
SpatialPlot <- function(tdwg.map, vector, vector.name, vector.size = NULL,
                        title = NULL, subtitle = NULL) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units).
#         Data must be a continuous numeric variable.
# vector.name: character string giving name for the vector (used in legend)
  tdwg.map$vector <- vector
  if (!is.null(vector.size)) {
    tdwg.map$vector.size <- vector.size
    subset <- tdwg.map[!is.na(tdwg.map$vector), ]
    size <- rescale(subset$vector.size) * 10
  } else {
    subset <- tdwg.map[!is.na(tdwg.map$vector), ]
    size <- 3
  }
  if (min(vector, na.rm = TRUE) >= 0) {
    color.scale <- scale_fill_gradient(low = "yellow", high = "red")
  } else {
    if (max(vector, na.rm = TRUE) <= 0) {
      color.scale <- scale_fill_gradient(low = "blue", high = "cyan")
    } else {
      color.scale <-
        scale_fill_gradientn(colours = c("blue", "cyan", "white", "yellow", "red"),
                             values = rescale(c(min(vector, na.rm = TRUE),
                                                -1e-10,
                                                0,
                                                1e-10,
                                                max(vector, na.rm = TRUE)
                                                )
                                              )
                             )
    }
  }

  ggplot(data = tdwg.map) + 
         geom_sf(size = 0.15, color = "black") +
         # This magically only adds axes:
         geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
         # While this does NOT add axes but does add points:
         geom_point(data = subset,
                    mapping = aes(x = Long, y = Lat, fill = vector),
                    size = size,
                    shape = 21) +
         labs(fill = vector.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              ) +
         color.scale
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
                 xlab(NULL) +
                 ylab(NULL)
  tdwg.plot
}


##########################################################
# Plots of FRic and FDis: observed and null model z-scores
##########################################################
cat("Plotting Functional Diversity distributions...\n")

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
                                  title = paste(model, index.long),
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
                            title = paste("Community mean", trait),
                            subtitle = paste(source, "filled data")
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

# Community endemism proportion
# -----------------------------
ggsave(plot = SpatialPlot(tdwg.map,
                          vector = tdwg3.info[, "endemism"],
                          vector.name = "endemism",
                          vector.size = tdwg3.info[, "experimental.richness"],
                          title = "Proportion of endemic palm species",
                          subtitle = "(points scaled by species richness)"
                          ),
       filename = paste0(plot.dir,
                         "TDWG3_community_endemism_scaled",
                         ".png"
                         ),
       width = 8,
       height = 4
       )

ggsave(plot = SpatialPlot(tdwg.map,
                          vector = tdwg3.info[, "endemism"],
                          vector.name = "endemism",
                          title = "Proportion of endemic palm species",
                          subtitle = " "
                          ),
       filename = paste0(plot.dir,
                         "TDWG3_community_endemism",
                         ".png"
                         ),
       width = 8,
       height = 4
       )


cat("Done.\n")

