##########################################################
# Palm FD project: Summary and visualization of FD results
##########################################################
#
# In which the calculated FD indices (both observed and null model z-scores)
# are summarized and visualized.


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(sf)
library(ggplot2)
library(scales)

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")

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

# There is also data for an extra run, where MDG (Madagascar) was removed, because
# it has a large influence on the OWWest realm.
fd.regional.stochastic.noMDG <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_z.scores_noMDG.csv",
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
# Create and save a dataframe with all FRic and FDis values for each TDWG3 unit,
# including null model means and SDs. Using data from stochastic gapfilling.

# Load null model means and SDs
# -----------------------------
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
regional.means.noMDG <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_means_noMDG.csv",
           row.names = 1
           ) %>%
  Temp()
regional.sds.noMDG <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_sds_noMDG.csv",
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

# A function to do the hard work:
# -------------------------------
FdSummary <- function(index) {
# Function to build a summary dataframe for each FD index, with all the null model
# data for that index. The dataframe is also saved to disk.
#
# x: Fd Index for which to build the summary. A character string,
#    e.g. "FRic.all.traits"
  summary.index <-
    data.frame(tdwg3.name       = tdwg3.info[, "tdwg3.name"],
               tdwg3.realm      = tdwg3.info[, "realm"],
               palm.richness    = tdwg3.info[, "experimental.richness"],
               percent.endemism = tdwg3.info[, "endemism"],
               observed         = fd.observed.stochastic[, index],
               global.SES       = fd.global.stochastic[, index],
               global.mean      = global.means[, index],
               global.sd        = global.sds[, index],
               realm.SES        = fd.regional.stochastic[, index],
               realm.mean       = regional.means[, index],
               realm.sd         = regional.sds[, index],
               realm.SES.noMDG  = fd.regional.stochastic.noMDG[, index],
               realm.mean.noMDG = regional.means.noMDG[, index],
               realm.sd.noMDG   = regional.sds.noMDG[, index],
               adf.SES          = fd.local.stochastic[, index],
               adf.mean         = local.means[, index],
               adf.sd           = local.sds[, index],
               row.names        = fd.observed.stochastic[, "tdwg3.code"]
               )
  write.csv(summary.index,
          file = paste0("output/FD_summary_", index, ".csv"),
          eol = "\r\n",
          row.names = TRUE
          )
  summary.index
}

# Apply function to all FD indices
for (index in c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
                "FRic.fruit.length", "FDis.all.traits", "FDis.stem.height",
                "FDis.blade.length", "FDis.fruit.length"
                )
     ) {
  cat(index, "\n")
  FdSummary(index)
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


# Community similarity based on ADF weight matrix
# -----------------------------------------------
adf.weights <-
  read.csv(file = "output/adf_weight_matrix.csv",
           header = TRUE,
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()

# Plotting all ADF connections is too much. We need a cutoff, eliminating
# connections with only a few shared species. 
# Method: subset each ADF to only those communities accounting for x% of the total
#         probability mass.
# In addition, plotting all ADFs at once results in a complete mess, so subset
# to only a few chosen TDWG3 units.
include <- c("BOR", "NWG", "CLM", "MDG", "BZN", "ZAI")
adf.subset <- adf.weights[rownames(adf.weights) %in% include, ]

adf.weights.pruned <-
  aaply(adf.subset,
        1,
        function(x, coverage = 0.9) {
          # keep n highest values, by setting all other values to 0.
          # We exploit the presence of (column) names.
          # The output columns will be in alphabetical order.
          # In case of tied values at the cutoff, ALL of these values are
          # included. I.e., the rowSums will be at least equal to coverage, but
          # could be higher.
          sorted <- sort(x, decreasing = TRUE)
          sorted.sum <- cumsum(sorted)
          n <- min(which(sorted.sum > coverage))
          min.to.keep <- min(x[WhichMax(x, n = n)])
          sorted[sorted < min.to.keep] <- 0
          sorted[order(names(sorted))]
       }
       )

# Compile dataframe with the data for segments, based on the ADF matrix.
# We need x, y, xend and yend, all in lon/lat coordinates.
tdwg3.ends <- alply(adf.weights.pruned, 1, function(x) {
                                             ind <- x > 0
                                             names(x)[ind]
                                             }
                    )
end.coords <-
  ldply(tdwg3.ends,
        function(x) {
          tdwg3.info[tdwg3.info[, "tdwg3.code"] %in% x, c("lat", "lon")]
        }
        )
xy.indices <- match(end.coords[, "X1"], tdwg3.info[, "tdwg3.code"])
segments <- data.frame(x    = tdwg3.info[xy.indices, "lon"],
                       y    = tdwg3.info[xy.indices, "lat"],
                       xend = end.coords[, "lon"],
                       yend = end.coords[, "lat"]
                       )

# And create the plot.
# Here, the segment widths are NOT weighted. It's just the connections between
# TDWG3 units.
if (FALSE) {
# UPDATE: this code is obsolete. See presentation_graphics.R for an updated version.
ggsave(plot = SpatialPlotSegments(tdwg.map,
                                  segments = segments,
                                  title = "Spatial structure of Assemblage Dispersion Fields"
                                  ),
       filename = paste0(plot.dir,
                         "TDWG3_adf_segments",
                         ".png"
                         ),
       width = 8,
       height = 4
       )
}


cat("Done.\n")

