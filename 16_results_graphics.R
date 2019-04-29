########################################
# Palm FD project: Graphs of all results
########################################
#
# One single script to generate all relevant graphs of the final results.
#
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(plyr)
library(magrittr)
library(car)
library(sf)
library(ggplot2)
library(scales)
library(viridis)
library(corrplot)
library(leaps)
library(gridExtra)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")
source(file = "functions/OLS_regression_functions.R")

# Graphics output directory (with trailing slash!)
plot.dir <- "graphs/main_results/"


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


# Functional diversity indices
# ----------------------------
# A list with all the FD index data.
# Note these dataframes include richness and endemism variables.
fd.indices <- vector("list", length = 0)
fric <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
          "FRic.fruit.length")
fdis <- c("FDis.all.traits", "FDis.stem.height", "FDis.blade.length",
          "FDis.fruit.length")
fd.names <- c(fric, fdis)
for (index in fd.names) {
  fd.indices [[index]] <-
    read.csv(file = paste0("output/FD_summary_", index, ".csv"), row.names = 1)
}


# Community trait means
# ---------------------
Temp <- function(x) {
  x[, "tdwg3.code"] <- rownames(x)
  merge(x = tdwg3.info[, "tdwg3.code", drop = FALSE],
        y = x,
        by = "tdwg3.code",
        all.x = TRUE,
        sort = TRUE
        )
}
cwm.observed <- 
  read.csv(file = "output/observed_FD/community_trait_means_stochastic_mean_of.csv",
           row.names = 1
           ) %>%
  Temp()
# Keep in mind these are log10-transformed (after adding one):
#   trans = log10(data + 1)
# The reverse operation would be
#   original = (10 ^ trans) - 1


# Environmental data
# ------------------
env.complete <- read.csv(file = "output/tdwg3_predictors_complete.csv",
                         row.names = 1
                         )



########################################################
cat("Plotting correlograms of predictor variables...\n")
########################################################

# 1. Full set of predictors, including species, endemism and realm
GraphSVG(Correlogram(env.complete),
         file = paste0(plot.dir, "Correlogram_predictors_ALL.svg"),
         width = 7,
         height = 7
         )

# 2. All environmental predictors
exclude <- c("palm.richness", "endemism", "realm")
GraphSVG(Correlogram(env[, !colnames(env.complete) %in% exclude]),
         file = paste0(plot.dir, "Correlogram_predictors_ENV.svg"),
         width = 7,
         height = 7
         )

# 3. Environmental predictors excluding collinear variables,
#    i.e. without alt_range
exclude <- c("palm.richness", "endemism", "realm", "alt_range")
GraphSVG(Correlogram(env[, !colnames(env.complete) %in% exclude]),
         file = paste0(plot.dir, "Correlogram_predictors_ENV_noncoll.svg"),
         width = 7,
         height = 7
         )



#######################################################
cat("Plotting community mean trait distributions...\n")
#######################################################

MultiTraitPlot <- function(tdwg.map, cwm.observed) {

  LabTrans <- function(breaks) {
    round((10 ^ breaks) - 1, digits = 1)
  }
  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  height <- SpatialPlot(tdwg.map = tdwg.map[no.ANT, ],
                        vector = cwm.observed[no.ANT, "stem.height"],
                        vector.name = "Stem height (m)",
                        title = "Mean stem height in palm communities",
                        legend.position = "bottom",
                        labels = LabTrans
                        ) +
            theme(legend.key.width = unit(1.5, "cm"))
  blade <- SpatialPlot(tdwg.map = tdwg.map[no.ANT, ],
                        vector = cwm.observed[no.ANT, "blade.length"],
                        vector.name = "Blade length (m)",
                        title = "Mean blade length in palm communities",
                        legend.position = "bottom",
                        labels = LabTrans
                        ) +
            theme(legend.key.width = unit(1.5, "cm"))
  fruit <- SpatialPlot(tdwg.map = tdwg.map[no.ANT, ],
                        vector = cwm.observed[no.ANT, "fruit.length"],
                        vector.name = "Fruit length (cm)",
                        title = "Mean fruit length in palm communities",
                        legend.position = "bottom",
                        labels = LabTrans
                        ) +
            theme(legend.key.width = unit(1.5, "cm"))

  arrangeGrob(height, blade, fruit, ncol = 1)
}

ggsave(plot = MultiTraitPlot(tdwg.map, cwm.observed),
       filename = paste0(plot.dir, "_TDWG3_palm_cwm_traits.png"),
       width = 7,
       height = 12
       )



#######################################################
cat("Plotting functional diversity distributions...\n")
#######################################################

# First the all-traits data. 
# A single huge grid. Two columns: FDis and FRic.
# Five rows: observed, global, realm, realmnoMDG, ADF.
MultiFDPlot <- function(tdwg.map, fd.indices) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"
  fdis.observed <- SpatialPlot(tdwg.map = tdwg.map[no.ANT, ],
                               vector = GetFD(fd.indices["FDis.all.traits"],
                                              "observed"
                                              ),
                               vector.name = "FDis (observed)",
                               title = "Functional Dispersion (observed)",
                               legend.position = "bottom",
                               ) +
                   theme(legend.key.width = unit(1.5, "cm"))




cat("Done.\n")

