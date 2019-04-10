###############################################################
# Palm FD project: Generating graphs and plots for presentation
###############################################################
#
# In which we generate plots and figures of the project, suitable for presentation
#
# Input files:
#  stuff
# Generated output files:
#  pictures


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(sf)
library(ggplot2)
library(scales)
library(viridis)

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")

# Directory for saving plots (with a trailing slash):
plot.dir = "graphs/Presentation_plots/"
# Create this directory
if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
}



#######################
# Prepare plot datasets
#######################
cat("Loading data...\n")

# TDWG3 level info
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")
tdwg3.info[, "palm.richness"][tdwg3.info[, "palm.richness"] == 0] <- NA

# Spatial polygons
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")


# Parsing function for palm data
Temp <- function(x) {
  x[, "tdwg3.code"] <- rownames(x)
  merge(x = tdwg3.info[, "tdwg3.code", drop = FALSE],
        y = x,
        by = "tdwg3.code",
        all.x = TRUE,
        sort = TRUE
        )
}

# Community mean trait values
cwm <- read.csv(file = "output/observed_FD/community_trait_means_stochastic_mean_of.csv",
                                    row.names = 1
                                    ) %>%
  Temp()

# Assemblage dispersion field weights
adf.weights <-
  read.csv(file = "output/adf_weight_matrix.csv",
           header = TRUE,
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()


# Calculated functional diversity indices
fd.global <-
  read.csv(file = "output/null_model_global/FD_z_scores_global_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()

fd.realm <-
  read.csv(file = "output/null_model_regional/FD_z_scores_regional_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()

fd.ADF <-
  read.csv(file = "output/null_model_local/FD_z_scores_local_stochastic_mean_of_z.scores.csv",
           row.names = 1
           ) %>%
  Temp()
# Clean infinities from the data (there can be infinities for HAW, NFK, NWC and SEY
# because these have 100% endemic palm communities)
fd.ADF <-
  lapply(fd.ADF, function(x) { replace(x, is.infinite(x), NA) } ) %>%
  as.data.frame()









#####################################################
cat("Defining universal spatial plotting function\n")
#####################################################
# A function which combines filling polygons and plotting points on top of
# polygons.
MultiSP <- function(tdwg.map, fill.vector, fill.name, point.vector, point.name,
                    point.size = NULL, title = NULL, subtitle = NULL) {
# Multi Spatial Plot.
# Args:
#  tdwg.map: sf object with spatial polygon information
#  fill.vector: vector with data by which to fill polygons
#  fill.name: character string naming the 'fill.vector' variable (for legend)
#  point.vector: vector with data by which to color points plotted on polygons
#  point.name: character string naming the 'point.vector' variable (for legend)
#  Point.size: optional vector of the same length as point.vector, giving the
#              the relative size of each point.
#  title: overall plot title.
#  subtitle: overall plot subtitle.
  tdwg.map[, "fill.vector"] <- fill.vector
  tdwg.map[, "points"] <- point.vector
  if (!is.null(point.size)) {
    tdwg.map[, "point.size"] <- point.size
    subset <- tdwg.map[!is.na(tdwg.map[, "points"]), "point.size"]
    size <- rescale(subset) * 10
  } else {
    size <- 3
  }
  subset <- tdwg.map[!is.na(point.vector), ]
  fill.vector.cols <-
    viridis_pal(option = "viridis", begin = 0.25) (length(fill.vector))
  fill.vector.cols[is.na(fill.vector)] <- "#C0C0C0"

  ggplot(data = tdwg.map) + 
    geom_sf(size = 0.15,
            color = "black",
            show.legend = FALSE
            ) +
    scale_fill_manual(values = fill.vector.cols)

    # This magically only adds axes:
    geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
    # While this does NOT add axes but does add points:
    geom_point(data = subset,
               mapping = aes(x = Long, y = Lat, fill = points),
               size = size,
               shape = 21
               ) +
    scale_fill_viridis(discrete = is.discrete(point.vector),
                       option = "plasma",
                       direction = -1
                       ) +
    labs(fill = point.name,
         x = NULL,
         y = NULL,
         title = title,
         subtitle = subtitle
         )
}

if (FALSE) {
fill.vector <- tdwg3.info[, "experimental.richness"] > 0
#fill.vector <- runif(nrow(tdwg3.info))

MultiSP(tdwg.map,
        fill.vector = fill.vector,
        fill.name = "filler",
        point.vector = log(tdwg3.info[, "palm.richness"]),
        point.name = "log(richness)",
        point.size = NULL,
        title = "Palm species richness",
        subtitle = "In botanical countries (TDWG3 units)"
        )

# magma, inferno, plasma, viridis (default), cividis

# You don't actually need ggplot:
plot(st_geometry(tdwg.map), col = c("green", "red"), lwd = 0.5)
# the underlying plotting method is 'plot_sf'
points(Lat ~ Long, data = tdwg.map[1:100, ], pch = 21, bg = "blue")
axis(1, at = c(-70, 70), lwd = 0, lwd.ticks = 1)

}




###################################
cat("Palm species richness plots\n")
###################################

richness.all <-
  SpatialPlotFill(tdwg.map,
                  vector = log(tdwg3.info[, "palm.richness"]),
                  vector.name = "log(richness)",
                  title = "Palm species richness",
                  subtitle = "In botanical countries (TDWG3 units)"
                  )
ggsave(richness.all,
       filename = paste0(plot.dir, "01_palm_richness_all.png"),
       width = 8,
       height = 4
       )

################################################
cat("Palm community trait distribution plots\n")
################################################

MakePlots <- function() {
  for (trait in c("stem.height", "blade.length", "fruit.length")) {
    p <- SpatialPlot(tdwg.map,
                     vector = cwm[, trait],
                     vector.name = paste0("log10(", trait, ")"),
                     title = NULL,
                     subtitle = NULL
                     )
    ggsave(p,
           filename = paste0(plot.dir, "TDWG3_cwm_", trait, ".png"),
           width = 8,
           height = 4
           )
  }
}

MakePlots()


################################
cat("Three realms definition\n")
################################

realm <- tdwg3.info$realm
realm[is.na(tdwg3.info$experimental.richness)] <- NA

threerealm <- 
  SpatialPlotFill(tdwg.map,
                  vector = realm,
                  vector.name = "Realm",
                  title = "Division of botanical countries by realm",
                  subtitle = NULL,
                  colors = "plasma",
                  begin = 0.25
                  )
ggsave(threerealm,
       filename = paste0(plot.dir, "three_realm.png"),
       width = 8,
       height = 4
       )


#####################################
cat("Species pools of null models\n")
#####################################

# Global map
global <- as.character(realm)
global[!is.na(global)] <- "global"

global.map <-
  SpatialPlotFill(tdwg.map,
                  vector = global,
                  vector.name = "global",
                  title = NULL,
                  subtitle = NULL,
                  colors = "plasma",
                  begin = 0.5
                  )
ggsave(global.map,
       filename = paste0(plot.dir, "species_pool_global.png"),
       width = 8,
       height = 4
       )

# Realm map is what we have shown before

# ===================================
# ADF maps: we need multiple examples
# ===================================

# Functions to extract telemetry:
# -------------------------------
ParseADF <- function(tdwg3.code, coverage) {
# specify a single tdwg3 unit code
  subset <- adf.weights[rownames(adf.weights) %in% tdwg3.code, ]
  if (all.equal(coverage, 1)) {
    coverage <- 0.9999
  }

  fun <- function(x, cover = coverage) {
    # keep n highest values, by setting all other values to 0.
    # We exploit the presence of (column) names.
    # The output columns will be in alphabetical order.
    # In case of tied values at the cutoff, ALL of these values are
    # included. I.e., the rowSums will be at least equal to coverage, but
    # could be higher.
    sorted <- sort(x, decreasing = TRUE)
    sorted.sum <- cumsum(sorted)
    n <- min(which(sorted.sum > cover))
    min.to.keep <- min(x[WhichMax(x, n = n)])
    sorted[sorted < min.to.keep] <- 0
    sorted[order(names(sorted))]
  }
  adf.weights.pruned <- fun(subset)

  ind <- adf.weights.pruned > 0
  tdwg3.ends <- names(adf.weights.pruned)[ind]
  end.coords <- tdwg3.info[tdwg3.info[, "tdwg3.code"] %in% tdwg3.ends,
                           c("lat", "lon")
                           ]
  xy.ind <- match(tdwg3.code, tdwg3.info[, "tdwg3.code"])
  segments <- data.frame(x    = rep(tdwg3.info[xy.ind, "lon"], nrow(end.coords)),
                         y    = rep(tdwg3.info[xy.ind, "lat"], nrow(end.coords)),
                         xend = end.coords[, "lon"],
                         yend = end.coords[, "lat"]
                         )
  segments
}

GetFill <- function(segments, tdwg3.code) {
  fill.vector <- as.numeric(rep(NA, nrow(tdwg3.info)))
  indices <- tdwg3.info[, "lon"] %in% segments[, "xend"]
  tdwg <- as.character(tdwg3.info[indices, "tdwg3.code"])
  fill.vector[indices] <- adf.weights[tdwg3.code, tdwg]
  fill.vector
}

# Plotting the dispersion fields
# ------------------------------
tdwg3.codes <- c("ZAI", "MDG", "BZL", "CLM", "IND", "NWG")
for (area in tdwg3.codes) {
  name <- as.character(tdwg3.info[match(area, tdwg3.info$tdwg3.code), "tdwg3.name"])
  p <- SpatialPlotSegments(tdwg.map,
                           segments = NULL,
                           fill.vector = GetFill(ParseADF(area, 1), area),
                           fill.name = "Sampling chance",
                           title = paste("Dispersion field of", name),
                           subtitle = NULL,
                           begin = 0.25
                           )
  ggsave(p,
         filename = paste0(plot.dir, "species_pool_adf_", area, ".png"),
         width = 8,
         height = 4
         )
}


#################################################
cat("Maps of global palm functional diversity\n")
#################################################

MakePlots <- function() {
  for (model in c("global", "realm", "ADF")) {
    for (index in c("FRic", "FDis")) {

      if (index == "FRic") {
        index.long <- "Functional Richness"
      } else {
        index.long <- "Functional Dispersion"
      }
      objname <- paste0("fd.", model)

      ggsave(plot = SpatialPlot(tdwg.map,
                                vector = get(objname)[, paste0(index,
                                                               ".all.traits"
                                                               )
                                                      ],
                                vector.name = index,
                                title = paste(index.long, "- null model", model),
                                subtitle = NULL
                                ),
             filename = paste0(plot.dir,
                               "FD_global_",
                               index,
                               "_",
                               model,
                               ".png"
                               ),
             width = 8,
             height = 4
             )
    }
  }
  return (0)
}

MakePlots()




cat("Done.\n")

