############################################
# Palm FD project: custom plotting functions
############################################
#
# In which custom plotting functions using the ggplot2 package are defined.
# Many of these functions are designed for plotting spatial data.
#
# NOTE: 'simple features' type spatial data is used, with the R package 'sf'.
# With sf and ggplot2, instead of using geom_polygon(), we use geom_sf() with the
# 'sf' object. And instead of ggplot2::fortify(), we pretend the 'sf' object
# is a dataframe and simply add our data to it as a column.
# See the function SpatialPlot() below for an example.
# For info and tutorial see https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
#
#
# FUNCTION LIST:
# SpatialPlot
# SpatialPlotFill
# SpatialPlotSegments
# SpatialPlotFactor
# SpatialPlotNB
# TraitTrans
# NBreaks
# SpatialPlotBins


library(magrittr)
library(plyr)
library(sf)
library(ggplot2)
library(scales)
library(viridis)

source(file="functions/base_functions.R")


SpatialPlot <- function(tdwg.map, vector, vector.name, vector.size = NULL,
                        title = NULL, subtitle = NULL, legend.position = "right",
                        labels = waiver(), colors = "viridis") {
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
    color.scale <-
#      scale_fill_gradient(low = "yellow", high = "red", labels = labels)
       scale_fill_viridis(discrete = is.discrete(vector),
                          option = colors,
                          begin = 0.5,
                          direction = 1,
                          labels = labels
                          )
  } else {
    if (max(vector, na.rm = TRUE) <= 0) {
      color.scale <-
#        scale_fill_gradient(low = "blue", high = "cyan", labels = labels)
       scale_fill_viridis(discrete = is.discrete(vector),
                          option = colors,
                          begin = 0.5,
                          direction = -1,
                          labels = labels
                          )
    } else {
      color.scale <-
#        scale_fill_gradientn(colours = c("blue", "cyan", "white", "yellow", "red"),
#                             values = rescale(c(min(vector, na.rm = TRUE),
#                                                -1e-10,
#                                                0,
#                                                1e-10,
#                                                max(vector, na.rm = TRUE)
#                                                )
#                                              ),
#                             labels = labels
#                             )
       scale_fill_viridis(discrete = is.discrete(vector),
                          option = colors,
                          begin = 0.0,
                          direction = 1,
                          labels = labels
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
                    shape = 21
                    ) +
         labs(fill = vector.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              ) +
         color.scale +
         theme(legend.position = legend.position)
}


SpatialPlotFill <- function(tdwg.map, vector, vector.name, title = NULL,
                            subtitle = NULL, colors = "viridis",
                            direction = 1, begin = 0, legend.position = "right",
                            poly.col = "black", legend.labels = waiver()) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units)
# vector.name: character string giving name for the vector (used in legend)
  tdwg.map$vector <- vector
  tdwg.plot <- ggplot(data = tdwg.map) + 
                 geom_sf(size = 0.15, color = poly.col, aes(fill = vector)) +
                 # This magically only adds axes:
                 geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
                 labs(fill = vector.name,
                      x = NULL,
                      y = NULL,
                      title = title,
                      subtitle = subtitle
                      ) +
                 scale_fill_viridis(na.value = "#B0B0B0",
                                    discrete = is.discrete(vector),
                                    option = colors,
                                    direction = direction,
                                    begin = begin,
                                    labels = legend.labels
                                    ) +
                 theme(legend.position = legend.position)
  tdwg.plot
}


SpatialPlotSegments <- function(tdwg.map, segments = NULL, fill.vector, fill.name, 
                                title = NULL, subtitle = NULL, colors = "viridis",
                                direction = 1, begin = 0) {
# Segments: data.frame with coordinates for segments. See geom_segment().
#           The order of columns in the dataframe should be
#           c("x", "y", "xend", "yend")
  tdwg.plot <- ggplot(data = tdwg.map) + 
                 geom_sf(size = 0.15, color = "black", aes(fill = fill.vector)) +
                 # This magically only adds axes:
                 geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
                 labs(fill = fill.name,
                      x = NULL,
                      y = NULL,
                      title = title,
                      subtitle = subtitle
                      ) +
                 scale_fill_viridis(na.value = "#B0B0B0",
                                    discrete = is.discrete(fill.vector),
                                    option = colors,
                                    direction = direction,
                                    begin = begin
                                    )
  if (!is.null(segments)) {
    tdwg.plot <- tdwg.plot +
                 # Make segments
                 geom_segment(data = segments,
                              aes(x = segments[, 1],
                                  y = segments[, 2],
                                  xend = segments[, 3],
                                  yend = segments[, 4]
                                  ),
                              size = 0.5,
                              colour = "red",
                              show.legend = FALSE
                              )
  }
  tdwg.plot
}


SpatialPlotFactor <- function(tdwg.map, factor, factor.name, title = NULL,
                              subtitle = NULL) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units).
#         Data must be a continuous numeric variable.
# vector.name: character string giving name for the vector (used in legend)
  tdwg.map$factor <- as.factor(factor)
  subset <- tdwg.map[!is.na(tdwg.map$factor), ]
  ggplot(data = tdwg.map) + 
         geom_sf(size = 0.15, color = "black") +
         # This magically only adds axes:
         geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
         # While this does NOT add axes but does add points:
         geom_point(data = subset,
                    mapping = aes(x = Long,
                                  y = Lat,
                                  fill = factor
                                  ),
                    size = 3,
                    shape = 21) +
         labs(fill = factor.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              )
}



SpatialPlotNB <- function(tdwg.map, presence, segments, title = NULL,
                          subtitle = NULL) {
# Function to visualize a neighbourhood.
# Args:
#   tdwg.map: the 'sf' object with all polygons to be plotted in the background
#   presence: vector by which to color (fill) the polygons
#   Segments: data.frame with coordinates for segments. See geom_segment().
#             You likely want to use the function Nb2segments()

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
# Points disabled because they clutter the graph too much. Lines should be enough
#                 geom_point(data = tdwg.map[which(presence), ],
#                            aes(x = Long, y = Lat),
#                            pch = 21,
#                            size = 2
#                            ) +
                 geom_segment(data = segments,
                              aes(x = segments[, 1],
                                  y = segments[, 2],
                                  xend = segments[, 3],
                                  yend = segments[, 4]
                                  ),
                              size = 0.35,
                              show.legend = FALSE
                              ) +
                 scale_fill_viridis(na.value = "#B0B0B0",
                                    discrete = is.discrete(presence),
                                    option = "inferno",
                                    begin = 0.75
                                    )
  tdwg.plot
}



TraitTrans <- function(breaks) {
# Transform legend or axis tick marks of log-transformed trait values
# back to real values.
  round((10 ^ breaks) - 1, digits = 1)
}


NBreaks <- function(limits, n = 5) {
# Create a custom number of equidistant breaks.
# For use with e.g. scale_x_continuous()
# IF the range encompasses 0, 0 is added as a break because the origin is important,
# and breaks too close to 0 will be removed
  range <- max(limits) - min(limits)
  margin <- range * 0.05
  dist <- (range - 2 * margin) / (n - 1)
  breaks <- (min(limits) + margin) + (0:(n - 1)) * dist
  if (any(breaks < 0) & any(breaks > 0)) {
    remove.ind <- breaks > (0 - margin * 2) & breaks < (0 + margin * 2)
    breaks <- breaks[!remove.ind]
    breaks <- sort(c(breaks, 0))
  }
  round(breaks, digits = 2)
}


SpatialPlotBins <- function(tdwg.map, vector, vector.name, vector.size = NULL,
                            title = NULL, subtitle = NULL, legend.position = "right",
                            labels = "default", colors = "viridis", gfx.scale = 1) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units).
#         Data must be a continuous numeric variable.
# vector.size: vector of length(vector) giving size of points. If NULL,
#              size will be scaled by the values of 'vector'
# vector.name: character string giving name for the vector (used in legend)
# labels: Legend labels. By default, labels will be names of bins.
#         Can be a character vector giving legend labels.
#         Can be a function that takes a numeric vector as input, which will be
#         applied to vector quantiles.
  tdwg.map$vector <- vector
  vector.quantiles <-
    quantile(vector, probs = seq(0, 1, 0.125), na.rm = TRUE)
  vector.bins <- cut(vector, vector.quantiles, include.lowest = TRUE)
  tdwg.map$bins <- vector.bins

  if (!is.null(vector.size)) {
    tdwg.map$vector.size <- vector.size
    subset <- tdwg.map[!is.na(tdwg.map$vector), ]
    size <- rescale(subset$vector.size, to = c(0.25, 1)) * 3.5 * gfx.scale
  } else {
    subset <- tdwg.map[!is.na(tdwg.map$vector), ]
    size <- 3 * gfx.scale
  }
  legend.size <- rescale(vector.quantiles[-1], to = c(0.25, 1)) * 4.5 * gfx.scale

  if (identical(labels, "default")) {
    x <- vector.quantiles
    lower <- round(x[1:(length(x) - 1)], digits = 2)
    upper <- round(x[2:length(x)], digits = 2)
    plot.labels <- paste(lower, "-", upper)
  } else {
    if (is.function(labels)) {
      x <- do.call(labels, list(vector.quantiles))
      lower <- round(x[1:(length(x) - 1)], digits = 2)
      upper <- round(x[2:length(x)], digits = 2)
      plot.labels <- paste(lower, "-", upper)
    } else {
      plot.labels <- labels
    }
  }

  color.scale <- scale_fill_viridis(discrete = TRUE,
                                    option = colors,
                                    labels = plot.labels,
                                    breaks = levels(vector.bins)
                                    )

  ggplot(data = tdwg.map) + 
         geom_sf(size = 0.15 * gfx.scale, color = "white", fill = "#B0B0B0") +
         # This magically only adds axes:
         geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
         # While this does NOT add axes but does add points:
         geom_point(data = subset,
                    mapping = aes(x = Long, y = Lat, fill = bins),
                    size = size,
                    shape = 21
                    ) +
         labs(fill = vector.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              ) +
         color.scale +
         theme(legend.position = legend.position,
               legend.margin = margin(t = 0 * gfx.scale,
                                      r = 1 * gfx.scale,
                                      b = 1 * gfx.scale,
                                      l = 1 * gfx.scale,
                                      unit = "pt"
                                      ),
               legend.box.spacing = unit(0.1 * gfx.scale, "cm")
               ) +
         guides(fill = guide_legend(override.aes = list(size = legend.size)))
}

