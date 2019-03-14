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
# The geom_sf() function is included in the github version of ggplot2, but not
# in the CRAN version. 
# For info and tutorial see https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
#
#
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
# SpatialPlot
# SpatialPlotFill
# SpatialPlotSegments
# SpatialPlotFactor


library(magrittr)
library(plyr)
library(sf)
library(ggplot2)
library(scales)

source(file="functions/base_functions.R")


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


SpatialPlotSegments <- function(tdwg.map, segments, title = NULL, subtitle = NULL) {
# Segments: data.frame with coordinates for segments. See geom_segment().
#           The order of columns in the dataframe should be
#           c("x", "y", "xend", "yend")
  tdwg.plot <- ggplot(data = tdwg.map) + 
                 geom_sf(size = 0.15, color = "black") +
                 # This magically only adds axes:
                 geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
                 geom_segment(data = segments,
                              aes(x = segments[, 1],
                                  y = segments[, 2],
                                  xend = segments[, 3],
                                  yend = segments[, 4]
                                  ),
                              size = 0.3,
                              colour = "red",
                              show.legend = FALSE
                              ) +
                 labs(x = NULL,
                      y = NULL,
                      title = title,
                      subtitle = subtitle
                      )
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

