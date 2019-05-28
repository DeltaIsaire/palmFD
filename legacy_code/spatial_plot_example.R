# R code example: plotting a worldmap of TDWG3 units, with points to indicate some
# sort of data.
# Written with R version 3.6.0, and I think you need at least R version 3.5.2.



# Required libraries:
library(sf)  # for handling the shapefiles
library(ggplot2)  # for generating the plots

# Optional: set a neat and simple black-on-white theme for the spatial plots
theme_set(theme_bw())


# ---------------------
# Loading the shapefile
# ---------------------
# Read the shapefile using sf::read_sf()
# Argument 'dsn' is the folder where the shapefiles are stored
# Argument 'layer' is the name of the shapefile, without filetype extension
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")
# The resulting object is of a complex type, but for the most part you can
# pretend it is a data frame.


# ------------------------------
# Basics of making spatial plots
# ------------------------------
# Plotting the polygon data (the 'shapes') is fairly simple:
p <- ggplot(tdwg.map) +
     geom_sf()
p  # to render the plot
# Rendering spatial plots is computionally expensive...
# This is just a featureless plot. Now we can expand on it in various ways, using
# the plotting methods of the ggplot package. Learning to use ggplot is a bit
# tricky. You should search for a good tutorial.

# Improve the visuals of the polygons:
p <- ggplot(tdwg.map) +
     geom_sf(size = 0.15, color = "white", fill = "#B0B0B0")
# (that last value is a hexadecimal rgb color value, in this case light grey)

# If necessary, you can add axes indicating latitude and longitude using this
# workaround:
p <- ggplot(tdwg.map) +
     geom_sf(size = 0.15, color = "white", fill = "#B0B0B0") +
     geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white")
# This works by plotting fake invisible points.
# Normally, variable names inside the aes() function should NOT be quoted.


# ------------------------------
# Adding data points to the plot
# ------------------------------
# To visualize data of the TDWG3 units on the map, we can add points to the
# centroids of the TDWG3 units with geom_point. This involves some clever
# data management, so it takes a few steps.

# (1) First we need some data.
# ---------------------------
# Let's make a simple dataframe with mock data, based on the spatial data object
cat <- ifelse(tdwg.map$Long < -20,
              "West",
              ifelse(tdwg.map$Long > 60,
                     "East",
                     "Central"
                     )
              )
my.data <- data.frame(area = tdwg.map$LEVEL_3_CO,
                      num  = sqrt(abs(tdwg.map$Lat)),
                      cat  = as.factor(cat)
                      )
# Make sure categorical data is a factor.

# Palms do not occur in every TDWG3 unit, so some of your real data will be NA.
# Here, let's pretend every third data point is missing.
my.data[c(FALSE, FALSE, TRUE), -1] <- NA

# (2) Merge your data with the spatial object
# -------------------------------------------
tdwg.data <- merge(x = tdwg.map,
                   y = my.data,
                   by.x = "LEVEL_3_CO",
                   by.y = "area"
                   )
# Do not overwrite the original spatial object. We will need it later

# (3) Subset the merged data to complete cases
# --------------------------------------------
# In other words, remove the TDWG3 units for which you do not want to plot points
subset <- tdwg.data[!is.na(tdwg.data$num), ]

# (4) Create the spatial map with our data points
# -----------------------------------------------
# We expand the basic plotting code by adding a geom_point, where the coordinates
# of the points are given by the latitude and longitude data.
# The color of the points is used to indicate the value of our chosen data vector.
# Note we use the original spatial object 'tdwg.map' for plotting the shapes,
# and the 'subset' object for plotting our points
p <- ggplot(tdwg.map) +
     geom_sf(size = 0.15, color = "white", fill = "#B0B0B0") +
     geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
     geom_point(data = subset,
                size = 3,  # Size of the points. Play around with this
                shape = 21,  # Black circle to be filled with color specified below
                aes(x = Long,  # Do not use quotes for the variable name
                    y = Lat,  # Do not use quotes for the variable name
                    fill = cat  # define fill colors by the 'cat' variable
                    )
                )


# ------------------------------
# Saving a plot as an image file
# ------------------------------
# Saving a ggplot plot is done with ggsave(). It works like this:
ggsave(plot = p,  # the plot object created above
       filename = "graphs/spatial_plot_example.png",
       width = 8,
       height = 4,
       dpi = 300
       )
# The filename should include a file extension. I recommend .png.
# Width and height specify the dimensions of the graph in inches, while 'dpi'
# specifies the quality as dots (pixels) per inch. So here, the resulting image
# will be 8*300 by 4*300 = 2400 by 1200 pixels.


# ----------------------------
# Advanced plotting techniques
# ----------------------------
# You now know the basics of plotting spatial data. But there are a number of
# tweaks that can be made to enhance the plot: adding titles, choosing your
# own colors, editing the legend, and so on. Read a ggplot tutorial to find out
# how most of it works.

# Here I will show one of my custom plotting functions, as an example of what is
# possible. This function does the following, for a given numerical vector
# and spatial object:
# - merge the vector into the spatial object
# - transform the vector into quantile bins, which are plotted instead of the
#   raw data, to improve contrast
# - the raw data is instead used for the size of the points
# - remove missing values
# - tweak the displayed legend
# - generate a colourblind-friendly color scale
# - Parse titles for the plots, axes and legend
# - create the spatial plot
# The advantage of a custom function like this is that it saves you a lot of work
# if you want to create more than one plot of the same type.

# This function needs additional libraries:
library(viridis)
library(scales)

# Function definition, with added annotation
# ------------------------------------------
SpatialPlotBins <- function(tdwg.map, vector, vector.name, vector.size = NULL,
                            title = NULL, subtitle = NULL,
                            legend.position = "right", labels = "default",
                            colors = "viridis") {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map. Data must be a continuous numeric variable.
#         For best results the order of values should match the order of areas in
#         'tdwg.map'.
# vector.size: vector of the same length as 'vector' giving size of points. If NULL,
#              size will be scaled by the values of 'vector'
# vector.name: character string giving name for the vector (used in legend)
# labels: Legend labels. By default, labels will be names of bins.
#         Can be a character vector giving legend labels.
#         Can be a function that takes the vector quantiles as input and returns
#         labels as output.
# colors: character string giving the viridis colormap to use, see argument 'option'
#         to viridis::scale_color_viridis().

  # Add the vector to the spatial object
  tdwg.map$vector <- vector

  # Compute quantile bins and add those to the spatial object.
  # Note we are making 1 / 0.125 = 8 bins.
  vector.quantiles <-
    quantile(vector, probs = seq(0, 1, 0.125), na.rm = TRUE)
  vector.bins <- cut(vector, vector.quantiles, include.lowest = TRUE)
  tdwg.map$bins <- vector.bins

  # Subset the data (removing NAs) and define the plotting point size
  if (!is.null(vector.size)) {
    tdwg.map$vector.size <- vector.size
    subset <- tdwg.map[!is.na(tdwg.map$vector), ]
    size <- rescale(subset$vector.size, to = c(0.25, 1)) * 3.5
  } else {
    subset <- tdwg.map[!is.na(tdwg.map$vector), ]
    size <- rescale(subset$vector, to = c(0.25, 1)) * 3.5
  }

  # Define size of quantile bin points in the legend
  legend.size <- rescale(vector.quantiles[-1], to = c(0.25, 1)) * 3.5 * 1.3

  # Define labels for the quantile bins in the legend.
  # The default is a character string giving the value range of the quantile bin.
  if (identical(labels, "default")) {
    x <- vector.quantiles
    lower <- round(x[1:(length(x) - 1)], digits = 2)
    upper <- round(x[2:length(x)], digits = 2)
    legend.labels <- paste(lower, "-", upper)
  } else {
    if (is.function(labels)) {
      legend.labels <- do.call(labels, list(vector.quantiles))
    } else {
      legend.labels <- labels
    }
  }

  # Define the colors for the points in the plot.
  # Note that this is also where the legend labels are set.
  color.scale <- scale_fill_viridis(discrete = TRUE,
                                    option = colors,
                                    labels = legend.labels,
                                    breaks = levels(vector.bins)
                                    )

  # Finally, create the plot
  ggplot(data = tdwg.map) +
         # plot the shapes
         geom_sf(size = 0.15, color = "white", fill = "#B0B0B0") +
         # This magically only adds axes:
         geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
         # While this does NOT add axes but does add points:
         geom_point(data = subset,
                    mapping = aes(x = Long, y = Lat, fill = bins),
                    size = size,
                    shape = 21
                    ) +
         # Specify all plot labels, including titles. NULL means no label.
         labs(fill = vector.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              ) +
         # Add the defined color scale. ggplot is smart enough to know where
         # these colors should be used
         color.scale +
         # Tweak the legend display
         theme(legend.position = legend.position,
               legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
               legend.box.spacing = unit(0.1, "cm")
               ) +
         guides(fill = guide_legend(override.aes = list(size = legend.size)))

  # By default, functions return the last evaluated expression, in this case
  # returning the ggplot object
}

# Applying the plotting function
# ------------------------------
# Remember we have the spatial object:
tdwg.map
# and our mock dataset:
my.data

# Now we can simply plot the 'num' variable like this:
p <- SpatialPlotBins(tdwg.map,
                     vector = my.data$num,
                     vector.name = "sqrt(lat)",
                     title = "Square-root of Latitude"
                     )

# And for the last example, remove Antarctica from the plot because it takes up
# a lot of space and isn't very interesting.
# First we make a logical vector specifying our subset: which TDWG3 units are
# NOT antarctica.
no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"
# Then make the plot as follows. Let's put the legend at the bottom where
# Antarctica used to be.
p <- SpatialPlotBins(tdwg.map[no.ANT, ],
                     vector = my.data$num[no.ANT],
                     vector.name = "sqrt(lat)",
                     title = "Square-root of Latitude",
                     legend.position = "bottom"
                     )

