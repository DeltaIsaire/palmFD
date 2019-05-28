## map plot example


## Prereqs
# --------
library(sf)
library(ggplot2)
library(scales)
library(viridis)
library(plyr)

theme_set(theme_bw())


# Data
# ----
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")
tdwg3.info[, "palm.richness"][tdwg3.info[, "palm.richness"] == 0] <- NA

tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")


# Basic map
# ---------
map <- ggplot(tdwg.map) + geom_sf()

# Map without Antarctica:
# -----------------------
index <- which(tdwg.map[["LEVEL_NAME"]] == "Antarctica")
sub.map <- tdwg.map[-index, ]
map2 <- ggplot(sub.map) + geom_sf()


# Map of species richness:
# ------------------------
SpatialPlotFill <- function(tdwg.map, vector, vector.name, title = NULL,
                            subtitle = NULL, colors = "viridis",
                            direction = 1, begin = 0) {
# tdwg.map: the spatial map, as an object of class 'sf'
# vector: vector with data to plot on the map. Must be of the same length as the
#         number of rows in tdwg.map (i.e. length 368 for the 368 tdwg3 units)
# vector.name: character string giving name for the vector (used in legend)
  tdwg.map$vector <- vector
  tdwg.plot <- ggplot(data = tdwg.map) + 
                 geom_sf(size = 0.15, color = "black", aes(fill = vector)) +
                 # This magically only adds axes:
                 geom_point(aes(x = "Long", y = "Lat"), size = 0, color = "white") +
                 labs(fill = vector.name,
                      x = NULL,
                      y = NULL,
                      title = title,
                      subtitle = subtitle
                      ) +
                 scale_fill_viridis(na.value = "#C0C0C0",
                                    discrete = is.discrete(vector),
                                    option = colors,
                                    direction = direction,
                                    begin = begin
                                    )
  tdwg.plot
}

rich.map <-
  SpatialPlotFill(sub.map,
                  vector = log(tdwg3.info[, "palm.richness"])[-index],
                  vector.name = "log(richness)",
                  title = "Palm species richness",
                  subtitle = "In botanical countries (TDWG3 units)"
                  )

