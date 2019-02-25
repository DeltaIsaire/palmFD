cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(sf)
library(ggplot2)
library(scales)
library(mclust)
library(NbClust)

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")

# Directory for saving plots (with a trailing slash):
plot.dir = "graphs/ADF_kmeans/"


###################
# Load spatial data
###################
cat("Preparing data...\n")

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
                    mapping = aes(x = Long, y = Lat, fill = factor),
                    size = 3,
                    shape = 21) +
         labs(fill = factor.name,
              x = NULL,
              y = NULL,
              title = title,
              subtitle = subtitle
              ) +
         scale_colour_discrete(drop = TRUE, limits = levels(subset$factor))
}


#########################################
# K-means clustering on ADF probabilities
#########################################
cat("Running K-means clustering on ADF probabilities:\n")
# Clustering TDWG3 units by 'community similarity', as defined by proportion
# of shared species.

# Load adf weights matrix
# -----------------------
# Rows are focal TDWG3 units, columns are ADF TDWG3 units.
# You can tell because the rowSums will add to 1, but the colSums won't.
# This is also how functions expect the data to be formatted.
adf.weight.matrix <- read.csv(file = "output/adf_weight_matrix.csv",
                              header = TRUE,
                              row.names = 1,
                              check.names = FALSE
                              )

# Argument 'centers' in kmeans() is k, i.e. how many clusters to identify.
# The real trick is choosing the right value for k.
# The following is based on 
# https://www.r-bloggers.com/finding-optimal-number-of-clusters/

set.seed(999)

if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
}

# 'Elbow method'
# --------------
cat("Finding optimal k: \nelbow method...\n")
# Run clustering for differing k, and compare the fit using
# total within-clusters sum of squares.
# with 131 TDWG3 units to cluster, 100 clusters is clearly overkill,
# so that's a good upper bound.
ssq <- numeric(length = length(2:100))
names(ssq) <- 2:100
for (i in seq_along(ssq)) {
  clust <- kmeans(adf.weight.matrix,
                  centers = (i + 1),
                  nstart = 100,
                  iter.max = 100
                  )
  ssq[i] <- clust[["tot.withinss"]]
}
GraphSVG(plot(x = 2:100,
              y = ssq,
              type = "b",
              pch = 19,
              frame = FALSE,
              xlab = "Number of clusters K",
              ylab = "Total within-clusters sum of squares"
              ),
          file = paste0(plot.dir, "kmeans_k_ssq.svg"),
          width = 8,
          height = 4
          )
# That's a neat saturation curve. The problem is, defining the saturation point
# is more or less arbitrary. At best, this graph narrows the search space.
# The optimal k is likely within 5 - 20.

# Bayesian method
# ---------------
cat("Bayesian method...\n")
bayes.clust <- Mclust(adf.weight.matrix, G = 5:40)
bayes.clust$BIC
# First I picked max G of 20, but then 20 was the best fit. If the best fit
# is the biggest chosen k, then checking larger k is warranted.
# Based on BIC, the optimal k is 24 followed by 23 or 20.

# NbClust method
# --------------
cat("Using methods from package 'Nbclust'...\n")
# Method using the NbClust package.
# Using only the clustering methods that work for our data, i.e.
# do not result in an error.
clust.indices <- c("kl", "ch", "hartigan",
                   "cindex", "db", "silhouette",
                   "duda", "pseudot2", "beale", "ratkowsky", "ball",
                   "ptbiserial", "frey", "mcclain", "dunn", "sdindex", "sdbw"
                   )
result <- matrix(nrow = length(clust.indices),
                 ncol = 2,
                 dimnames = list(clust.indices, c("number_clusters", "Value_Index")),
                 data = NA
                 )
for (i in seq_along(clust.indices)) {
  cat(clust.indices[i], "\n")
  nb.clust <- NbClust(adf.weight.matrix,
                      distance = "euclidean",
                      method = "kmeans",
                      min.nc = 2,
                      max.nc = 30,
                      index = clust.indices[i]
                      )
  result[i, ] <- nb.clust$Best.nc
}
result
# These results are not conclusive.
# We can try the Dindex method:
GraphSVG(nb.clust <- NbClust(adf.weight.matrix,
                             distance = "euclidean",
                             method = "kmeans",
                             min.nc = 2,
                             max.nc = 30,
                             index = "dindex"
                             ),
         file = paste0(plot.dir, "kmeans_k_dindex.svg"),
         width = 8,
         height = 4
         )
# No clear winners, but strongest peaks at k = 3 and k = ~12.
# K = 3 makes sense if it corresponds to realms, but is obviously too broad.

# These NbClust methods are hard to interpret without really getting into it.
# I have more faith in the Bayesian clustering and the elbow graph.
#
# Based on that, I consider k = 24 to be a reasonable number.

set.seed(999)

adf.clusters <-
  kmeans(adf.weight.matrix, centers = 24, nstart = 100, iter.max = 100) %>%
  .[["cluster"]]
clusters <- data.frame(tdwg3.code = tdwg3.info$tdwg3.code,
                       cluster = rep(NA, nrow(tdwg3.info))
                       )
rows <- CrossCheck(x = clusters$tdwg3.code,
                   y = names(adf.clusters),
                   presence = TRUE,
                   value = FALSE
                   )
clusters[rows, "cluster"] <- adf.clusters


###################
# Plotting Clusters
###################

ggsave(plot = SpatialPlotFactor(tdwg.map, clusters$cluster, "ADF clusters"),
       filename = paste0(plot.dir, "adf_kmeans_clustering_24.png"),
       width = 8,
       height = 4
       )


cat("Done.\n")

