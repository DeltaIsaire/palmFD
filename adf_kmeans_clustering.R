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
source(file = "functions/weighted_ADF_null_model_functions.R")

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

# We can also get full species richness from the original distribution data
palmdist <- read.csv(file = "data/palms_in_tdwg3.csv")
palm.richness <- count(palmdist$Area_code_L3)
names(palm.richness) <- c("tdwg3.code", "palm.richness")
palm.richness <- merge(palm.richness,
                       tdwg3.info[, "tdwg3.code", drop = FALSE],
                       all.y = TRUE,
                       by = "tdwg3.code"
                       )
# And the full presence / absence matrix:
pres.abs.all <- 
  table(palmdist) %>%
  as.data.frame.matrix() %>%
  as.matrix()


# Spatial map of tdwg3 units
# --------------------------
tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")
# Note I am using the 'sf' spatial format instead of 'SpatialPolygonsDataFrame'.
# The sf package with function sf::read_sf is more modern than using rgdal::readOGR.
# That modernity translates into improved performance.


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

# For all palm data:
communities <- rownames(pres.abs.all)
adf <- ADF(communities, pres.abs.all)
weights.all <-
  matrix(nrow = nrow(pres.abs.all),
         ncol = nrow(pres.abs.all),
         dimnames = list(focal = rownames(pres.abs.all),
                         adf = rownames(pres.abs.all)
                         ),
         data = 0  # countries not in the ADF have weight 0, not NA.
         )
for (i in seq_along(adf)) {
  # Get species in focal community
  focal.species <-
    pres.abs.all[match(names(adf)[i], rownames(pres.abs.all)),
                     ,
                    drop = FALSE
                    ] %>%
    { colnames(.)[. > 0] }
  # Get species in each ADF community
  community.species <-
    alply(adf[[i]],
          .margins = 1,
          function(x) {
            pres.abs.all[match(x, rownames(pres.abs.all)),
                             ,
                            drop = FALSE
                            ] %>%
              { colnames(.)[. > 0] }
          }
          )
    names(community.species) <- adf[[i]]
  # Define weights as number of shared species, and convert to probabilities
  community.weights <-
    llply(community.species, function(x) { sum(x %in% focal.species) } ) %>%
    unlist()
  community.chances <- community.weights / sum(community.weights)
  # Add this data to the matrix
  columns <- CrossCheck(x = colnames(weights.all),
                        y = names(community.chances),
                        presence = TRUE,
                        value = FALSE
                        )
  weights.all[i, columns] <- community.chances
}


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

# For all palms:
ssq.all <- numeric(length = length(2:100))
names(ssq.all) <- 2:100
for (i in seq_along(ssq.all)) {
  clust <- kmeans(weights.all,
                  centers = (i + 1),
                  nstart = 100,
                  iter.max = 100
                  )
  ssq.all[i] <- clust[["tot.withinss"]]
}
GraphSVG(plot(x = 2:100,
              y = ssq.all,
              type = "b",
              pch = 19,
              frame = FALSE,
              xlab = "Number of clusters K",
              ylab = "Total within-clusters sum of squares"
              ),
          file = paste0(plot.dir, "kmeans_k_ssq_all_palms.svg"),
          width = 8,
          height = 4
          )
# Looks like somewhere in the 20-30 area is optimal for this one


# Bayesian method
# ---------------
cat("Bayesian method...\n")
bayes.clust <- Mclust(adf.weight.matrix, G = 5:40)
bayes.clust$BIC
# First I picked max G of 20, but then 20 was the best fit. If the best fit
# is the biggest chosen k, then checking larger k is warranted.
# Based on BIC, the optimal k is 24 followed by 23 or 20.

# for all palms:
bayes.clust.all <- Mclust(weights.all, G = 5:40)
bayes.clust.all$BIC
# Based on BIC, the optimal k here is 12, followed by 13 or 11.
# But for k > 13 all values are NA, so we cannot be completely certain.


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


# For all palms:
result2 <- matrix(nrow = length(clust.indices),
                  ncol = 2,
                  dimnames = list(clust.indices, c("number_clusters", "Value_Index")),
                  data = NA
                  )
for (i in seq_along(clust.indices)) {
  cat(clust.indices[i], "\n")
  nb.clust <- NbClust(weights.all,
                      distance = "euclidean",
                      method = "kmeans",
                      min.nc = 2,
                      max.nc = 30,
                      index = clust.indices[i]
                      )
  result2[i, ] <- nb.clust$Best.nc
}
result2
# These results are also not conclusive.
# We can try the Dindex method:
GraphSVG(nb.clust <- NbClust(weights.all,
                             distance = "euclidean",
                             method = "kmeans",
                             min.nc = 2,
                             max.nc = 30,
                             index = "dindex"
                             ),
         file = paste0(plot.dir, "kmeans_k_dindex_all_palms.svg"),
         width = 8,
         height = 4
         )
# 25 scores pretty high here. 
# All in all, that is consistent with the Elbow graph, and many of the indices
# in result 2 being in the 20-30 range. 

# We should do two graphs, with k = 12 and k = 25.
# Actually, let's do 24 because the other one is also 24.


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


###################
# Plotting Clusters
###################
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

ggsave(plot = SpatialPlotFactor(tdwg.map, clusters$cluster, "ADF clusters"),
       filename = paste0(plot.dir, "adf_kmeans_clustering_24.png"),
       width = 8,
       height = 4
       )

clust.all.3 <-
  kmeans(weights.all, centers = 3, nstart = 100, iter.max = 100) %>%
  .[["cluster"]]
clust.all.3 %<>% .[-match("VNA", names(.))]
clusters <- data.frame(tdwg3.code = tdwg3.info$tdwg3.code,
                       cluster = rep(NA, nrow(tdwg3.info))
                       )
rows <- CrossCheck(x = clusters$tdwg3.code,
                   y = names(clust.all.3),
                   presence = TRUE,
                   value = FALSE
                   )
clusters[rows, "cluster"] <- clust.all.3

ggsave(plot = SpatialPlotFactor(tdwg.map, clusters$cluster, "ADF clusters"),
       filename = paste0(plot.dir, "adf_kmeans_clustering_all_palms_3.png"),
       width = 8,
       height = 4
       )

clust.all.5 <-
  kmeans(weights.all, centers = 5, nstart = 100, iter.max = 100) %>%
  .[["cluster"]]
clust.all.5 %<>% .[-match("VNA", names(.))]
clusters <- data.frame(tdwg3.code = tdwg3.info$tdwg3.code,
                       cluster = rep(NA, nrow(tdwg3.info))
                       )
rows <- CrossCheck(x = clusters$tdwg3.code,
                   y = names(clust.all.5),
                   presence = TRUE,
                   value = FALSE
                   )
clusters[rows, "cluster"] <- clust.all.5

ggsave(plot = SpatialPlotFactor(tdwg.map, clusters$cluster, "ADF clusters"),
       filename = paste0(plot.dir, "adf_kmeans_clustering_all_palms_5.png"),
       width = 8,
       height = 4
       )

clust.all.7 <-
  kmeans(weights.all, centers = 7, nstart = 100, iter.max = 100) %>%
  .[["cluster"]]
clust.all.7 %<>% .[-match("VNA", names(.))]
clusters <- data.frame(tdwg3.code = tdwg3.info$tdwg3.code,
                       cluster = rep(NA, nrow(tdwg3.info))
                       )
rows <- CrossCheck(x = clusters$tdwg3.code,
                   y = names(clust.all.7),
                   presence = TRUE,
                   value = FALSE
                   )
clusters[rows, "cluster"] <- clust.all.7

ggsave(plot = SpatialPlotFactor(tdwg.map, clusters$cluster, "ADF clusters"),
       filename = paste0(plot.dir, "adf_kmeans_clustering_all_palms_7.png"),
       width = 8,
       height = 4
       )


cat("Done.\n")

