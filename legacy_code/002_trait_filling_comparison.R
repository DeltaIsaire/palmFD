###################################################
# Palm FD project: Trait Filling Methods Comparison
###################################################
#
# In which the results of genus-mean and BHPMF gap-filling are compared.
#
# Input files:
#   output/palm_traits.csv
#   output/traits_filled_genus_mean.csv
#   output/traits_filled_BHPMF.csv
#
# Generated output files:
#   graphs/gapfilling_estimates_species_counts.svg
#   graphs/ < several graphs comparing trait estimates.
#             Filenames all start with 'gapfilling_estimates_' >


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")

if (!dir.exists("graphs/")) {
  cat("creating directory: graphs/\n")
  dir.create("graphs/")
}


# ---------------------------------------
# Data preparation: filled trait matrices
# ---------------------------------------
cat("Preparing data...\n")

# Trait filling begins with the unimputed trait matrix, which has for each
# species a single 'observed' value for each trait, with gaps (NAs).
palm.traits <- read.csv(file = "output/palm_traits.csv")
# Subsequently we have filled this trait matrix in the following ways:
#   genus means
filled.mean <- read.csv(file = "output/traits_filled_genus_mean.csv")
#   BHPMF
filled.BHPMF <- read.csv(file = "output/traits_filled_BHPMF.csv")

all.filled <- list(original = palm.traits,
                   mean = filled.mean,
                   BHPMF = filled.BHPMF
                   )
filled.names <- names(all.filled)

# First, extract the number of species retained in each method:
species.counts <- ldply(all.filled,
                        function (x) { length(x[, 1]) }
                        )
names(species.counts) <- c("data", "species.count")

# Second, the different methods had to exclude different species, so subset
# the trait matrices to the list of shared species.
species.list <- llply(all.filled, function(x) { x[, "species"] } )
shared.species <- 
  MultiCheck(species.list) %>%
  as.character()
all.filled <- lapply(all.filled,
                     function(x) {
                       indices <- CrossCheck(x = x[, "species"],
                                             y = shared.species,
                                             value = FALSE
                                             )
                       x[indices, ]
                     }
                     )
# Third, we need for each trait the subset of values that were estimated
# (filled in), because that's most interesting to compare.
trait.names <- c("stem.height", "blade.length", "fruit.length")
all.estimates <- 
  llply(as.list(trait.names),
        function(trait) { 
          missing <- which(is.na(palm.traits[, trait]))
          species <- CrossCheck(x = palm.traits[missing, "species"],
                                y = shared.species
                                )
          estimates <- llply(all.filled,
                             function(df) {
                               indices <- CrossCheck(x = df$species,
                                                     y = species,
                                                     value = FALSE
                                                     )
                               df[indices, trait]
                             }
                             )
          data.frame(species = species,
                     estimates %>%
                       simplify2array() %>%
                       as.data.frame()
                     )
        }
        )
names(all.estimates) <- paste0(trait.names, ".estimates")


# ---------------------------------
# Comparing the gap-filling methods
# ---------------------------------
cat("Comparing the gap-filling methods:\n")
cat("(1) Creating barplot of species counts...\n")
# Barplot of the number of species retained in each method.
GraphSVG(BarPlot(species.counts[, 2],
                 names = species.counts[, 1],
                 ylab = "Species count",
                 ),
         file = "graphs/gapfilling_estimates_species_counts.svg",
         width = 8,
         height = 4
         )

# Histograms of estimated trait values
# ------------------------------------
cat("(2) Creating histograms of estimated trait values...\n")
# For each trait, a stacked histogram of the trait value estimates for each method
ids <- names(all.estimates[[1]])[-c(1, 2)]
for (i in seq_along(all.estimates)) {
  GraphSVG(MultiHist(all.estimates[[i]][, ids],
                     id = ids,
                     xlab = paste("Estimated", trait.names[i])
                     ),
           file = paste0("graphs/gapfilling_estimates_distributions_",
                         trait.names[i],
                         ".svg"
                         ),
           width = 6,
           height = 1 * max(seq_along(all.estimates[[i]]))
           )
}

# Scatterplots of estimated trait values
# --------------------------------------
cat("(3) Creating scatterplots comparing trait estimates...\n")
# For each trait, compare genus means to standard BHPMF.
y.names <-
  names(all.estimates[[1]]) %>%
  .[!. %in% c("species", "original", "mean")]
for (df in seq_along(all.estimates)) {
  for (y in y.names) {
    GraphSVG(MultiScatter(all.estimates[[df]][, "mean"],
                          all.estimates[[df]][, y],
                          x.name = paste("Est.",
                                         trait.names[df],
                                         "(genus-mean)"
                                         ),
                          y.name = paste0("Est. ",
                                          trait.names[df],
                                          " (",
                                          y,
                                          ")"
                                          )
                          ),
             file=paste0("graphs/gapfilling_estimates_scatter_",
                         trait.names[df],
                         "_genus_mean_vs_",
                         y,
                         ".svg"
                         ),
             width = 12,
             height = 4
             )
  }
}

cat("Done.\n")

