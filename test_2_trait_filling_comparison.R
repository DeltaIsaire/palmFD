###################################################
# Palm FD project: Trait Filling Methods Comparison
###################################################
#
# In which the results of different trait gap-filling trials are compared.
# These trials were generated in 'test_1_trait_filling.R'.
#
# Input files:
#   output/palm.traits.csv
#   output/traits_filled_genus_mean.csv
#   output/traits_filled_BHPMF_one.csv
#   output/traits_filled_BHPMF_two.csv
#   output/traits_filled_BHPMF_three.csv
#   output/traits_filled_BHPMF_four.csv
#   output/traits_filled_BHPMF_five.csv
# Generated output files:
#   graphs/test_estimates_species_counts.svg
#   graphs/ <several graphs comparing trait estimates.
#            Filenames all start with 'test_estimates_'>


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


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
#   BHPMF, in different ways
filled.BHPMF.one <- read.csv(file = "output/traits_filled_BHPMF_one.csv")
filled.BHPMF.two <- read.csv(file = "output/traits_filled_BHPMF_two.csv")
filled.BHPMF.three <- read.csv(file = "output/traits_filled_BHPMF_three.csv")
filled.BHPMF.four <- read.csv(file = "output/traits_filled_BHPMF_four.csv")
filled.BHPMF.five <- read.csv(file = "output/traits_filled_BHPMF_five.csv")
all.filled <- list(original = palm.traits,
                   mean = filled.mean,
                   b.one = filled.BHPMF.one,
                   b.two = filled.BHPMF.two,
                   b.three = filled.BHPMF.three,
                   b.four = filled.BHPMF.four,
                   b.five = filled.BHPMF.five
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
all.estimates <- llply(as.list(trait.names),
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


# -----------------------------------
# Comparing all gap-filling scenarios
# -----------------------------------
cat("Comparing all gap-filling scenarios:\n")
cat("(1) Creating barplot of species counts...\n")
# Barplot of the number of species retained in each method.
GraphSVG(BarPlot(species.counts[, 2],
                 names = species.counts[, 1],
                 ylab = "Species count",
                 ),
         file = "graphs/test_estimates_species_counts.svg",
         width = 6,
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
           file = paste0("graphs/test_estimates_distributions_",
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
# For each trait, compare standard BHPMF (BHPMF.one) to all other estimates.
y.names <-
  names(all.estimates[[1]]) %>%
  .[!. %in% c("species", "original", "b.one")]
for (df in seq_along(all.estimates)) {
  for (y in y.names) {
    GraphSVG(MultiScatter(all.estimates[[df]][, "b.one"],
                          all.estimates[[df]][, y],
                          x.name = paste("Est.",
                                         trait.names[df],
                                         "(BHPMF.one)"
                                         ),
                          y.name = paste0("Est. ",
                                          trait.names[df],
                                          " (",
                                          y,
                                          ")"
                                          )
                          ),
             file=paste0("graphs/test_estimates_scatter_",
                         trait.names[df],
                         "_b.one_vs_",
                         y,
                         ".svg"
                         ),
             width = 12,
             height = 4
             )
  }
}

# BHPMF: estimates vs observed
# ----------------------------
cat("(4) Plotting observed traits vs BHPMF estimates...\n")
# BHPMF also output estimates for 'observed' trait values.
# How accurate are these?
# For each trait, plot 'observed' values vs estimates from filled.BHPMF.two.
# First prepare the data:
observed.estimates <- 
  llply(as.list(trait.names),
        function(trait) {
          observed <- palm.traits[!is.na(palm.traits[, trait]),
                                  c("species", trait)
                                  ]
          estimate <- 
            CrossCheck(filled.BHPMF.two$species,
                       observed$species,
                       value = FALSE
                       ) %>%
            filled.BHPMF.two[., trait]
          data.frame(observed, estimate = estimate)
        }
        )
names(observed.estimates) <- trait.names
# Then make Scatterplots:
for (i in seq_along(observed.estimates)) {
  df <- observed.estimates[[i]]
  GraphSVG(MultiScatter(x = df[, 2],
                        y = df[, 3],
                        x.name = paste("Observed",
                                       names(df)[2]
                                       ),
                        y.name = paste("Estimated",
                                       names(df)[2]
                                       )
                        ),
           file = paste0("graphs/test_estimates_original_vs_estimated_",
                         names(df)[2],
                         ".svg"
                         ),
           width = 12,
           height = 4
           )
}

# BHPMF: Effect of adding dummy traits
# ------------------------------------
cat("(5) Assessing effect of adding dummy traits...\n")
# in BHPMF.three we added one dummy trait.
# in BHPMF.five we added ten dummy traits.
# Including a dummy allows inclusion in BHPMF of species with NA for all
# real traits. But does the dummy variable distort estimates?
# Scatterplots made above already compare b.three and b.five to b.one.
# Here, we first compare b.three to b.five, for each trait:
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "b.three"],
                        all.estimates[[df]][, "b.five"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "(BHPMF.three)"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "BHPMF.five"
                                       )
                        ),
             file=paste0("graphs/test_estimates_scatter_",
                         trait.names[df],
                         "_b.three_vs_b.five.svg"
                         ),
             width = 12,
             height = 4
             )
}

cat("Done.\n")

