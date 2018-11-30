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
#   graphs/test_estimates_mean_BHPMF_height.svg
#   graphs/test_estimates_mean_BHPMF_blade.svg
#   graphs/test_estimates_mean_BHPMF_fruit.svg
#   graphs/test_estimates_completeness.svg


library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


# ---------------------------------------
# Data preparation: filled trait matrices
# ---------------------------------------
cat("Comparing trait-filling methods:\n")
cat("Preparing data...\n")
# Trait filling begins with the unimputed trait matrix, which has for each
# species a single 'observed' value for each trait, with gaps (NAs).
palm.traits <- read.csv(file = "output/palm.traits.csv")
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
lengths <- ldply(all.filled,
                 function (x) { length(x[, 1]) }
                 )
names(lengths) <- c("data", "species.count")

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


# Comparing genus means with standard BHPMF (trial 1):
# Scatterplots of gap-filled trait values
# ----------------------------------------------------
cat("(1) Comparing genus-mean with BHPMF:\n")
cat("Creating scatterplots...\n")
# Scatterplots of BHPMF trait estimates ~ genus mean estimates.
Temp <- function(data, title) {
# data: dataframe with the estimates for each trait
# title: title for the overall graph. Should include the trait name.
  MultiScatter(x = data[, "genus.mean"],
               y = data[, "BHPMF"],
               x.name = "Genus Mean",
               y.name = "BHPMF",
               title = title
               )
}
GraphSVG(Temp(height.estimates,
              title = "Estimated stem height (m)"
              ),
         file = "graphs/test_estimates_mean_BHPMF_height.svg",
         width = 12,
         height = 4
         )
GraphSVG(Temp(blade.estimates,
              title = "Estimated blade length (m)"
              ),
         file = "graphs/test_estimates_mean_BHPMF_blade.svg",
         width = 12,
         height = 4
         )
GraphSVG(Temp(fruit.estimates,
              title = "Estimated fruit length (cm)"
              ),
         file = "graphs/test_estimates_mean_BHPMF_fruit.svg",
         width = 12,
         height = 4
         )
cat("Done.\n")


# -----------------------------------
# Comparing all gap-filling scenarios
# -----------------------------------
cat("(2) Comparing all gap-filling scenarios:\n")
cat("Preparing data...\n")
# Trait filling begins with the unimputed trait matrix, which has for each
# species a single 'observed' value for each trait, with gaps (NAs).
# Subsequently we have filled this trait matrix in the following ways:
#   genus means
filled.mean <- read.csv(file = "output/traits_filled_genus_mean.csv")
#   BHPMF, in different ways
filled.BHPMF.one <- read.csv(file = "output/traits_filled_BHPMF_one.csv")
filled.BHPMF.two <- read.csv(file = "output/traits_filled_BHPMF_two.csv")
filled.BHPMF.three <- read.csv(file = "output/traits_filled_BHPMF_three.csv")
filled.BHPMF.four <- read.csv(file = "output/traits_filled_BHPMF_four.csv")
filled.BHPMF.five <- read.csv(file = "output/traits_filled_BHPMF_five.csv")
all.filled <- list(palm.traits,
                   filled.mean,
                   filled.BHPMF.one,
                   filled.BHPMF.two,
                   filled.BHPMF.three,
                   filled.BHPMF.four,
                   filled.BHPMF.five
                   )
filled.names <- c("original", "mean", "one", "two", "three", "four", "five")

# First, extract the number of species retained in each method:
lengths <- ldply(all.filled,
                 function (x) { length(x[, 1]) }
                 )
rownames(lengths) <- filled..names
names(lengths) <- "species.count"
# Next, we need for each trait the subset of values that were estimated
# (filled in), because that's what we should compare.


# Completeness
# ------------
cat("Creating barplot of completeness...\n")
# Barplot of the number of species retained in each method.
GraphSVG(barplot(lengths[, 1],
                 names.arg = rownames(lengths),
                 ylab = "species count"
                 ),
         file = "graphs/test_estimates_completeness.svg",
         width = 6,
         height = 4
         )

# The different methods had to exclude different species, so first we must
# subset the filled trait matrices to the list of shared species.
