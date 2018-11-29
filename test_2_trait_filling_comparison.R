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


# ---------------------------------------------------
# Comparing genus means with standard BHPMF (trial 1)
# ---------------------------------------------------
cat("Comparing genus-mean with BHPMF:", "\n")
cat("Preparing data...", "\n")
palm.traits <- read.csv(file = "output/palm.traits.csv")
filled.mean <- read.csv(file = "output/traits_filled_genus_mean.csv")
filled.BHPMF <- read.csv(file = "output/traits_filled_BHPMF_one.csv")

# The different methods had to exclude different species, so first we must
# subset the filled trait matrices to the list of shared species.
shared.species <- CrossCheck(filled.mean$species,
                             filled.BHPMF$species,
                             presence = TRUE,
                             value = TRUE
                             )
filled.mean %<>% .[CrossCheck(.$species,
                              shared.species,
                              presence = TRUE,
                              value = FALSE
                              ),
                    ]
filled.BHPMF %<>% .[CrossCheck(.$species,
                               shared.species,
                               presence = TRUE,
                               value = FALSE
                               ),
                     ]
# Next, we need for each trait the subset of values that were estimated
# (filled in), because that's what we should compare.
Temp <- function(trait) {
# trait: character vector of length 1 giving the trait column name
  missing <- which(is.na(palm.traits[, trait]))
  species <- CrossCheck(palm.traits[missing, "species"],
                        shared.species,
                        presence = TRUE,
                        value = TRUE
                        )
  genus.mean <- CrossCheck(filled.mean$species,
                           species,
                           presence = TRUE,
                           value = FALSE
                           )
  BHPMF <- CrossCheck(filled.BHPMF$species,
                      species,
                      presence = TRUE,
                      value = FALSE
                      )
  data.frame(species = species,
             genus.mean = filled.mean[, trait][genus.mean],
             BHPMF = filled.BHPMF[, trait][BHPMF]
             )
}
height.estimates <- Temp("stem.height")
blade.estimates <- Temp("blade.length")
fruit.estimates <- Temp("fruit.length")


# Scatterplots of gap-filled trait values
# ---------------------------------------
cat("Creating scatterplots...", "\n")
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

cat("Done.", "\n")


# -----------------------------------
# Comparing all gap-filling scenarios
# -----------------------------------
cat("Comparing all gap-filling scenarios:", "\n")
cat("Preparing data...", "\n")
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

# The different methods had to exclude different species, so first we must
# subset the filled trait matrices to the list of shared species.
all.filled <- list(palm.traits,
                   filled.mean,
                   filled.BHPMF.one,
                   filled.BHPMF.two,
                   filled.BHPMF.three,
                   filled.BHPMF.four,
                   filled.BHPMF.five
                   )
all.names <- c("original", "mean", "one", "two", "three", "four", "five")

# Completeness
# ------------
# Barplot of the number of species retained in each method.
lengths <- ldply(all.filled,
                 function (x) { length(x[, 1]) }
                 )
rownames(lengths) <- all.names
names(lengths) <- "species.count"
GraphSVG(barplot(lengths[, 1],
                 names.arg = rownames(lengths),
                 ylab = "species count"
                 ),
         file = "graphs/test_estimates_completeness.svg",
         width = 6,
         height = 4
         )


