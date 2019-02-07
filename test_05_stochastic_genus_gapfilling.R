###############################################################
# Palm FD Project: Stochastic gapfilling of trait data by Genus
###############################################################

# Previous gapfilling of the trait matrix has used Genus means, where gaps are 
# filled with a single genus mean trait value. This does not account for the 
# uncertainty in the Genus mean estimate.
# Here, we do gapfilling by stochastic sampling of the Genus-level trait
# distributions. For each Genus, find the mean and sd of each trait, then 
# fill trait gaps with a random sample from a normal distribution with that mean
# and sd. The procedure is repeated 100 times, resulting in 100 versions of the
# Genus-gapfilled dataset.
#
# This script supersedes all previous trait-filling scripts.

# Input files:
#   data/palms_in_tdwg3.csv
#   data/PalmTraits_10.csv
# Generated output files:
#   output/test/palm_traits.csv
#   output/test/stochastic_gapfilling_excluded_genera.csv
#   output/test/stochastic_gapfilling_genus_distributions_stem.height.csv
#   output/test/stochastic_gapfilling_genus_distributions_blade.length.csv
#   output/test/stochastic_gapfilling_genus_distributions_fruit.length.csv
#   dir 'output/test/stochastic_gapfilled/' with stochastically filled
#     trait matrices
#   graphs/test/genus_distribution_completeness.svg
#   graphs/test/genus_distribution_height.png
#   graphs/test/genus_distribution_blade.png
#   graphs/test/genus_distribution_fruit.png

cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


# Set random seed for reproducible output.
# Using 9-digit integer generated with runif(1, max = 1e9)
set.seed(336421650)


####################
# Prepare trait data
####################
cat("Preparing trait data...\n")
palm.dist <- read.csv(file = "data/palms_in_tdwg3.csv")
# Long-list Data frame with two columns:
#	1. Area_code_L3 - 3-letter code for each botanical country
#	2. SpecName - Species name
trait.data <- 
  read.csv(file = "data/PalmTraits_10.csv") %>%
  .[order(.$SpecName), ]
trait.data$SpecName %<>% sub(pattern = " ",
                             replacement = "_",
                             x = .
                             )
# Data frame with 31 columns:
#	1. SpecName - binomial species name [sorted alphabetically]
#	2. accGenus - only the genus name
#	3. accSpecies - only the species name
#	4. PalmTribe
#	5. PalmSubfamily
#	6-30 - information on 25 traits
#	31. extra reference notes

# Do the datasets agree on which palm species exist?
species.list <- list(palm.dist  = levels(palm.dist$SpecName),
                     trait.data = trait.data$SpecName 
                     )
species.shared <-
  MultiCheck(species.list, presence=TRUE, value=TRUE, unique = TRUE) %>%
  sort()

# Create palm traits dataframe and subset to species we have occurrence data
# for.
# Note: there should be no missing species, so the subset is for posterity.
palm.traits <- 
  data.frame(species      = trait.data$SpecName,
             genus        = trait.data$accGenus,
             stem.height  = trait.data$MaxStemHeight_m,
             blade.length = trait.data$Max_Blade_Length_m,
             fruit.length = trait.data$AverageFruitLength_cm
             ) %>%
  .[CrossCheck(x = .$species,
               y = species.shared,
               presence = TRUE,
               value = FALSE,
               unique = TRUE
               ),
    ]
# stem.height = maximum stem height in meters.
# blade.length = maximum blade length in meters.
# fruit.length = average fruit length in centimeters.
# Genus is included because we need it to calculate genus trait distributions
write.csv(palm.traits, 
          file = "output/test/palm_traits.csv",
          eol = "\r\n",
          row.names = FALSE
          )
# This will be useful:
trait.names <- c("stem.height", "blade.length", "fruit.length")


#######################################################
# Calculate observed trait distributions for each Genus
#######################################################
cat("Determining observed Genus-level trait distributions...\n")

# First, subset to Genera with at least 2 observed values for each trait.
# -----------------------------------------------------------------------
# Find genera to exclude:
genus.counts <- ddply(palm.traits, "genus", numcolwise(CountObserved))
exclude.indices <-
  adply(genus.counts[, -1], .margins = 1, function(x) { any(x < 2) } ) %>%
  { which(.[, "V1"]) }
exclude.genera <- unique(palm.traits[, "genus"]) [exclude.indices]

# Verify visually and save excluded data:
excluded <- palm.traits[palm.traits[, "genus"] %in% exclude.genera, ]
# 78 species in 63 Genera are excluded.
# These are mainly genera with only 1 species.
write.csv(excluded, 
          file = "output/test/stochastic_gapfilling_excluded_genera.csv",
          eol = "\r\n",
          row.names = FALSE
          )

# Subset palm traits dataset
palm.traits.subset <- palm.traits[!palm.traits[, "genus"] %in% exclude.genera, ]


# Find the mean, sd, min and max of each trait.
# ---------------------------------------------
genus.means <- ddply(palm.traits.subset, "genus", numcolwise(mean, na.rm = TRUE))
genus.sd <- ddply(palm.traits.subset, "genus", numcolwise(sd, na.rm = TRUE))
genus.min <- ddply(palm.traits.subset, "genus", numcolwise(min, na.rm = TRUE))
genus.max <- ddply(palm.traits.subset, "genus", numcolwise(max, na.rm = TRUE))
# let's include total count and missing counts as well, for completeness
genus.species <- ddply(palm.traits.subset, "genus", nrow)
genus.observed <- ddply(palm.traits.subset, "genus", numcolwise(CountObserved))

# Combine distribution data into a dataframe for each trait
genus.distributions <- vector("list", length = length(trait.names))
names(genus.distributions) <- trait.names
for (i in seq_along(trait.names)) {
  genus.distributions[[i]] <-
    data.frame(genus = genus.means[, "genus"],
               mean = genus.means[, (i + 1)],
               sd   = genus.sd[, (i + 1)],
               min  = genus.min[, (i + 1)],
               max  = genus.max[, (i + 1)],
               species.count = genus.species[, 2],
               species.observed = genus.observed[, (i + 1)]
               )
}

# Save distribution data for later reference
for (i in seq_along(genus.distributions)) {
  write.csv(genus.distributions[[i]],
            file = paste0("output/test/stochastic_gapfilling_genus_distributions_",
                          names(genus.distributions)[i],
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = FALSE
            )
}


##################################################
# Generate 100 stochastic gapfilled trait matrices
##################################################
cat("Generating 100 stochastic gapfilled trait matrices...\n")

if (!dir.exists("output/test/stochastic_gapfilled/")) {
  cat("creating directory: output/test/stochastic_gapfilled/\n")
  dir.create("output/test/stochastic_gapfilled/")
}

for (sample in seq_len(100)) {
  cat("sample", sample, "\n")
  traits.filled <- palm.traits.subset
  for (i in seq_along(traits.filled[, -c(1, 2)])) {
    traits.filled[, (i + 2)] <-
      StochasticFill(trait = traits.filled[, (i + 2)],
                     genus = traits.filled[, "genus"],
                     distribution = genus.distributions[[i]]
                     )
  }
  write.csv(traits.filled,
            file = paste0("output/test/stochastic_gapfilled/",
                          "stochastic_gapfilled_traits_",
                          sample,
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = FALSE
            )
}


###################################################
# Statistics of the observed genus-level trait data
###################################################
cat("Generating plots of observed genus-level trait values...\n")
# The subsetted dataset contains 122 Genera with 2-379 species (total 2479 species)

# Histogram of proportion of species with known trait values
# ----------------------------------------------------------
completeness <- llply(genus.distributions,
                      function(x) {
                        (x$species.observed / x$species.count)
                      }
                      )
GraphSVG(MultiHist(completeness, id = names(completeness), xlab = "completeness"),
         file = "graphs/test/genus_distribution_completeness.svg",
         width = 3,
         height = 4
         )
# Most genera have relatively complete data.


# Shape of genus-level trait distributions
# ----------------------------------------
# Can we make histograms for 122 genera, for three traits?
# I would, but my supervisor isn't going to look at all of them...
# We can make a list of observed trait values by genus:
genus.traits <- dlply(palm.traits.subset, "genus", function(x) { x } )
names(genus.traits) <- unique(palm.traits.subset[, "genus"])

# Then subset to genera with more than one growthform.
# First, find growthforms for all species
growthform <- numeric(length = length(trait.data$SpecName))
growthform[which(trait.data$Climbing == 1)] <- "climbing"
growthform[which(trait.data$Acaulescence == 1)] <- "acaulescent"
growthform[which(trait.data$Errect == 1)] <- "freestanding"
growthform[which(growthform == 0)] <- NA
# Second, subset to species in palm.traits.subset
indices <- CrossCheck(x = palm.traits[, "species"],
                      y = palm.traits.subset[, "species"],
                      presence = TRUE,
                      value = FALSE
                      )
growthform %<>% .[indices]
# Third, combine with genus names and find uniques
genus.growthform <- data.frame(genus = palm.traits.subset$genus,
                               growthform = growthform
                               )
genus.growthform %<>% count()
# Fourth, count growthforms per genus and subset to counts > 1
growthform.counts <- ddply(genus.growthform, "genus", nrow)
growthform.counts %<>% { .[which(.[, 2] > 1), ] }
# And finally, subset genus.traits
genus.traits %<>% .[names(.) %in% growthform.counts[, "genus"]]
# This limits us to a mere 33 genera.

# There is no obvious way to subset it further, so on with the plotting!
# First obtain for each trait, a list of values for each genus,
# with the NAs removed,
# and log10-transformed to compress the x-axis
Temp <- function(x, trait) {
  log10(x[!is.na(x[, trait]), trait] + 1)
}
genus.traits.height <- llply(genus.traits, Temp, trait = "stem.height")
genus.traits.blade <- llply(genus.traits, Temp, trait = "blade.length")
genus.traits.fruit <- llply(genus.traits, Temp, trait = "fruit.length")
# Then make the plots:
GraphPNG(MultiHist(genus.traits.height,
                   id = names(genus.traits.height),
                   xlab = "Observed log10(stem height)"
                   ),
         file = "graphs/test/genus_distribution_height.png",
         width = 480,
         height = (120 * 33)
         )
GraphPNG(MultiHist(genus.traits.blade,
                   id = names(genus.traits.blade),
                   xlab = "Observed log10(blade length)"
                   ),
         file = "graphs/test/genus_distribution_blade.png",
         width = 480,
         height = (120 * 33)
         )
GraphPNG(MultiHist(genus.traits.fruit,
                   id = names(genus.traits.fruit),
                   xlab = "Observed log10(fruit length)"
                   ),
         file = "graphs/test/genus_distribution_fruit.png",
         width = 480,
         height = (120 * 33)
         )




# Reset RNG
set.seed(NULL)

cat("Done.\n")
