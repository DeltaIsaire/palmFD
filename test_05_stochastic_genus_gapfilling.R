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
genus.missing <- ddply(palm.traits.subset, "genus", numcolwise(CountObserved))

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
               species.missing = genus.missing[, (i + 1)]
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


# Reset RNG
set.seed(NULL)

cat("Done.\n")
