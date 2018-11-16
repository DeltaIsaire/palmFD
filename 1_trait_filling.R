###########################################
# Palm FD project: gap-filling trait matrix
###########################################
#
# In which <description>
# Input files:
#   data/palms_in_tdwg3.csv
#   data/PalmTraits_10.csv
#   data/TREE.nex
# Generated output files:
#   output/palm.traits.csv
#   output/palm.traits.genus.mean.csv

library(ape)
library(plyr)


# Read palm distribution data
palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")
# Long-list Data frame with two columns:
#	1. Area_code_L3 - 3-letter code for each botanical country
#	2. SpecName - Species name

# Read palm trait data
trait.data <- read.csv(file="data/PalmTraits_10.csv")
# Data frame with 31 columns:
#	1. SpecName - binomial species name [sorted alphabetically]
#	2. accGenus - only the genus name
#	3. accSpecies - only the species name
#	4. PalmTribe
#	5. PalmSubfamily
#	6-30 - information on 25 traits
#	31. extra reference notes

# Read phylogenetic tree
palm.tree <- read.nexus(file="data/TREE.nex")
# It's in nexus format. see str(palm.tree) for components


# Explore if these three sources agree on which palm species exist
# 1. Extract species names
# Ensure that each variable uses underscores instead of spaces
species <- list(palm.dist  = levels(palm.dist$SpecName),
                trait.data = sub(pattern=" ", replacement="_",
                                 x=trait.data$SpecName), 
                palm.tree  = palm.tree$tip.label)
# 2. cross-check
missing.table <- multiCC(species, presence=FALSE, value=FALSE)
missing.list <- multiCC(species, presence=FALSE, value=TRUE)
# palm.dist and trait.data are the same, but palm.tree differs from the others.
# 3. Find the list of species upon which all three sources agree
present.list <- multiCC(species, presence=TRUE, value=TRUE)
# Just to verify:
new.species <- list(palm.dist = present.list$palm.dist$palm.tree,
                    palm.tree = present.list$palm.tree$palm.dist)
multiCC(new.species, presence=FALSE, value=FALSE)
# Yup, nothing missing.
species.agreed <- sort(present.list$palm.dist$palm.tree)


# Subset palm traits dataset
palm.traits <- data.frame(species      = sub(pattern=" ", replacement="_",
                                             x=trait.data$SpecName),
                          genus        = trait.data$accGenus,
                          stem.height  = trait.data$MaxStemHeight_m,
                          blade.length = trait.data$Max_Blade_Length_m,
                          fruit.length = trait.data$AverageFruitLength_cm)
# stem.height = maximum stem height in meters.
# blade.length = maximum blade length in meters.
# fruit.length = average fruit length in centimeters.
palm.traits <- palm.traits[crossCheck(palm.traits$species,
                                      species.agreed, presence=TRUE,
                                      value=FALSE), ]
write.csv(palm.traits, file="output/palm.traits.csv", eol="\r\n")


# Calculate Genus-level mean trait values and fill trait matrix
genus.means <- ddply(palm.traits, "genus", numcolwise(mean, na.rm=TRUE))
traits.filled.means <- gapFill(palm.traits, genus.means, by="genus",
                               fill=c("stem.height", "blade.length",
                                      "fruit.length"))
# Beware NAs for genus means where all species of the genushave NA for the
# same trait.
# Solution: delete all cases where genus-mean filling is not possible
# That is a drawback BHPMF *may* not have.
# TODO: discuss this solution
genus.mean.missing <- traits.filled.means[!complete.cases(traits.filled.means), ]
traits.filled.means <- traits.filled.means[complete.cases(traits.filled.means), ]
write.csv(traits.filled.means, file="output/palm.traits.genus.mean.csv",
          eol="\r\n")
write.csv(genus.mean.missing, file="output/genus.mean.missing.csv", eol="\r\n")






# TODO: Subset all datasets to the list of agreed species.
#   Can we remove species from the phylogenetic tree? How?



