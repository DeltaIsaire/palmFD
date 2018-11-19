###########################################
# Palm FD project: gap-filling trait matrix
###########################################
#
# In which the the three trait variables we will consider in this study
# (stem height, blade length, fruit length) are extracted from the trait
# dataset, and missing values in these three variables are filled using
# estimates based on (1) Genus-level means or (2) BHPMF.
# Input files:
#   data/palms_in_tdwg3.csv
#   data/PalmTraits_10.csv
#   data/TREE.nex
# Generated output files:
#   output/palm.traits.csv
#   output/palm.traits.genus.mean.csv

library(ape)
library(plyr)
library(BHPMF)

# -----------------
# Read raw datasets
# -----------------
palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")
# Long-list Data frame with two columns:
#	1. Area_code_L3 - 3-letter code for each botanical country
#	2. SpecName - Species name
trait.data <- read.csv(file="data/PalmTraits_10.csv")
# Ensure alphabetical ordering by species:
trait.data <- trait.data[order(trait.data$SpecName), ]
# Data frame with 31 columns:
#	1. SpecName - binomial species name [sorted alphabetically]
#	2. accGenus - only the genus name
#	3. accSpecies - only the species name
#	4. PalmTribe
#	5. PalmSubfamily
#	6-30 - information on 25 traits
#	31. extra reference notes
palm.tree <- read.nexus(file="data/TREE.nex")
# The tree file is nexus format. palm.tree is a phylo object.
# see str(palm.tree) for components.


# ------------------------------------------------
# Compare the list of palm species of each dataset
# ------------------------------------------------
# Do the three datasets agree on which palm species exist?
# 1. Extract species names.
# Ensure that each variable uses underscores instead of spaces.
species <- list(palm.dist  = levels(palm.dist$SpecName),
                trait.data = sub(pattern=" ", replacement="_",
                                 x=trait.data$SpecName), 
                palm.tree  = palm.tree$tip.label)
# 2. cross-check
missing.table <- MultiCC(species, presence=FALSE, value=FALSE)
missing.list <- MultiCC(species, presence=FALSE, value=TRUE)
# palm.dist and trait.data are the same, but palm.tree differs from the others.
# 3. Find the list of species upon which all three sources agree.
present.list <- MultiCC(species, presence=TRUE, value=TRUE)
# Verification test:
new.species <- list(palm.dist = present.list$palm.dist$palm.tree,
                    palm.tree = present.list$palm.tree$palm.dist)
MultiCC(new.species, presence=FALSE, value=FALSE)
# Yup, nothing missing.
species.agreed <- sort(present.list$palm.dist$palm.tree)


# -------------------------------------------------------
# Subset palm traits dataset to the three selected traits
# -------------------------------------------------------
palm.traits <- data.frame(species      = sub(pattern=" ", replacement="_",
                                             x=trait.data$SpecName),
                          genus        = trait.data$accGenus,
                          stem.height  = trait.data$MaxStemHeight_m,
                          blade.length = trait.data$Max_Blade_Length_m,
                          fruit.length = trait.data$AverageFruitLength_cm)
# stem.height = maximum stem height in meters.
# blade.length = maximum blade length in meters.
# fruit.length = average fruit length in centimeters.
# Genus column is included because we need it to calculate genus mean trait values

# DISABLED SUBSETTING TO AGREED SPECIES LIST
# palm.traits <- palm.traits[CrossCheck(palm.traits$species,
#                                       species.agreed, presence=TRUE,
#                                       value=FALSE), ]

write.csv(palm.traits, file="output/palm.traits.csv", eol="\r\n")


# ---------------------------------------
# Gap-fill trait matrix using genus means
# ---------------------------------------
genus.means <- ddply(palm.traits, "genus", numcolwise(mean, na.rm=TRUE))
traits.filled.means <- GapFill(palm.traits, genus.means, by="genus",
                               fill=c("stem.height", "blade.length",
                                      "fruit.length"))
# Beware NAs for genus means where all species of the genus have NA for the
# same trait.
# Solution: delete all cases where genus-mean filling is not possible.
# First extract these cases to a seperate dataframe:
genus.mean.missing <- traits.filled.means[!complete.cases(traits.filled.means), ]
# Then remove them from the filled trait matrix:
traits.filled.means <- traits.filled.means[complete.cases(traits.filled.means), ]

write.csv(traits.filled.means, file="output/palm.traits.genus.mean.csv",
          eol="\r\n")
write.csv(genus.mean.missing, file="output/genus.mean.missing.csv", eol="\r\n")


# ---------------------------------
# Gap-fill trait matrix using BHPMF
# ---------------------------------
# First assemble the prerequisites data matrices.
# Source dataframes are sorted by species name, so rows will correspond
# automatically.
trait.matrix <- as.matrix(palm.traits[, c("stem.height", "blade.length",
                                          "fruit.length")]
                          )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
hierarchy.matrix <- as.matrix(data.frame(species   = sub(pattern=" ",
                                                         replacement="_",
                                                         x=trait.data$SpecName),
                                         genus     = trait.data$accGenus,
                                         tribe     = trait.data$PalmTribe,
                                         subfamily = trait.data$PalmSubfamily)
                              )
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
# First extract these cases to a seperate dataframe:
to.remove <- which(rowSums(is.na(trait.matrix)) == ncol(trait.matrix))
BHPMF.missing <- palm.traits[to.remove, ]
# Then remove them from the matrices:
trait.matrix <- trait.matrix[-to.remove, ]
hierarchy.matrix <- hierarchy.matrix[-to.remove, ]
rm(to.remove)
# Run BHPMF. The function is BHPMF::GapFilling
# GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
#            used.num.hierarchy.levels=3,
#  mean.gap.filled.output.path="output/",
# std.gap.filled.output.path="output/")



# TODO: complete gap-filling with BHPMF
# This will be all complicated because the GapFilling function does not directly
# accept the phylogenetic tree in Nexus format.
# The tree has to be converted to a hierarchy matrix.
# see also ape::makeNodeLabel; ape::nodelabels; ape::plot.phylo
# Or, we use the taxonomic data from the trait dataset.

# TODO: compare genus-mean filling with BHPMF filling.
# One thing to note is that genus.mean.missing is a lot shorter than
# BHPMF.missing


# TODO: Subset all datasets to the list of agreed species.
#   Can we remove species from the phylogenetic tree? How?
#   Look into ape::drop.tip



