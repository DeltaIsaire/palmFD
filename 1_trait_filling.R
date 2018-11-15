###########################################
# Palm FD project: gap-filling trait matrix
###########################################
#
# In which (...)
# Input files:
#   data/palms_in_tdwg3.csv
#   data/PalmTraits_10.csv
#   data/TREE.nex
# Generated output files:
#   none

library(ape)


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
tree <- read.nexus(file="data/TREE.nex")


# Explore if these three sources agree on which palm species exist
# 1. Extract species names
# Ensure that each variable uses underscores instead of spaces
species <- list(
  palm.dist  = levels(palm.dist$SpecName), 
  trait.data = sub(pattern=" ", replacement="_", x=trait.data$SpecName), 
  tree       = tree$tip.label)
# 2. cross-reference
missing.table <- multiCC(species, presence=FALSE, value=FALSE)
missing.list <- multiCC(species, presence=FALSE, value=TRUE)
# palm.dist and trait.data are the same, but tree differs from the others.
# 3. Find the list of species upon which all three sources agree
present.list <- multiCC(species, presence=TRUE, value=TRUE)
# Just to verify:
new.species <- list(
  palm.dist = present.list$palm.dist$tree,
  tree      = present.list$tree$palm.dist)
multiCC(new.species, presence=FALSE, value=FALSE)
# Yup, nothing missing
species.agreed <- sort(present.list$palm.dist$tree)


# TODO: Subset all datasets to the list of agreed species.
# Can we remove species from the phylogenetic tree? How?



