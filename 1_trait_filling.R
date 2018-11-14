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
species <- list(
  palm.dist  = levels(palm.dist$SpecName), 
  trait.data = trait.data$SpecName, 
  tree       = tree$tip.label)
# 2. cross-reference, using crossReference function
crossReference(species$tree, species$palm.dist, presence=FALSE, value=FALSE)


# TODO: implement wrapper function for crossReference, which takes a list and cross-references all list elements with each other. Maybe lapply or plyr::**ply can help.

a <- plyr::llply(species, crossReference, y=species$palm.dist,
  presence=FALSE, value=FALSE)

# TODO: what you want to end up with is a table listing, for each element,
# the number of entries present and/or missing from the list it is compared to.
str(as.vector(a[3]))
# it's a list, and lists are special. Use unlist to extract the values:
b <- unlist(a[3],use.names=FALSE)
# and then you can do stuff like
length(b)

