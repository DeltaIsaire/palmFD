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
#   output/BHPMF_preprocessing (directory containing many files)
#   output/BHPMF.mean.gap.filled.txt
#   output/BHPMF.std.gap.filled.txt
#   output/palm.traits.BHPMF.csv
#   output/BHPMF.missing.csv


library(ape)
library(plyr)
library(BHPMF)

source(file="functions/base.functions.R")


# -----------------
# Read raw datasets
# -----------------
print("preparing data...",quote=FALSE)
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
# First extract species names.
# Ensure that each variable uses underscores instead of spaces.
species <- list(palm.dist  = levels(palm.dist$SpecName),
                trait.data = sub(pattern=" ", replacement="_",
                                 x=trait.data$SpecName), 
                palm.tree  = palm.tree$tip.label)
# For each dataset, extract the list of species that are present in the other
# datasets
present.list <- MultiCC(species, presence=TRUE, value=TRUE)
# palm.dist and trait.data are the same, but palm.tree differs from the others.
# Therefore, the list of species present in all 3 datasets is:
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
# Genus is included because we need it to calculate genus mean trait values
write.csv(palm.traits, file="output/palm.traits.csv", eol="\r\n")


# ---------------------------------------
# Gap-fill trait matrix using genus means
# ---------------------------------------
print("gap-filling with genus means...",quote=FALSE)
genus.means <- ddply(palm.traits, "genus", numcolwise(mean, na.rm=TRUE))
traits.filled.mean <- GapFill(palm.traits, genus.means, by="genus",
                               fill=c("stem.height", "blade.length",
                                      "fruit.length")
                               )
# Beware NAs for genus means where all species of the genus have NA for the
# same trait.
# Solution: delete all cases where genus-mean filling is not possible.
# First extract these cases to a seperate dataframe:
genus.mean.missing <- traits.filled.mean[!complete.cases(traits.filled.mean), ]
# Then remove them from the filled trait matrix:
traits.filled.mean <- traits.filled.mean[complete.cases(traits.filled.mean), ]

write.csv(traits.filled.mean, file="output/palm.traits.genus.mean.csv",
          eol="\r\n")
write.csv(genus.mean.missing, file="output/genus.mean.missing.csv", eol="\r\n")


# ---------------------------------
# Gap-fill trait matrix using BHPMF
# ---------------------------------
print("gap-filling with BHPMF... (this may take a while)",quote=FALSE)
# First assemble the prerequisites data matrices.
# Source dataframes are sorted by species name, so rows will correspond
# automatically.
trait.matrix <- as.matrix(palm.traits[, c("stem.height", "blade.length",
                                          "fruit.length")
                                      ]
                          )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
# columns of the hierarchy matrix should be from low-level to high-level:
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
# BHPMF is unable to handle observed trait values of exactly 0.
# These do occur in our dataset, for stem height.
# Solution: set these values to 0.01 (corresponding to stem height of 1 cm,
# which is sufficiently close to zero for biological purposes)
trait.matrix[, "stem.height"][which(trait.matrix[, "stem.height"] == 0)] <- 0.01
# BHPMF wants a preprocessing directory, where it saves pre-processed files.
# To avoid errors or erronous output when re-running the code, this directory
# needs to be emptied.
unlink("output/BHPMF_preprocessing", recursive=TRUE)
dir.create("output/BHPMF_preprocessing")

# Run BHPMF. The function is BHPMF::GapFilling.
# This step is computationally intensive!
GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
           prediction.level=4, used.num.hierarchy.levels=3,
           mean.gap.filled.output.path="output/BHPMF.mean.gap.filled.txt",
           std.gap.filled.output.path="output/BHPMF.std.gap.filled.txt",
           tmp.dir="output/BHPMF_preprocessing", 
           rmse.plot.test.data=FALSE, verbose=FALSE)

# Extract BHPMF output into a neat dataframe:
mean.BHPMF <- read.table(file="output/BHPMF.mean.gap.filled.txt",
                         header=TRUE, sep="	")
traits.filled.BHPMF <- data.frame(species = as.character(hierarchy.matrix[, 1]),
                                  genus   = as.character(hierarchy.matrix[, 2]),
                                  mean.BHPMF)
write.csv(traits.filled.BHPMF, file="output/palm.traits.BHPMF.csv",
          eol="\r\n")
write.csv(BHPMF.missing, file="output/BHPMF.missing.csv", eol="\r\n")



# TODO:we use the taxonomic data from the trait dataset.
# Will we do nothing with the phylogenetic tree?
# Also, consider tuning (see ?GapFilling)

# TODO: compare genus-mean filling with BHPMF filling.
# One thing to note is that genus.mean.missing is a lot shorter than
# BHPMF.missing
# NOTE: it may be possible to re-run BHPMF with additional traits from the trait
# dataset. While we will not investigate those traits, they might allow BHPMF
# to estimate missing values for our three traits of interest, even for species
# where all three trait values are missing.


# TODO: Subset all datasets to the list of agreed species.
#   Can we remove species from the phylogenetic tree? How?
#   Look into ape::drop.tip



