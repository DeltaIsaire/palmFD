#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we explore how to properly implement calculation of FD indices.
#
# Input files:
#
# Generated output files:
#


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)

source(file = "functions/base_functions.R")


# --------------------------------------------------
# Data preparation: trait matrix and pres/abs matrix
# --------------------------------------------------
cat("Preparing data...\n")
# First thing we need is the (filled) trait matrix.
# It needs 'species labels'. I'm assuming that means rownames.
trait.matrix <- read.csv(file = "output/palm_traits.csv")
rownames(trait.matrix) <- trait.matrix$species
trait.matrix %<>% 
  .[, c("stem.height", "blade.length", "fruit.length")] %>%
  as.matrix()

# Second thing we need is a presence/absence matrix, where rows are TDWG3
# units and columns are species.
palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")
# TODO: implement subsetting (removal of palm species) as necessary.
# Verify: do species in trait.matrix and palm.dist match?
if (!length(rownames(trait.matrix)) == length(unique(palm.dist$SpecName))) {
  stop("trait.matrix and palm.dist have differing number of species")
}
dist.species <-
  unique(palm.dist$SpecName) %>%
  sort() %>%
  as.character()
if (!identical(rownames(trait.matrix), dist.species)) {
  stop("species in trait.matrix and palm.dist do not match")
}
# with palm.dist subsetted correctly, we transform it to pres/abs matrix:
pres.abs.matrix <-
  table(palm.dist) %>%
  as.data.frame.matrix() %>%  # yes, that's a function, and yes, we need it.
  as.matrix()
# Subset to TDWG3 units with richness > 3 as requirement for FD indices:
pres.abs.matrix %<>% .[which(rowSums(.) > 3), ]
# Remove TDWG3 unit 'VNA' because we have no environmental data for it:
pres.abs.matrix %<>% .[-which(rownames(.) == "VNA"), ]
# Make sure all species still occur in at least 1 community:
orphaned.species <-
  which(colSums(pres.abs.matrix) == 0) %>%
  colnames(pres.abs.matrix)[.]
if (length(orphaned.species) > 0) {
  pres.abs.matrix %<>% .[, -which(colSums(.) == 0)]
  indices <- CrossCheck(x = rownames(trait.matrix),
                        y = orphaned.species,
                        presence = TRUE,
                        value = FALSE
                        )
  trait.matrix %<>% .[-indices, ]
}


# --------------------------------
# Calculating Functional Diversity
# --------------------------------
cat("Calculating Functional Diversity indices...\n")

output <- dbFD(trait.matrix,
               pres.abs.matrix,
               stand.x = TRUE,  # Standardize traits?
               corr = "cailliez",  # Is this the best option?
               calc.FRic = TRUE,
               m = "max",
               stand.FRic = FALSE,  # Standardize FRic to the global volume?
               calc.CWM = TRUE,  # calculate community-weighted mean trait values
               calc.FDiv = TRUE
               )


cat("Done.\n")
