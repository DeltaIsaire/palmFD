#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we explore how to properly implement calculation of FD indices.
#
# Input files:
#   output/traits_filled_genus_mean.csv
#   data/palms_in_tdwg3.csv
# Generated output files:
#   output/test_fd.indices.csv
#   output/test_fd.indices_single_traits.csv


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
# For now, we'll use genus-mean filled data.
# It needs 'species labels'. I'm assuming that means rownames.
trait.names <- c("stem.height", "blade.length", "fruit.length")
traits.filled <- read.csv(file = "output/traits_filled_genus_mean.csv")
trait.matrix <- traits.filled[, trait.names]

# The traits should be log10-transformed.
# The log of 0 is undefined (-Inf), so we cannot have stem height == 0.
# In addition log10 for values below 1 is negative, which I feel better avoiding.
# Solution:
# Transform stem height from m to cm,
# Transform blade length from m to cm,
# Transform fruit length from cm to mm.
# After log-transformation, substitute -Inf with 0.
trait.matrix <- 
  data.frame(stem.height  = log10(trait.matrix[, 1] * 100),
             blade.length = log10(trait.matrix[, 2] * 100),
             fruit.length = log10(trait.matrix[, 3] * 10)
             ) %>%
  as.matrix()
trait.matrix[which(trait.matrix == -Inf)] <- 0
rownames(trait.matrix) <- traits.filled$species

# Second thing we need is a presence/absence matrix, where rows are TDWG3
# units and columns are species.
palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")
# Subset to the species for which we have completely filled trait data
species.shared <- CrossCheck(x = unique(palm.dist$SpecName),
                             y = traits.filled$species,
                             presence = TRUE,
                             value = TRUE
                             )
indices <- CrossCheck(x = palm.dist$SpecName,
                      y = species.shared,
                      presence = TRUE,
                      value = FALSE
                      )
palm.dist %<>% .[indices, ]
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
cat("Using trait dataset with",
    dim(pres.abs.matrix)[2],
    "species occurring in",
    dim(pres.abs.matrix)[1],
    "botanical countries.\n"
    )


# --------------------------------
# Calculating Functional Diversity
# --------------------------------
cat("Calculating Functional Diversity indices... (this may take a moment)\n")

# Subset data to expedite testing. This stuff is computationally intensive.
#pres.abs.matrix <- pres.abs.matrix[1:10, ]
#orphaned.species <-
#  which(colSums(pres.abs.matrix) == 0) %>%
#  colnames(pres.abs.matrix)[.]
#pres.abs.matrix %<>% .[, -which(colSums(.) == 0)]
#indices <- CrossCheck(x = rownames(trait.matrix),
#                      y = orphaned.species,
#                      presence = TRUE,
#                      value = FALSE
#                      )
#trait.matrix <- trait.matrix[-indices, ]

# FD using all traits
# -------------------
cat("(1) FD with all traits...\n")
output.all <- dbFD(x = trait.matrix,
                   a = pres.abs.matrix,
                   w.abun = FALSE,
                   stand.x = TRUE,
                   corr = "cailliez",  # Is this the best option?
                   calc.FRic = TRUE,
                   m = "max",
                   stand.FRic = TRUE,  # Would this be useful to do?
                   calc.CWM = TRUE,
                   calc.FDiv = FALSE,
                   messages = TRUE
                   )
# Structure of output.all:
# $nbsp      - number of species in each community.
# $sing.sp   - number of species with a unique trait combination in each community.
# $FRic      - Functional Richness of each community. WE WANT THIS.
# $qual.FRic - Quality of reduced-space representation of FD. In our case, space
#              reduction should not be necessary, so this value should be 1.
# $FEve      - Functional Evenness for each community.
# $FDis      - Functional Dispersion of each community. WE WANT THIS.
# $RaoQ      - Rao's quadratic entropy (Q) of each community.
# $CWM       - Data frame with community-weighted mean trait values. We use
#              presence/absence data so there is no weighing in our case.
#              So for us, this is the mean trait value in each community.
#              COULD BE USEFUL.

# FD for single traits
# --------------------
cat("(2) FD for single traits...\n")

Temp <- function(trait) {
  cat(trait, "...\n")
  subset.matrix <-
     trait.matrix[, trait] %>%
     as.matrix()
  colnames(subset.matrix) <- trait
  fd.trait <- dbFD(x = subset.matrix,
                   a = pres.abs.matrix,
                   w.abun = FALSE,
                   stand.x = TRUE,
                   corr = "cailliez",  # Is this the best option?
                   calc.FRic = TRUE,
                   m = "max",
                   stand.FRic = TRUE,  # Would this be useful to do?
                   calc.CWM = TRUE,
                   calc.FDiv = FALSE,
                   messages = TRUE
                   )
  fd.trait
}

output.stem.height  <- suppressWarnings(Temp(trait.names[1]))
output.blade.length <- suppressWarnings(Temp(trait.names[2]))
output.fruit.length <- suppressWarnings(Temp(trait.names[3]))
# Zero distance(s)' warning means some species within the same community
# had identical coordinates in functional space. This is not a problem
# so we can ignore these warnings.


# ----------------------
# Parse and save results
# ----------------------
cat("Parsing and saving output...\n")

# FD indices for all traits
fd.indices <- data.frame(TDWG3 = names(output.all$nbsp),
                         FRic  = output.all$FRic,
                         FDis  = output.all$FDis
                         )
write.csv(fd.indices,
          file = "output/test_fd.indices.csv",
          eol="\r\n",
          row.names=FALSE
          )

# FD indices for single traits
single.fd <- data.frame(TDWG3       = names(output.all$nbsp),
                        height.FRic = output.stem.height$FRic,
                        height.FDis = output.stem.height$FDis,
                        blade.FRic  = output.blade.length$FRic,
                        blade.FDis  = output.blade.length$FDis,
                        fruit.FRic  = output.fruit.length$FRic,
                        fruit.FDis  = output.fruit.length$FDis
                        )
write.csv(single.fd,
          file = "output/test_fd.indices_single_traits.csv",
          eol="\r\n",
          row.names=FALSE
          )

# Community weighted mean trait values
# Not weighted in our case, so just the mean
# KEEP IN MIND THESE ARE LOG10-TRANSFORMED VALUES
community.means <- output.all$CWM
write.csv(community.means,
          file = "output/test_tdwg3_trait_means.csv",
          eol="\r\n",
          row.names=FALSE
          )

cat("Done.\n")

