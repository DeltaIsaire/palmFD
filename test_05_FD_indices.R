#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we explore how to properly implement calculation of FD indices.
# Keep in mind the whole procedure is done twice, for two different datasets:
# The gap-filled palm traits data, and the unfilled palm traits data with complete
# cases only.
#
# Input files:
#   output/test/traits_filled_genus_mean.csv
#   output/test/palm_traits.csv
#   data/palms_in_tdwg3.csv
# Generated output files:
#   output/test/test_palm_traits_transformed_gapfilled.csv
#   output/test/test_palm_traits_transformed_unfilled.csv
#   output/test/test_palm_tdwg3_pres_abs_gapfilled.csv
#   output/test/test_palm_tdwg3_pres_abs_unfilled.csv
#   output/test/test_fd_indices_gapfilled.csv
#   output/test/test_fd_indices_unfilled.csv
#   output/test/test_tdwg3_trait_means_gapfilled.csv
#   output/test/test_tdwg3_trait_means_unfilled.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")

# Subset FD input for quick code testing?
subset <- FALSE


####################################################
# Data preparation: trait matrix and pres/abs matrix
####################################################
cat("Preparing data...\n")
# -----------------------------------------------------------------
# First thing we need is the (gapfilled and unfilled) trait matrix.
# -----------------------------------------------------------------
# The gapfilled dataset we use is the one from genus-mean filling.
# The matrices need 'species labels'. I'm assuming that means rownames.
traits.gapfilled <- read.csv(file = "output/test/traits_filled_genus_mean.csv")
traits.unfilled <-
   read.csv(file = "output/test/palm_traits.csv") %>%
  .[complete.cases(.), ]

trait.names <- c("stem.height", "blade.length", "fruit.length")

ToMatrix <- function(x) {
# Function to transform trait dataset x into a ready-to-use trait matrix,
# where x is the two datasets we just loaded.
  data <- x[, trait.names]
    # The traits should be log10-transformed.
  # The log of 0 is undefined (-Inf), so we cannot have stem height == 0.
  # In addition, The log10 of values below 1 is negative, and negative trait
  # values are weird. Or put more accurately, log transformation is symmetrical
  # around unity, not around zero.
  # Solution:
  # Add +1 to all observed trait values prior to log10-transformation
  data <-
    data.frame(stem.height  = log10(data[, "stem.height"] + 1),
               blade.length = log10(data[, "blade.length"] + 1),
               fruit.length = log10(data[, "fruit.length"] + 1)
               ) %>%
    as.matrix()
  rownames(data) <- x[, "species"]
  data
}

matrix.gapfilled <- ToMatrix(traits.gapfilled)
matrix.unfilled <- ToMatrix(traits.unfilled)

# ---------------------------------------------------
# Second thing we need is the presence/absence matrix
# ---------------------------------------------------
# Rows should be TDWG3 units and columns should be species.
palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")

PresAbs <- function(trait.matrix, dist.data = palm.dist) {
# Function to create a presence/absence matrix for a given trait matrix and
# distribution data.
# Output is a list with [1] the new trait matrix and [2] the pres/abs matrix.

  # Subset palm.dist to the species in trait.matrix
  species.shared <- CrossCheck(x = unique(dist.data[, "SpecName"]),
                               y = rownames(trait.matrix),
                               presence = TRUE,
                               value = TRUE
                               )
  indices <- CrossCheck(x = dist.data[, "SpecName"],
                        y = species.shared,
                        presence = TRUE,
                        value = FALSE
                        )
  new.dist.data <- dist.data[indices, ]
  # Verify: do species in trait.matrix and palm.dist match?
  if (!identical(length(rownames(trait.matrix)),
                 length(unique(new.dist.data[, "SpecName"]))
                 )
      ) {
    stop("trait.matrix and dist.data have differing number of species")
  }
  dist.species <-
    unique(new.dist.data[, "SpecName"]) %>%
    sort() %>%
    as.character()
  if (!identical(rownames(trait.matrix), dist.species)) {
    stop("species in trait.matrix and dist.data do not match")
  }
  # with palm.dist subsetted correctly, we transform it to pres/abs matrix:
  pres.abs.matrix <-
    table(new.dist.data) %>%
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
  list(trait.matrix, pres.abs.matrix)
}

gapfilled.data <- PresAbs(matrix.gapfilled)
unfilled.data <- PresAbs(matrix.unfilled)

# Parse and save results
# ----------------------
matrix.gapfilled <- gapfilled.data[[1]]
write.csv(matrix.gapfilled,
          file = "output/test/test_palm_traits_transformed_gapfilled.csv",
          eol = "\r\n",
          row.names = TRUE
          )

matrix.unfilled <- unfilled.data[[1]]
write.csv(matrix.unfilled,
          file = "output/test/test_palm_traits_transformed_unfilled.csv",
          eol = "\r\n",
          row.names = TRUE
          )

pres.abs.gapfilled <- gapfilled.data[[2]]
write.csv(pres.abs.gapfilled,
          file = "output/test/test_palm_tdwg3_pres_abs_gapfilled.csv",
          eol = "\r\n",
          row.names = TRUE
          )

pres.abs.unfilled <- unfilled.data[[2]]
write.csv(pres.abs.unfilled,
          file = "output/test/test_palm_tdwg3_pres_abs_unfilled.csv",
          eol = "\r\n",
          row.names = TRUE
          )


##################################
# Calculating Functional Diversity
##################################
cat("Calculating Functional Diversity indices... (this may take a while)\n")

cat("(1) for the gap-filled dataset:\n")
output.all <- RunFD(matrix.gapfilled,
                    pres.abs.gapfilled,
                    subset = subset,
                    verbose = TRUE
                    )
output.single <- SingleFD(matrix.gapfilled,
                          pres.abs.gapfilled,
                          subset = subset,
                          verbose = TRUE
                          )
fd.indices.gapfilled <- c(list(all.traits = output.all), output.single)

cat("(2) For the unfilled dataset:\n")
output.all <- RunFD(matrix.unfilled,
                    pres.abs.unfilled,
                    subset = subset,
                    verbose = TRUE
                    )
output.single <- SingleFD(matrix.unfilled,
                          pres.abs.unfilled,
                          subset = subset,
                          verbose = TRUE
                          )
fd.indices.unfilled <- c(list(all.traits = output.all), output.single)

# Structure of FD output list:
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


########################
# Parse and save results
########################
cat("Parsing and saving output...\n")

ParseOutput <- function(output, name) {
# Where output is the output list from RunFD.
# Name is an identifier character string, should be either "gapfilled" or
# "unfilled".

  # Turns out it is best to combine all FD data into a single dataframe / file.
  fd.indices <- data.frame(FRic.all.traits   = output[[1]]$FRic,
                           FRic.stem.height  = output[[2]]$FRic,
                           FRic.blade.length = output[[3]]$FRic,
                           FRic.fruit.length = output[[4]]$FRic,
                           FDis.all.traits   = output[[1]]$FDis,
                           FDis.stem.height  = output[[2]]$FDis,
                           FDis.blade.length = output[[3]]$FDis,
                           FDis.fruit.length = output[[4]]$FDis,
                           row.names = names(output[[1]]$nbsp)
                           )
  write.csv(fd.indices,
            file = paste0("output/test/test_fd_indices_", name, ".csv"),
            eol = "\r\n",
            row.names = TRUE
            )
  # Community weighted mean trait values
  # Not weighted in our case, so just the mean
  # KEEP IN MIND THESE ARE LOG10-TRANSFORMED VALUES
  community.means <- output[[1]]$CWM
  write.csv(community.means,
            file = paste0("output/test/test_tdwg3_trait_means_", name, ".csv"),
            eol = "\r\n",
            row.names = TRUE
            )
  return (0)
}

ParseOutput(fd.indices.gapfilled, "gapfilled")
ParseOutput(fd.indices.unfilled, "unfilled")

cat("Done.\n")

