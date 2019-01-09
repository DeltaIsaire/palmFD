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
#   output/test/test_palm_traits_transformed_filled.csv
#   output/test/test_palm_traits_transformed_unfilled.csv
#   output/test/test_palm_tdwg3_pres_abs_filled.csv
#   output/test/test_palm_tdwg3_pres_abs_unfilled.csv
#   output/test/test_fd_indices_filled.csv
#   output/test/test_fd_indices_unfilled.csv
#   output/test/test_fd_indices_single_traits_filled.csv
#   output/test/test_fd_indices_single_traits_unfilled.csv
#   output/test/test_tdwg3_trait_means_filled.csv
#   output/test/test_tdwg3_trait_means_unfilled.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)

source(file = "functions/base_functions.R")


####################################################
# Data preparation: trait matrix and pres/abs matrix
####################################################
cat("Preparing data...\n")
# --------------------------------------------------------------
# First thing we need is the (filled and unfilled) trait matrix.
# --------------------------------------------------------------
# The filled dataset we use is the one from genus-mean filling.
# The matrices need 'species labels'. I'm assuming that means rownames.
traits.filled <- read.csv(file = "output/test/traits_filled_genus_mean.csv")
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
  # Transform stem height from m to cm,
  # Transform blade length from m to cm,
  # Transform fruit length from cm to mm.
  # After log-transformation, substitute -Inf with 0.
  data <-
    data.frame(stem.height  = log10(data[, "stem.height"] * 100),
               blade.length = log10(data[, "blade.length"] * 100),
               fruit.length = log10(data[, "fruit.length"] * 10)
               ) %>%
    as.matrix()
  data[which(data == -Inf)] <- 0
  rownames(data) <- x[, "species"]
  data
}

matrix.filled <- ToMatrix(traits.filled)
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

filled.data <- PresAbs(matrix.filled)
unfilled.data <- PresAbs(matrix.unfilled)

# Parse and save results
# ----------------------
matrix.filled <- filled.data[[1]]
write.csv(matrix.filled,
          file = "output/test/test_palm_traits_transformed_filled.csv",
          eol = "\r\n",
          row.names = TRUE
          )

matrix.unfilled <- unfilled.data[[1]]
write.csv(matrix.unfilled,
          file = "output/test/test_palm_traits_transformed_unfilled.csv",
          eol = "\r\n",
          row.names = TRUE
          )

pres.abs.filled <- filled.data[[2]]
write.csv(pres.abs.filled,
          file = "output/test/test_palm_tdwg3_pres_abs_filled.csv",
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

RunFD <- function(trait.matrix, pres.abs.matrix, subset = FALSE) {
# Wrapper function to run FD for the given dataset.
# Output is a list, where the first element is the output of FD using all traits,
# and further elements are the output of FD using single traits.
# If subset is true, FD is calculated only for the first 10 communities,
# which reduces the function runtime. Obviously that is for testing purposes.
  if (subset) {
    pres.abs <- pres.abs.matrix[1:10, ]
    orphaned.species <-
      which(colSums(pres.abs) == 0) %>%
      colnames(pres.abs)[.]
    pres.abs %<>% .[, -which(colSums(.) == 0)]
    indices <- CrossCheck(x = rownames(trait.matrix),
                          y = orphaned.species,
                          presence = TRUE,
                          value = FALSE
                          )
    traits <- trait.matrix[-indices, ]
  } else {
  pres.abs <- pres.abs.matrix
  traits <- trait.matrix
  }
  cat("Using trait dataset with",
      dim(pres.abs)[2],
      "species occurring in",
      dim(pres.abs)[1],
      "botanical countries.\n"
      )

  # FD using all traits
  # -------------------
  cat("(1) FD with all traits...\n")
  output.all <- suppressWarnings(dbFD(x = traits,
                                      a = pres.abs,
                                      w.abun = FALSE,
                                      stand.x = TRUE,
                                      corr = "cailliez",  # Is this the best option?
                                      calc.FRic = TRUE,
                                      m = "max",
                                      stand.FRic = FALSE,  # Would this be useful?
                                      calc.CWM = TRUE,
                                      calc.FDiv = FALSE,
                                      messages = TRUE
                                      )
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
  cat("(2) FD for single traits:\n")
  Temp <- function(trait) {
    cat(trait, "...\n")
    subset.matrix <-
      traits[, trait] %>%
      as.matrix()
    colnames(subset.matrix) <- trait
    fd.trait <- suppressWarnings(dbFD(x = subset.matrix,
                                      a = pres.abs,
                                      w.abun = FALSE,
                                      stand.x = TRUE,
                                      corr = "cailliez",  # Is this the best option?
                                      calc.FRic = TRUE,
                                      m = "max",
                                      stand.FRic = FALSE,  # Would this be useful?
                                      calc.CWM = FALSE,
                                      calc.FDiv = FALSE,
                                      messages = TRUE
                                      )
                                 )
    fd.trait
  }
  output.stem.height  <- suppressWarnings(Temp(trait.names[1]))
  output.blade.length <- suppressWarnings(Temp(trait.names[2]))
  output.fruit.length <- suppressWarnings(Temp(trait.names[3]))
  # Zero distance(s)' warning means some species within the same community
  # had identical coordinates in functional space. This is not a problem
  # so we can ignore these warnings.

  # Gather results:
  # ---------------
  list(output.all,
       output.stem.height,
       output.blade.length,
       output.fruit.length
       )
}

# ----------------------
# Run the RunFD function
# ----------------------
cat("(A) for the gap-filled dataset:\n")
fd.indices.filled <- RunFD(matrix.filled, pres.abs.filled, subset = FALSE)
cat("(B) for the unfilled dataset:\n")
fd.indices.unfilled <- RunFD(matrix.unfilled, pres.abs.unfilled, subset = FALSE)


########################
# Parse and save results
########################
cat("Parsing and saving output...\n")

ParseOutput <- function(output, name) {
# Where output is the output list from RunFD.
# Name is an identifier character string, should be either "filled" or "unfilled".

  # FD indices for all traits
  fd.indices <- data.frame(TDWG3 = names(output[[1]]$nbsp),
                           FRic  = output[[1]]$FRic,
                           FDis  = output[[1]]$FDis
                           )
  write.csv(fd.indices,
            file = paste0("output/test/test_fd_indices_", name, ".csv"),
            eol = "\r\n",
            row.names = FALSE
            )
  # FD indices for single traits
  single.fd <- data.frame(TDWG3       = names(output[[1]]$nbsp),
                          height.FRic = output[[2]]$FRic,
                          height.FDis = output[[2]]$FDis,
                          blade.FRic  = output[[3]]$FRic,
                          blade.FDis  = output[[3]]$FDis,
                          fruit.FRic  = output[[4]]$FRic,
                          fruit.FDis  = output[[4]]$FDis
                          )
  write.csv(single.fd,
            file = paste0("output/test/test_fd_indices_single_traits_",
                          name,
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = FALSE
            )
  # Community weighted mean trait values
  # Not weighted in our case, so just the mean
  # KEEP IN MIND THESE ARE LOG10-TRANSFORMED VALUES
  community.means <- output[[1]]$CWM
  write.csv(community.means,
            file = paste0("output/test/test_tdwg3_trait_means_", name, ".csv"),
            eol = "\r\n",
            row.names = FALSE
            )
  return (0)
}

ParseOutput(fd.indices.filled, "filled")
ParseOutput(fd.indices.unfilled, "unfilled")

cat("Done.\n")

