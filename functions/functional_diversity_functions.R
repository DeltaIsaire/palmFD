########################################################
# Palm FD project: custom functional diversity functions
########################################################
#
# In which custom (wrapper) functions for calculating functional richness (FRic)
# and functional dispersion (FDis) are defined.
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
# RunFD
# SingleFD


library(plyr)
library(magrittr)
library(FD)



RunFD <- function(trait.matrix,
                  pres.abs.matrix,
                  subset = FALSE,
                  verbose = TRUE
                  ) {
# A convenient wrapper for the function 'FD::dbFD', which automatically selects
# arguments suitable for our data.
#
# Args:
#   trait.matrix: A matrix with observed numerical trait values in columns,
#                 without NAs, and rownames with species names.
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be tdwg3 labels, and colnames should be species names.
#                    The colnames should match the rownames of the trait.matrix!
#  subset: Logical indicating whether to subset the data to only the first 10
#          tdwg3 units. This speeds up code execution, for testing purposes.
#  verbose: Logical indicating whether to print progress messages.
#
# Returns: A list with the output of function 'FD:dbFD' applied to the 
#          provided dataset.
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
  if (verbose) {
    cat("Using trait dataset with",
        dim(pres.abs)[2],
        "species occurring in",
        dim(pres.abs)[1],
        "botanical countries.\n"
        )
  }
  # function 'dbFD' can output lots of warnings with our data.
  # The 'Zero distance(s)' warning means some species within the same community
  # had identical coordinates in functional space. This is not a problem,
  # so we can ignore these warnings.
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
                                      messages = verbose
                                      )
                                 )
  output.all
}


SingleFD <- function(trait.matrix,
                     pres.abs.matrix,
                     subset = FALSE,
                     verbose = TRUE
                     ) {
# A convenient wrapper for the function 'FD::dbFD', which automatically selects
# arguments suitable for our data. SingleFD calculates FD indices for each trait
# individually (single-trait FD) and combines the result into a list.
#
# Args:
#   trait.matrix: A matrix with observed numerical trait values in columns,
#                 without NAs, and rownames with species names.
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be tdwg3 labels, and colnames should be species names.
#                    The colnames should match the rownames of the trait.matrix!
#  subset: Logical indicating whether to subset the data to only the first 10
#          tdwg3 units. This speeds up code execution, for testing purposes.
#  verbose: Logical indicating whether to print progress messages.
#
# Returns: A list of the same length as the number of traits in the trait matrix,
#          where each element is a list with the output of function 'FD:dbFD'
#          applied to the dataset subsetted to a single trait. 
  if (is.null(colnames(trait.matrix))) {
    stop("trait matrix should have colnames to identify traits")
  }
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
  if (verbose) {
    cat("Using trait dataset with",
        dim(pres.abs)[2],
        "species occurring in",
        dim(pres.abs)[1],
        "botanical countries.\n"
        )
  }
  output.list <- rep(list(NULL), ncol(trait.matrix))
  for (i in 1:ncol(trait.matrix)) {
    if (verbose) {
      cat("Calculating FD for trait", names(trait.matrix)[i], "...\n")
    }
    subset.matrix <- traits[, i, drop = FALSE]
    output.list[[i]] <- RunFD(trait.matrix = subset.matrix,
                              pres.abs.matrix = pres.abs,
                              subset = FALSE,
                              verbose = verbose
                              )
  }
  names(output.list) <- colnames(trait.matrix)
  output.list
}

