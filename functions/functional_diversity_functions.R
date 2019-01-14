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
# RandomSpecies


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


RandomSpecies <- function(communities, richness, species, trait.matrix,
                          verbose = TRUE) {
# A helpful function for generating null models. 
# For each community, get the species richness and sample that many
# species from the (group) species pool, without replacement,
# and combine results into a list.
#
# The FD requirement species > traits applies to functionally unique species.
# That means it does NOT suffice to have richness > 3 (our number of traits).
# Rather, we need > 3 species with unique trait combinations. 
# This happens to be the case for our real dataset, but is not automatically
# guaranteed for the null model samples.
# So we have to enforce this explicitly, by resampling until this condition
# is satisfied.
#
# If grouping is used, the groups in 'communities' and 'species' should be
# equal in number, and group names should match.
#
# Args:
#   communities: character vector giving community names for which to sample
#                species. Or, a list where each element gives the community
#                names belonging to a group of communities.
#   richness: dataframe where the first column is a vector of community names
#             and the second column is a vector of species richness for these
#             communities.
#   species: character vector giving the names of species in the species pool.
#            Or, a list where each element gives the species pool for a group
#            of communities.
#   trait.matrix: A matrix with observed numerical trait values in columns,
#                 without NAs, and rownames with species names. Used to
#                 verify the S > T criterion.
#  verbose: Logical indicating whether to print progress messages.
#
# Returns: 
#   A list giving the randomly sampled species names for each community

  # parse and check input
  if (is.list(communities)) {
    areas <- communities
  } else {
  areas <- list(group = communities)
  }
  if (is.list(species)) {
  list.species <- species
  } else {
  list.species <- list(group = species)
  }
  if (!identical(names(areas), names(list.species))) {
    stop("group names in lists 'communities' and 'species' do not match")
  }
  areas %<>% .[order(names(.))]
  list.species %<>% .[order(names(.))]
  # Run sampling routine
  do.sample <- TRUE
  sampling.iterations <- 0
  while (do.sample) {
    random.species <- rep(list(NULL), length = length(unlist(areas)))
    i <- 1
    for (group in seq_along(areas)) {
      for (area in seq_along(areas[[group]])) {
        area.richness <-
          richness %>% {
            .[which(.[1] == areas[[group]][area]), 2]
          }
        indices <-
          runif(n = area.richness * 3,
                min = 1,
                max = length(list.species[[group]])
                ) %>%
          round() %>%
          unique() %>%
          .[seq_len(area.richness)]
        random.species[[i]] <- list.species[[group]][indices]
        names(random.species)[i] <- areas[[group]][area]
        i <- i + 1
      }
    }
    random.species %<>% .[order(names(.))]
    # Check if uniqueness condition is satisfied:
    # Step 1: assemble trait values for species in each community
    sample.traits <-
      llply(random.species,
            function(area) {
              indices <- CrossCheck(x = rownames(trait.matrix),
                                    y = area,
                                    presence = TRUE,
                                    value = FALSE
                                    )
              trait.matrix[indices, ]
            }
            )
    # Step 2: Identify communities with < 4 unique trait combinations
    sample.counts <-
      llply(sample.traits,
            function(area) {
              count(area) %>%
                nrow()
            }
            ) %>%
      simplify2array()
    resample <- names(sample.counts)[sample.counts < 4]
    # Step 3: evaluate
    if (identical(length(resample), as.integer(0))) {
      do.sample <- FALSE
    } else {
      do.sample <- TRUE
    }
    sampling.iterations %<>% + 1
    if (sampling.iterations > 99) {
      stop("Unable to satisfy criterion S > T within 99 sampling attempts")
    }
  }
  if (verbose) {
    cat("Null model sampling required", sampling.iterations, "sampling attempts\n")
  }
  random.species
}

