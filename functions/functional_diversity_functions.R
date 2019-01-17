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
# NullModel
# ZScore
# NullTransform


library(plyr)
library(magrittr)
library(FD)
library(reshape2)
library(parallel)


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
# Returns:
#   A list with the output of function 'FD:dbFD' applied to the provided dataset.
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
# Returns:
#   A list of the same length as the number of traits in the trait matrix,
#   where each element is a list with the output of function 'FD:dbFD'
#   applied to the dataset subsetted to a single trait. 
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
      cat("Calculating FD for trait", colnames(trait.matrix)[i], "...\n")
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
          runif(n = area.richness * 100,
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


NullModel <- function(trait.matrix,
                      pres.abs.matrix,
                      groups = list(rownames(pres.abs.matrix)),
                      process.dir,
                      iterations = 100,
                      mc.cores = getOption("mc.cores", 2L) - 1,
                      subset,
                      verbose) {
# Generate a null model for the given dataset, returning functional richness
# (FRic) and functional dispersion (FDis) computed for all traits and for each
# individual trait.
#
# Args:
#   trait.matrix: A matrix with observed numerical trait values in columns,
#                 without NAs, and rownames with species names.
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be tdwg3 labels, and colnames should be species names.
#                    The colnames should match the rownames of the trait.matrix!
#   groups: A list of character vectors giving the community names belonging
#           to each group (list element). Community names must match the
#           rownames of the pres.abs.matrix. The list can be of length 1 (i.e.
#           only one group, also the default)
#   process.dir: character string giving the name of a processing directory
#                to store intermediate files (with a trailing /). Will be created
#                if it doesn't exist.
#   iterations: number of random samples for which to calculate the FD.
#   mc.cores: number of cores to use for parallel processing of iterations
#             (by passing the argument to 'parallel:mclapply')
#   subset: Logical indicating whether to subset the data to only the first 10
#           tdwg3 units. This speeds up code execution, for testing purposes.
#   verbose: output extended info about progress?
#
# Returns:
#   A nested list with 8 dataframes containing sampled FD indices for each
#   community.
  if (verbose) {
    cat("Generating null model with", iterations, "samples\n")
  }
  if (!dir.exists(process.dir)) {
    dir.create(process.dir, recursive = TRUE)
  }
  iteration.list <-
    seq_len(iterations) %>%
    as.array() %>%
    t() %>%
    as.list()
  mclapply(iteration.list,
           mc.cores = mc.cores,
           FUN = function(iteration) {
    if (!file.exists(paste0(process.dir,
                            "FD_null_model_iteration_",
                            iteration,
                            ".csv"
                            )
                     )
        ) {
      if (verbose) {
        cat("---------------------\n")
        cat("Iteration", iteration, "\n")
        cat("---------------------\n")
      }
      # Extract species list for each group
      # -----------------------------------
      if (verbose) {
        cat("Building species list for each group...\n")
      }
      # First, extract list of communities for each group that has palms,
      # based on the presence/absence matrix.
      richness <- data.frame(community = rownames(pres.abs.matrix),
                             richness  = rowSums(pres.abs.matrix)
                             )
      groups.subset <- {
        llply(groups,
              function(group) {
                CrossCheck(x = rownames(pres.abs.matrix),
                           y = group,
                           presence = TRUE,
                           value = TRUE
                           )
              }
              )
      }
      # Then, for each group, extract the present species
      group.species <- llply(groups.subset,
                             function(group) {
                               indices <- CrossCheck(x = rownames(pres.abs.matrix),
                                                     y = group,
                                                     presence = TRUE,
                                                     value = FALSE
                                                     )
                               pres.abs.subset <- pres.abs.matrix[indices, ]
                               species <-
                                 pres.abs.subset %>%
                                 { colnames(.)[which(colSums(.) > 0)] }
                               species
                             }
                             )
  
      # Randomly sampling null model communities
      # ----------------------------------------
      if (verbose) {
        cat("Randomly sampling null model communities from realm species pool...\n")
      }
      # We need a new presence/absence matrix of randomly sampled species for our
      # communities.
      nm.species <- RandomSpecies(communities = groups.subset,
                                  richness = richness,
                                  species = group.species,
                                  trait.matrix = trait.matrix,
                                  verbose = verbose
                                  )
      nm.species %<>% melt(., value.name = "species")
      # Melt is an awesome function but its output is a mess, so restructure:
      names(nm.species)[2] <- "community"
      nm.species[, 2] %<>% as.factor()
      nm.species %<>% .[, c(2, 1)]
      # Then transform:
      nm.pres.abs <-
        table(nm.species) %>%
        as.data.frame.matrix() %>%
        as.matrix()
      nm.pres.abs %<>% .[, order(colnames(.))]
      # Subset trait matrix to species in the null model
      indices <- CrossCheck(x = rownames(trait.matrix),
                            y = colnames(nm.pres.abs),
                            presence = TRUE,
                            value = FALSE
                            )
      trait.matrix.subset <- trait.matrix[indices, ]
  
      # Calculating Functional Diversity
      # --------------------------------
      if (verbose) {
        cat("Calculating Functional Diversity indices... (this may take a while)\n")
      }
      # FD using all traits
      if (verbose) {
        cat("(1) FD with all traits...\n")
      }
      output.all <- RunFD(trait.matrix.subset,
                          nm.pres.abs,
                          subset = subset,
                          verbose = verbose
                          )
      # FD for single traits
      if (verbose) {
        cat("(2) FD for single traits:\n")
      }
      output.single <- SingleFD(trait.matrix.subset,
                                nm.pres.abs,
                                subset = subset,
                                verbose = verbose
                                )
      # gather results
      fd.indices <- c(list(all.traits = output.all), output.single)
  
      # Parse and save results
      # ----------------------
      if (verbose) {
        cat("Parsing and saving output...\n")
      }
      # Minimize clutter, so generate one file from a single dataframe
      result <- data.frame(community          = names(fd.indices[[1]]$nbsp),
                           all.traits.FRic    = fd.indices[[1]]$FRic,
                           all.traits.FDis    = fd.indices[[1]]$FDis,
                           stem.height.FRic   = fd.indices[[2]]$FRic,
                           stem.height.FDis   = fd.indices[[2]]$FDis,
                           blade.length.FRic  = fd.indices[[3]]$FRic,
                           blade.length.FDis  = fd.indices[[3]]$FDis,
                           fruit.length.FRic  = fd.indices[[4]]$FRic,
                           fruit.length.FDis  = fd.indices[[4]]$FDis
                           )
      write.csv(result,
                file = paste0(process.dir,
                              "FD_null_model_iteration_",
                              iteration,
                              ".csv"
                              ),
                eol = "\r\n",
                row.names = FALSE
                )
    }
    return (0)
  })
  
  # With all iterations complete, gather up the results
  # ---------------------------------------------------
  if (verbose) {
    cat("---------------------\n")
    cat("Processing results...\n")
  }
  # First initialize a result list
  template <- read.csv(file = paste0(process.dir,
                                     "FD_null_model_iteration_1.csv"
                                     ),
                       header = TRUE
                       )
  mat <- matrix(ncol = nrow(template),
                nrow = iterations,
                dimnames = list(seq_len(iterations), template[, 1])
                )
  null.model <- list(FRic = rep(list(mat), 4),
                     FDis = rep(list(mat), 4)
                     )
  names <- c("all.traits", "stem.height", "blade.length", "fruit.length")
  names(null.model$FRic) <- names
  names(null.model$FDis) <- names
  # And then fill the list with the data
  for (i in seq_len(iterations)) {
    data <- read.csv(file = paste0(process.dir,
                                   "FD_null_model_iteration_",
                                   i,
                                   ".csv"
                                   ),
                     header = TRUE
                     )
    null.model[[1]] [[1]] [i, ] <- data[, 2]
    null.model[[1]] [[2]] [i, ] <- data[, 4]
    null.model[[1]] [[3]] [i, ] <- data[, 6]
    null.model[[1]] [[4]] [i, ] <- data[, 8]
    null.model[[2]] [[1]] [i, ] <- data[, 3]
    null.model[[2]] [[2]] [i, ] <- data[, 5]
    null.model[[2]] [[3]] [i, ] <- data[, 7]
    null.model[[2]] [[4]] [i, ] <- data[, 9]
  }
  null.model
}


ZScore <- function(x, y) {
# Given a vector of FD index values and a corresponding null model output,
# calculate z-scores for the FD index.
#
# Args:
#   x: A named numeric vector of FD index values, where names identify the
#      community
#   y: A dataframe with a column for each community, containing null model
#      FD index values corresponding to x. Column names should match names of
#      x (but do not need to be in the same order).
#
# Returns:
#   A vector with the same length and names as x, containing null model z-scores
#   for the FD index.
  result <- x
  for (i in seq_along(x)) {
    column <- match(names(x)[i], colnames(y))
    result[i] <- (x[i] - mean(y[, column])) / sd(y[, column])
  }
  result
}


NullTransform <- function(raw.fd, null.model) {
# Using the output from an FD null model, z-transform the raw FD indices
# using 'ZScore'.
#
# Args:
#   raw.fd: dataframe with raw FD indices. column names should match element names
#           of null.model
#   null.model: list of dataframes with null model output, see 'Zscore'.
#               Can be the direct output of 'NullModel'.
#
# Returns:
#   Dataframe with the same dimensions as raw.fd, containing the transformed
#   FD indices.

  # The default output of NullModel is a nested list, but we need an un-nested list.
  if (identical(class(null.model[[1]]), "list")) {
    null.list <- unlist(null.model, recursive = FALSE)
  } else {
    null.list <- null.model
  }
  # Apply ZScore to each combination of raw and null model data
  result <- raw.fd
  for (i in seq_along(raw.fd)) {
    result[, i] <- ZScore(x = structure(raw.fd[, i], names = rownames(raw.fd)),
                          y = null.list[[match(colnames(raw.fd)[i],
                                               names(null.list)
                                               )
                                         ]]
                          )
  }
  result
}


