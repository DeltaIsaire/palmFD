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
# ResampleTest
# StochasticMeans


library(plyr)
library(magrittr)
library(FD)
library(reshape2)
library(parallel)

source(file = "functions/base_functions.R")
source(file = "functions/faster_FD_function.R")


RunFD <- function(trait.matrix,
                  pres.abs.matrix,
                  subset = FALSE,
                  verbose = TRUE,
                  fast = FALSE
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
#  fast: use custom faster version of FD::dbFD, tailored to our data?
#
# Returns:
#   A list with the output of function 'FD:dbFD' applied to the provided dataset;
#   or if fast = TRUE, a list of length 2 with FRic and FDis.
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
  if (fast) {
    output.all <- FastFD(trait.matrix = traits, pres.abs.matrix = pres.abs)
  } else {
    output.all <- suppressWarnings(dbFD(x = traits,
                                        a = pres.abs,
                                        w.abun = FALSE,
                                        stand.x = FALSE,  # Do NOT use this!
                                        corr = "cailliez",
                                        calc.FRic = TRUE,
                                        m = "max",
                                        stand.FRic = FALSE,
                                        calc.CWM = TRUE,
                                        calc.FDiv = FALSE,
                                        messages = verbose
                                        )
                                   )
  }
  output.all
}


SingleFD <- function(trait.matrix,
                     pres.abs.matrix,
                     subset = FALSE,
                     verbose = TRUE,
                     fast = FALSE
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
#  fast: use custom faster version of FD::dbFD, tailored to our data?
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
                              verbose = verbose,
                              fast = fast
                              )
  }
  names(output.list) <- colnames(trait.matrix)
  output.list
}


RandomSpecies <- function(communities, richness, species, trait.matrix,
                          verbose = TRUE, groups = TRUE) {
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
#  groups: does 'communities' contain a list of groups? If FALSE, it is assumed
#          that each element of 'communities' corresponds to a community, and
#          contains a character vector with community names representing
#          the species pool. Element names should be community names.
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
    if (isTRUE(groups)) {
      random.species <- rep(list(NULL), length = length(unlist(areas)))
      i <- 1
      for (group in seq_along(areas)) {
        for (area in seq_along(areas[[group]])) {
          area.richness <-
            richness %>% {
              .[which(.[1] == areas[[group]][area]), 2]
            }
          indices <- sample(size = area.richness,
                            replace = FALSE,
                            x = sample.int(length(list.species[[group]]))
                            )
          random.species[[i]] <- list.species[[group]][indices]
          names(random.species)[i] <- areas[[group]][area]
          i <- i + 1
        }
      }
    } else {
      random.species <- rep(list(NULL), length = length(areas))
      for (area in seq_along(areas)) {
        area.richness <-
          richness %>% {
            .[which(.[1] == names(areas)[area]), 2]
          }
        indices <- sample(size = area.richness,
                          replace = FALSE,
                          x = sample.int(length(list.species[[area]]))
                          )
        random.species[[area]] <- list.species[[area]][indices]
        names(random.species)[area] <- names(areas)[area]
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
    # Step 4: check for trait invariance
    sample.trait.counts <-
      llply(sample.traits,
            function(area) {
              colwise(function(x) { length(unique(x)) } ) (as.data.frame(area)) %>%
                min()
            }
            ) %>%
      simplify2array()
    # Step 5: evaluate
    if (any(sample.trait.counts < 2)) {
      do.sample <- TRUE
    }
    # wrap up
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
                      groups = list(group = rownames(pres.abs.matrix)),
                      process.dir,
                      iterations = 100,
                      mc.cores = getOption("mc.cores", 2L),
                      subset = FALSE,
                      verbose = TRUE,
                      random.groups = TRUE,
                      single.traits = TRUE,
                      fast = FALSE
                      ) {
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
#   mc.cores: number of cores to use for parallel processing of iterations.
#   subset: Logical indicating whether to subset the data to only the first 10
#           communities. This speeds up code execution, for testing purposes.
#   verbose: output extended info about progress?
#   random.groups: logical indicating whether 'groups' contains groups. Also
#                  passed to argument 'groups' of function 'RandomSpecies'
#   single.traits: Include calculation of FD indices for single traits?
#   fast: use custom faster version of FD::dbFD, tailored to our data?
#
# Returns:
#   A nested list with dataframes containing sampled FD indices for each
#   community
  if (verbose) {
    cat("Generating null model with", iterations, "samples\n")
  }
  if (!dir.exists(process.dir)) {
    dir.create(process.dir, recursive = TRUE)
  }

  # Initiate cluster:
  cluster <- makeCluster(mc.cores)
  # Initiate list of iterations to run
  iteration.list <-
    seq_len(iterations) %>%
    as.array() %>%
    t() %>%
    as.list()
  # Export environment to cluster
  clusterExport(cluster, varlist = ls(), envir = environment())
  clusterExport(cluster, varlist = ls(globalenv()), envir = environment())
  clusterEvalQ(cluster,
               {
                 library(plyr)
                 library(magrittr)
                 library(FD)
                 library(reshape2)
               }
               )
  # parallel apply nullmodel code to each iteration
  parLapply(cluster,
            iteration.list,
            function(iteration) {
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
      if (!isTRUE(random.groups)) {
        groups.subset %<>% .[names(.) %in% rownames(pres.abs.matrix)]
      }

      # Then, for each group, extract the present species
      group.species <-
        llply(groups.subset,
              function(group) {
                if (!identical(group, character(0))) {
                  indices <- CrossCheck(x = rownames(pres.abs.matrix),
                                                     y = group,
                                                     presence = TRUE,
                                                     value = FALSE
                                                     )
                  pres.abs.subset <- pres.abs.matrix[indices, , drop = FALSE]
                  species <-
                    pres.abs.subset %>%
                    { colnames(.)[which(colSums(.) > 0)] }
                } else {
                  species <- character(0)
                }
                species
              }
              )

      # Randomly sampling null model communities
      # ----------------------------------------
      if (verbose) {
        cat("Randomly sampling null model communities from species pools...\n")
      }
      # We need a new presence/absence matrix of randomly sampled species for our
      # communities.
      nm.species <- RandomSpecies(communities = groups.subset,
                                  richness = richness,
                                  species = group.species,
                                  trait.matrix = trait.matrix,
                                  verbose = verbose,
                                  groups = random.groups
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
      nm.pres.abs %<>% .[, order(colnames(.)), drop = FALSE]
      # Subset trait matrix to species in the null model
      indices <- CrossCheck(x = rownames(trait.matrix),
                            y = colnames(nm.pres.abs),
                            presence = TRUE,
                            value = FALSE
                            )
      trait.matrix.subset <- trait.matrix[indices, ]
      trait.matrix.subset %<>% .[order(rownames(.)), ]
  
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
                          verbose = verbose,
                          fast = fast
                          )
      # FD for single traits
      if (single.traits) {
        if (verbose) {
          cat("(2) FD for single traits:\n")
        }
        output.single <- SingleFD(trait.matrix.subset,
                                  nm.pres.abs,
                                  subset = subset,
                                  verbose = verbose,
                                  fast = fast
                                  )
        # gather results
        fd.indices <- c(list(all.traits = output.all), output.single)
      } else {
        fd.indices <- list(all.traits = output.all)
      }

      # Parse and save results
      # ----------------------
      if (verbose) {
        cat("Parsing and saving output...\n")
      }
      # Minimize clutter, so generate one file from a single dataframe
      trait.names <- colnames(trait.matrix)
      result <- data.frame(community          = names(fd.indices[[1]]$FRic),
                           all.traits.FRic    = fd.indices[[1]]$FRic,
                           all.traits.FDis    = fd.indices[[1]]$FDis
                           )
      if (single.traits) {
        for (trait in seq_along(trait.names)) {
          result[, paste0(trait.names[trait], ".FRic")] <-
            fd.indices[[trait + 1]]$FRic
          result[, paste0(trait.names[trait], ".FDis")] <-
            fd.indices[[trait + 1]]$FDis
        }
      }
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
  }
  )
  # Gracefully end cluster
  stopCluster(cluster)
  
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
  if (single.traits) {
    null.model <- list(FRic = rep(list(mat), (ncol(template) - 1) / 2),
                       FDis = rep(list(mat), (ncol(template) - 1) / 2)
                       )
    names(null.model$FRic) <- c("all.traits", colnames(trait.matrix))
    names(null.model$FDis) <- c("all.traits", colnames(trait.matrix))
  } else {
    null.model <- list(FRic = list(all.traits = mat),
                       FDis = list(all.traits = mat)
                       )
  }
  # And then fill the list with the data
  for (i in seq_len(iterations)) {
    data <- read.csv(file = paste0(process.dir,
                                   "FD_null_model_iteration_",
                                   i,
                                   ".csv"
                                   ),
                     header = TRUE
                     )
    null.model[[1]] [[1]] [i, ] <- data[, 2]  # all.traits.FRic
    null.model[[2]] [[1]] [i, ] <- data[, 3]  # all.traits.FDis
    if (single.traits) {
      for (trait in seq_along(colnames(trait.matrix))) {
        null.model[[1]] [[trait + 1]] [i, ] <- data[, trait * 2 + 2]  # trait.FRic
        null.model[[2]] [[trait + 1]] [i, ] <- data[, trait * 2 + 3]  # trait.FDis
      }
    }
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
#   A dataframe with for each row (community) the z.score, nullmodel mean
#   and nullmodel sd.
  result <- data.frame(z.score     = numeric(length = length(x)),
                       sample.mean = numeric(length = length(x)),
                       sample.sd   = numeric(length = length(x)),
                       row.names   = names(x)
                       )
  result[, "z.score"] <- x
  for (i in seq_along(x)) {
    column <- match(names(x)[i], colnames(y))
    mean <- mean(y[, column])
    sd <- sd(y[, column])
    result[, "z.score"][i] <- (x[i] - mean) / sd
    result[, "sample.mean"][i] <- mean
    result[, "sample.sd"][i] <- sd
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
  result <- rep(list(raw.fd), 3)
  names(result) <- c("z.scores", "means", "sds")
  for (i in seq_along(raw.fd)) {
    model <- ZScore(x = structure(raw.fd[, i], names = rownames(raw.fd)),
                    y = null.list[[match(colnames(raw.fd)[i],
                                         names(null.list)
                                         )
                                   ]]
                    )
    result[["z.scores"]] [, i] <- model[, "z.score"]
    result[["means"]] [, i] <- model[, "sample.mean"]
    result[["sds"]] [, i] <- model[, "sample.sd"]
  }
  result
}


ResampleTest <- function(pres.abs.matrix, trait.matrix, groups,
                         random.groups = TRUE, verbose = TRUE, samples) {
# Test how often resampling the random species sample is necessary
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
#   random.groups: logical indicating whether 'groups' contains groups. Also
#                  passed to argument 'groups' of function 'RandomSpecies'
#   verbose: print info about progress?
#   samples: integer indicating how many times to run the species sampling routine
#
# Returns:
#
#
  # Define grouping of tdwg3 units
  if (verbose) {
    cat("Parsing community grouping...\n")
  }
  groups.subset <-
  llply(groups,
        function(group) {
          CrossCheck(x = rownames(pres.abs.matrix),
                     y = group,
                     presence = TRUE,
                     value = TRUE
                     )
        }
        )

  richness <- data.frame(community = rownames(pres.abs.matrix),
                         richness  = rowSums(pres.abs.matrix)
                         )

  if (!isTRUE(random.groups)) {
    groups.subset %<>% .[names(.) %in% rownames(pres.abs.matrix)]
  }

  # Then, for each group, extract the present species
  if (verbose) {
    cat("Building species pool for each community...\n")
  }
  group.species <-
    llply(groups.subset,
          function(group) {
            if (!identical(group, character(0))) {
              indices <- CrossCheck(x = rownames(pres.abs.matrix),
                                                 y = group,
                                                 presence = TRUE,
                                                 value = FALSE
                                                 )
              pres.abs.subset <- pres.abs.matrix[indices, , drop = FALSE]
              species <-
                pres.abs.subset %>%
                { colnames(.)[which(colSums(.) > 0)] }
            } else {
              species <- character(0)
            }
            species
          }
          )

  # Run sampling routine (edited copy of RandomSpecies)
  if (verbose) {
    cat("Preparing input for species sampling routine...\n")
  }
  communities <- groups.subset
  richness <- richness
  species <- group.species
  trait.matrix <- trait.matrix
  verbose <- verbose
  groups <- random.groups

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

  # Run the sampling routine n times
  if (verbose) {
    cat("Running species sampling routine", samples, "times\n")
  }
  # Initialize output matrix
  # Initial data is all zeroes. In each instance where resampling would be
  # required, the value will be changed to 1.
  result <- matrix(nrow = length(unlist(areas)),
                   ncol = samples,
                   dimnames = list(row = unlist(areas),
                                   col = seq_len(samples)
                                   ),
                   data = 0
                   )
  for (sample in seq_len(samples)) {
    if (verbose) {
      cat("sample", sample, "\n")
    }
    if (isTRUE(groups)) {
      random.species <- rep(list(NULL), length = length(unlist(areas)))
      i <- 1
      for (group in seq_along(areas)) {
        for (area in seq_along(areas[[group]])) {
          area.richness <-
            richness %>% {
              .[which(.[1] == areas[[group]][area]), 2]
            }
          indices <- sample(size = area.richness,
                            replace = FALSE,
                            x = sample.int(length(list.species[[group]]))
                            )
          random.species[[i]] <- list.species[[group]][indices]
          names(random.species)[i] <- areas[[group]][area]
          i <- i + 1
        }
      }
    } else {
      random.species <- rep(list(NULL), length = length(areas))
      for (area in seq_along(areas)) {
        area.richness <-
          richness %>% {
            .[which(.[1] == names(areas)[area]), 2]
          }
        indices <- sample(size = area.richness,
                          replace = FALSE,
                          x = sample.int(length(list.species[[area]]))
                          )
        random.species[[area]] <- list.species[[area]][indices]
        names(random.species)[area] <- names(areas)[area]
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
    resample.indices <- which(sample.counts < 4)
    # add results to output
    result[resample.indices, sample] <- 1
  }

  # Parse results: count frequency of resampling requirements
  if (verbose) {
    cat("Processing results...\n")
  }
  data.frame(community = rownames(result),
             resample.count = rowSums(result),
             resample.freq = rowSums(result) / samples,
             row.names = seq_len(nrow(result))
             )
}


StochasticMeans <- function(stat, samples, output.dir, header) {
# Function to assemble the mean over n null model output runs.
# For averaging the null model output of the 100 stochastic runs.
#
# Args:
#   stat: statistic to assemble the means for: a character string in
#         c("means", "sds", "z.scores")
#   samples: single integer giving the number of runs to average over.
#   output.dir: directory path (with a trailing slash) where the null model
#               output files are located.
#   header: character string giving the first part of the filename shared by all
#           output files.
# Returns:
#   A matrix with means. The matrix is also saved to the output.dir
  if (!identical(stat, "")) {
    if(!identical(substr(stat, 1, 1), "_")) {
      stat <- paste0("_", stat)
    }
  }

  # Load first sample as template, and construct an array
  template <- read.csv(file = paste0(output.dir,
                                     header,
                                     "1",
                                     stat,
                                     ".csv"
                                     ),
                       header = TRUE,
                       row.names = 1
                       )
  means.array <- array(data = NA,
                       dim = c(nrow(template),
                               ncol(template),
                               samples
                               ),
                       dimnames = list(row = rownames(template),
                                       col = colnames(template),
                                       sample = seq_len(samples)
                                       )
                       )
  # Fill array with the data
  for (i in seq_len(samples)) {
    means.array[, , i] <-
      read.csv(file = paste0(output.dir,
                             header,
                             i,
                             stat,
                             ".csv"
                             ),
               header = TRUE,
               row.names = 1
               ) %>%
      as.matrix()
  }
  # Collapse array into matrix by averaging over the samples
  mat <- means.array[, , 1]
  for (row in seq_len(nrow(mat))) {
    for (col in seq_len(ncol(mat))) {
      mat[row, col] <- mean(means.array[row, col, ])
    }
  }
  # Save result
  write.csv(mat,
            file = paste0(output.dir,
                          header,
                          "stochastic_mean_of",
                          stat,
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = TRUE
            )
  # return result
  mat
}

