#####################################################################
# Palm FD project: custom functions for null model using weighted ADF
#####################################################################
#
# Functions for (1) defining the assemblage dispersion field (ADF) for each
# community, (2) weighting the ADF by number of shared species with the focal
# community, (3) sample null communities using the weighted ADF species pool,
# and (4) run the null model for FD indices.
#
# Based on the probabilistic (weighted) ADF concept as described in
# Lessard et al (2012) and Carstensen et al 2013.
# NOTE: the original code by Lessard et al 2012 did NOT exclude the focal
# community from the ADF. We argue that the focal community SHOULD be excluded,
# because a site always shares the most species with itself and comparing a site
# with itself makes the null model z-score (SES) tend towards 0.
#
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
#   ADF
#   WeighADF
#   SampleADF
#   NullADF


library(magrittr)
library(plyr)

source(file = "functions/functional_diversity_functions.R")


ADF <- function(communities, pres.abs.matrix) {
# For a given set of communities, define the Assemblage Dispersion Field (ADF)
#
# Args:
#   communities: character vector giving the names of communities for which to
#                find the ADF.
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be tdwg3 labels, and colnames should be species names.
# Returns: 
#   A list giving, for each community, the names of the communities forming the
#   ADF of the focal community, EXCLUDING the focal community itself.
#   If exclusion results in an ADF of zero communities, an error is returned.
  adf <- vector("list", length = length(communities))
  names(adf) <- communities

  for (community in communities) {
    # 1. column indices of species in focal community
    col.indices <- 
      pres.abs.matrix %>%
      { .[match(community, rownames(.)), ] == 1 } %>%
      which()
    # 2. row indices of communities with 1 or more of these species
    row.indices <- which(rowSums(pres.abs.matrix[, col.indices, drop = FALSE]) > 0)
    # 3. Add corresponding community names to adf.
    adf[[match(community, communities)]] <-
      rownames(pres.abs.matrix)[row.indices]
  }
  # 4. Remove focal community from the adf
  # UPDATE: better to not exclude.
  # for (i in seq_along(adf)) {
  #   adf[[i]] %<>% .[!. %in% names(adf)[i]]
  # }

  # Integrity check
  adf.lengths <- unlist(llply(adf, length))
  if (any(adf.lengths < 1)) {
    names <- as.character(communities)[which(adf.lengths < 1)]
    stop(paste("Communities {",
               paste(names, collapse = " "),
               "} have an ADF of zero communities."
               )
         )
  }
  adf
}


WeighADF <- function(adf, pres.abs.matrix) {
# For a given list of ADFs, find the species in each ADF and their probability
# of being sampled for a WEIGHTED ADF null community.
#
# Args:
#   adf: List with for each community (element name), a character vector giving the
#        community names that make up the ADF. (i.e. the output of function 'ADF').
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be community names, and colnames should be species names.
# Returns:
#   List with, for each community in adf, a dataframe giving the species in the
#   focal communities' ADF and their overall sampling probability.
  output <- vector("list", length = length(adf))
  names(output) <- names(adf)
  for (i in seq_along(adf)) {
    # find species in ADF
    rows <- which(rownames(pres.abs.matrix) %in% adf[[i]])
    species <- 
      colnames(pres.abs.matrix)[colSums(pres.abs.matrix[rows, , drop = FALSE]) > 0]

    # find sample probabilities for ADF communities
    focal.species <-
      pres.abs.matrix[match(names(adf)[i], rownames(pres.abs.matrix)),
                       ,
                      drop = FALSE
                      ] %>%
      { colnames(.)[. > 0] }
    community.species <-
      alply(adf[[i]],
            .margins = 1,
            function(x) {
              pres.abs.matrix[match(x, rownames(pres.abs.matrix)),
                               ,
                              drop = FALSE
                              ] %>%
                { colnames(.)[. > 0] }
            }
            )
    names(community.species) <- adf[[i]]
    community.weights <-
      llply(community.species, function(x) { sum(x %in% focal.species) } ) %>%
      unlist()
    community.chances <- community.weights / sum(community.weights)

    # find sample probabilities for species in ADF species pool
    # 1. init matrix of species vs communities
    mat <- matrix(nrow = length(species),
                  ncol = length(adf[[i]]),
                  dimnames = list(row = species, col = adf[[i]]),
                  data = 0
                  )
    # 2. Populate matrix with probability a species is sampled IF a given community
    # is picked
    for (com in seq_len(ncol(mat))) {
      row.ind <- which(rownames(mat) %in% community.species[[com]])
      mat[row.ind, com] <- 1 / length(community.species[[com]])
      # 3. multiply these probabilities with the community.chances
      mat[, com] %<>% { . * community.chances[com] }
    }
    # 4. find species chances
    species.chances <- rowSums(mat)

    # parse output data
    output[[i]] <- data.frame(species = species,
                              weights = species.chances,
                              row.names = NULL
                              )
  }
  output
}


SampleADF <- function(species.pool, n, trait.matrix) {
# Draw a null community from a given species pool, without replacement.
#
# Args:
#   species.pool: a dataframe with two columns: a column 'species' with species
#                    names, and a column 'weights' with sampling weights for
#                    those species.
#   n: number of species to draw for the null community
#   trait.matrix: A matrix with observed numerical trait values in columns,
#                 without NAs, and rownames with species names. Used to
#                 verify the S > T criterion.
# Returns:
#   A character vector with the species names in the null community
  do.sample <- TRUE
  sampling.iterations <- 0
  while (do.sample) {
    null.community <- sample(x = species.pool[, "species"],
                             size = n,
                             replace = FALSE,
                             prob = species.pool[, "weights"]
                             ) %>%
                             as.character()
    # Hold your horses: check if the community has at least 4 functionally
    # unique species, as required for FRic calculation.
    # 1. subset trait matrix to sampled species
    sample.traits <-
      trait.matrix[rownames(trait.matrix) %in% null.community, ] %>%
      as.data.frame()                 
    # 2. Count unique trait combinations
    sample.count <- nrow(count(sample.traits))
    # 3. evaluate
    if (sample.count < 4) {
      do.sample <- TRUE
    }
    # 4. Similarly, check if each single trait has at least 2 unique values
    sample.trait.counts <-
      colwise(function(x) { length(unique(x)) } ) (sample.traits)
    # 5. Evaluate
    if (any(sample.trait.counts < 2)) {
      do.sample <- TRUE
    } else {
      do.sample <- FALSE
    }
    # wrap up
    sampling.iterations %<>% + 1
    if (sampling.iterations > 99) {
      stop("Unable to satisfy criterion S > T within 99 sampling attempts")
    }
  }
  null.community
}


NullADF <- function(trait.matrix,
                    pres.abs.matrix,
                    species.pools,
                    process.dir,
                    iterations = 1000,
                    mc.cores = getOption("mc.cores", 2L),
                    subset = FALSE,
                    verbose = TRUE,
                    single.traits = TRUE,
                    fast = FALSE
                    ) {
# Generate a null model for the given dataset, returning functional richness
# (FRic) and functional dispersion (FDis) computed for all traits and for each
# individual trait.
# This function uses null communities based on weighted species pools, derived
# from e.g. a weighted Assemblage Dispersion Field (ADF).
#
# Args:
#   trait.matrix: A matrix with observed numerical trait values in columns,
#                 without NAs, and rownames with species names.
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be tdwg3 labels, and colnames should be species names.
#                    The colnames should match the rownames of the trait.matrix!
#   species.pools: A list of dataframes giving the species pools for each
#                  community. Each dataframe should include two columns: a column
#                  'species' with species names, and a column 'weights' with
#                  sampling weights for those species.
#   process.dir: character string giving the name of a processing directory
#                to store intermediate files (with a trailing /). Will be created
#                if it doesn't exist.
#   iterations: number of random samples for which to calculate the FD.
#   mc.cores: number of cores to use for parallel processing of iterations.
#   subset: Logical indicating whether to subset the data to only the first 10
#           communities. This speeds up code execution, for testing purposes.
#   verbose: output extended info about progress?
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
      # Extract species richness for each community
      # -----------------------------------
      richness <- data.frame(community = rownames(pres.abs.matrix),
                             richness  = rowSums(pres.abs.matrix)
                             )

      # Randomly sampling null model communities
      # ----------------------------------------
      if (verbose) {
        cat("Randomly sampling null model communities from species pools...\n")
      }
      # We need a new presence/absence matrix of randomly sampled species for our
      # communities.
      # 1. Draw random sample
      nm.species <- vector("list", length = length(species.pools))
      names(nm.species) <- names(species.pools)
      for (i in seq_along(species.pools)) {
        nm.species[[i]] <- SampleADF(species.pools[[i]],
                                     richness[i, "richness"],
                                     trait.matrix
                                     )
      }
      # 2. Convert to presence/absence matrix
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

