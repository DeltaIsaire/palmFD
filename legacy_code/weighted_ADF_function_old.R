library(magrittr)
library(plyr)

WeightedADF <- function(communities, pres.abs.matrix) {
# For a given set of communities, define the Assemblage Diversity Field (ADF),
# Weigh it by the number of shared species with the focal community,
# Then sample a null community from the weighted ADF.
#
# Args:
#   communities: character vector giving the names of communities to generate
#                null communities for
#   pres.abs.matrix: A presence/absence matrix giving the presence (1) and
#                    absence (0) of palm species in tdwg3 units. Rownames should
#                    be tdwg3 labels, and colnames should be species names.
# Returns: 
#   A list giving the randomly sampled species names for each community

  # Initialize output list
  null.species <- vector("list", length = length(communities))
  names(null.species) <- communities

  for (community in communities) {
    # Define ADF for the community
    # ----------------------------
    # We can use the presence/absence matrix for this.
    # 1. column indices of species in focal community
    col.indices <- 
      pres.abs.matrix %>%
      { .[match(community, rownames(.)), ] == 1 } %>%
      which()

    # 2. row indices of communities with 1 or more of these species
    row.indices <- which(rowSums(pres.abs.matrix[, col.indices]) > 0)
    if (any(length(row.indices) < 2)) {
      stop(paste("All species in", community, "are endemic to this community"))
    }

    # 3. build a list of these these communities and their species
    adf <- vector("list", length = length(row.indices))
    names(adf) <- rownames(pres.abs.matrix)[row.indices]
    for (i in seq_along(adf)) {
      adf[[i]] <-
        pres.abs.matrix[match(names(adf)[i], rownames(pres.abs.matrix)), ] %>%
        { names(.)[. == 1] }
    }

    # 4. Get the sampling probability for each community in the ADF
    # Weight by number of shared species, then convert weights to probabilities
    weights <-
      llply(adf, function(x) { length(which(x %in% adf[[community]])) } ) %>%
      unlist()
    chances <- weights / sum(weights)
    
    # Construct null community from the ADF and chances
    # -------------------------------------------------
    # Initialize output: character of length = richness of the focal community
    null.community <- character(length = length(adf[[community]]))
    for (i in seq_along(null.community)) {
      do.sample <- TRUE
      while (do.sample) {
        # Select community by index
        index <- sample(x = sample.int(length(adf)), size = 1, prob = chances)
        # Select species from this community
        sample <- sample(x = adf[[index]], size = 1)
        if (!any(sample %in% null.community)) {
          null.community[i] <- sample
          do.sample <- FALSE
        }
      }
    }
    null.species[[match(community, names(null.species))]] <- null.community
  }
  null.species
}

