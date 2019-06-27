##############################################################
# Palm FD project: calculation of Functional Diversity indices
##############################################################
#
# For each of the 100 gapfilled datasets, calculate functional richness (FRic)
# and functional dispersion (FRic).
# FD is calculated (1) using all three traits, and (2) for each trait individually.
# The community mean (CWM) trait values are also calculated.
# For comparison, FD is also calculated for the genus-mean gapfilled dataset.
#
# Calculation of multidimensional FD is done using trait coordinates in PCoA-space.
# Here, the PcoA coordinates are pre-computed for all datasets and saved to disk.
# These can then be re-used downstream to save processing time, particularly in
# null model procedures.


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")

# Subset FD input for quick code testing?
subset <- FALSE
# Enable verbose reporting in the FD calculation?
verbose <- TRUE
# Directory for saved FD outputs (with trailing slash)
fd.dir <- "output/observed_FD/"
# Directory for saved trait matrices (with trailing slash)
trait.dir <- "output/trait_matrices/"
# Number of cores to use for parallel processing. Default is 80% of available cores.
num.cores <- 
  if (!is.na(detectCores())) {
    floor(detectCores() * 0.95)
  } else {
    1
    warning("unable to detect cores. Parallel processing is NOT used!")
    cat("unable to detect cores. Parallel processing is NOT used!")
  }


######################################################
# Data preparation: trait matrices and pres/abs matrix
######################################################
cat("Preparing data...\n")

# Trait matrices
# --------------
# The matrices need 'species labels'. I'm assuming that means rownames.
# First, load the data:
traits.gapfilled <- vector("list", length = 100)
for (i in seq_along(traits.gapfilled)) {
  traits.gapfilled[[i]] <-
    read.csv(file = paste0("output/stochastic_gapfilled/",
                           "stochastic_gapfilled_traits_",
                           i,
                           ".csv"
                           )
             )
}
trait.names <- c("stem.height", "blade.length", "fruit.length")

# Second, convert the dataframes to properly formatted matrices
ToMatrix <- function(x) {
# Function to transform trait dataset x into a ready-to-use trait matrix
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
traits.gapfilled <- llply(traits.gapfilled, ToMatrix)

# Same for the genus-mean filled dataset
traits.mean <- read.csv(file = "output/traits_filled_genus_mean.csv")
traits.mean <- ToMatrix(traits.mean)
# This dataset has slightly more species. To keep the comparison honest,
# subset traits.mean to the species in traits.gapfilled
traits.mean %<>% .[rownames(.) %in% rownames(traits.gapfilled[[1]]), ]


# presence/absence matrix
# -----------------------
# Rows should be TDWG3 units and columns should be species.
# We only need one pres/abs matrix, as its the same for all 100 gapfilled datasets.

palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")
trait.matrix <- traits.gapfilled[[1]]

# Subset palm.dist to the species in trait.matrix
species.shared <- CrossCheck(x = unique(palm.dist[, "SpecName"]),
                             y = rownames(trait.matrix),
                             presence = TRUE,
                             value = TRUE
                             )
indices <- CrossCheck(x = palm.dist[, "SpecName"],
                      y = species.shared,
                      presence = TRUE,
                      value = FALSE
                      )
palm.dist %<>% .[indices, ]
palm.dist %<>% droplevels()

# Verify: do species in trait.matrix and palm.dist match?
if (!identical(length(rownames(trait.matrix)),
               length(unique(palm.dist[, "SpecName"]))
               )
    ) {
  stop("trait.matrix and palm.dist have differing number of species")
}
dist.species <-
  unique(palm.dist[, "SpecName"]) %>%
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
# As a result of these subsets, some species no longer occur in any included
# tdwg3 unit. 
# Remove these species from the pres/abs matrix, as well as from the trait 
# matrices (a requirement for FD::dbFD)
orphaned.species <-
  which(colSums(pres.abs.matrix) == 0) %>%
  colnames(pres.abs.matrix)[.]
pres.abs.matrix %<>% .[, -which(colSums(.) == 0)]
indices <- CrossCheck(x = rownames(trait.matrix),
                      y = orphaned.species,
                      presence = TRUE,
                      value = FALSE
                      )
traits.gapfilled %<>% llply(., function(x) { x[-indices, ] } )
traits.mean %<>% .[-indices, ]
# NOTE: this subsetting excludes 20 species.
# The remaining dataset is 2539 species in 131 botanical countries.


# Standardize trait matrices
# --------------------------
# All three traits should be standardized to mean 0 and unit variance.
traits.mean.std <- apply(traits.mean, 2, scale, center = TRUE, scale = TRUE)
rownames(traits.mean.std) <- rownames(traits.mean)

traits.gapfilled.std <-
  llply(traits.gapfilled,
        function(x) {
          output <- apply(x, 2, scale, center = TRUE, scale = TRUE)
          rownames(output) <- rownames(x)
          output
        }
        )


# Save all resulting matrices
# ---------------------------
# We will need them later for the null models.
write.csv(pres.abs.matrix,
          file = "output/palm_tdwg3_pres_abs_gapfilled.csv",
          eol = "\r\n",
          row.names = TRUE
          )

if (!dir.exists(trait.dir)) {
  cat("creating directory:", trait.dir, "\n")
  dir.create(trait.dir)
}

write.csv(traits.mean,
          file = paste0(trait.dir, "palm_trait_matrix_genus_mean.csv"),
          eol = "\r\n",
          row.names = TRUE
          )

write.csv(traits.mean.std,
          file = paste0(trait.dir, "palm_trait_matrix_genus_mean_standardized.csv"),
          eol = "\r\n",
          row.names = TRUE
          )

for (i in seq_along(traits.gapfilled)) {
  write.csv(traits.gapfilled[[i]],
            file = paste0(trait.dir, 
                          "palm_trait_matrix_filled_",
                          i,
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = TRUE
            )
}

for (i in seq_along(traits.gapfilled.std)) {
  write.csv(traits.gapfilled.std[[i]],
            file = paste0(trait.dir, 
                          "palm_trait_matrix_filled_standardized_",
                          i,
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = TRUE
            )
}


########################################
# Precompute PcoA-traits for null models
########################################
  cat("Precomputing pcoa traits...\n")

if (!dir.exists(fd.dir)) {
  cat("creating directory:", fd.dir, "\n")
  dir.create(fd.dir)
}

# For genus mean data
# -------------------
if (!file.exists(paste0(fd.dir, "pcoa_traits_mean_fruit.length.csv"))) {
  # for all traits
  pcoa.traits.mean <-
    dist(traits.mean.std) %>%
    dudi.pco(., scannf = FALSE, full = TRUE) %>%
    .$li
  write.csv(pcoa.traits.mean,
            file = paste0(fd.dir, "pcoa_traits_mean.csv"),
            eol = "\r\n",
            row.names = TRUE
            )
  # for single traits
  for (trait in colnames(traits.mean.std)) {
    subset.matrix <- traits.mean.std[, trait, drop = FALSE]
    pcoa.single <-
      dist(subset.matrix) %>%
      dudi.pco(., scannf = FALSE, full = TRUE) %>%
      .$li
    write.csv(pcoa.single,
              file = paste0(fd.dir, "pcoa_traits_mean_", trait, ".csv"),
              eol = "\r\n",
              row.names = TRUE
              )
  }
}

# For stochastic gapfilled data
# -----------------------------
if (!file.exists(paste0(fd.dir, "pcoa_traits_gapfilled_100_fruit.length.csv"))) {
  cluster <- makeCluster(num.cores)
  sample.list <-
    seq_along(traits.gapfilled.std) %>%
    as.array() %>%
    t() %>%
    as.list()
  clusterExport(cluster, ls())
  clusterEvalQ(cluster,
               {
                 library(magrittr)
                 library(plyr)
                 library(FD)
               }
               )
  parLapply(cluster,
            sample.list,
            function(x) {
              # For all traits
              pcoa.stochastic <-
                dist(traits.gapfilled.std[[x]]) %>%
                dudi.pco(., scannf = FALSE, full = TRUE) %>%
                .$li
              write.csv(pcoa.stochastic,
                        file = paste0(fd.dir,
                                      "pcoa_traits_gapfilled_",
                                      x,
                                      ".csv"
                                      ),
                        eol = "\r\n",
                        row.names = TRUE
                        )
              # For single traits
              for (trait in colnames(traits.gapfilled.std[[1]])) {
                subset.matrix <- traits.gapfilled.std[[1]] [, trait, drop = FALSE]
                pcoa.single <-
                  dist(subset.matrix) %>%
                  dudi.pco(., scannf = FALSE, full = TRUE) %>%
                  .$li
                write.csv(pcoa.single,
                          file = paste0(fd.dir,
                                        "pcoa_traits_gapfilled_",
                                        x,
                                        "_",
                                        trait,
                                        ".csv"
                                        ),
                          eol = "\r\n",
                          row.names = TRUE
                          )
              }
              return (0)
            }
            )
  stopCluster(cluster)
}


##################################
# Calculating Functional Diversity
##################################
cat("Calculating Functional Diversity indices... (this may take a while)\n")

SampleFD <- function(trait.matrix, pres.abs.matrix, id, pcoa.traits,
                     pcoa.traits.single) {
# Function to calculate all FD output for a single gapfilled trait matrix.
#
# Args:
#   trait.matrix: the trait matrix version to calculate FD for
#   pres.abs.matrix: the presence/absence matrix
#   id: a unique identifier, used for naming the produced files.
#   pcoa.traits: matrix giving pcoa.traits
#   pcoa.traits.single: list of matrices giving pcoa.traits for each single trait
#
# Returns: nothing. The function saves one file with the FD indices and one file
#          with the community mean trait values.

  # Calculate FD indices
  output.all <- RunFD(trait.matrix,
                      pres.abs.matrix,
                      subset = subset,
                      verbose = verbose,
                      fast = TRUE,
                      pcoa.traits = pcoa.traits
                      )
  output.single <- SingleFD(trait.matrix,
                            pres.abs.matrix,
                            subset = subset,
                            verbose = verbose,
                            fast = TRUE,
                            pcoa.traits.single = pcoa.traits.single
                            )
  output.combined <- c(list(all.traits = output.all), output.single)

  # Structure of FD output list:
  # $nbsp      - number of species in each community.
  # $sing.sp   - number of species with a unique trait combination in each
  #              community.
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

  # Parse and save FRic and FDis results
  # Best to combine all FD data into a single dataframe / file.
  fd.indices <- data.frame(FRic.all.traits   = output.combined[[1]]$FRic,
                           FRic.stem.height  = output.combined[[2]]$FRic,
                           FRic.blade.length = output.combined[[3]]$FRic,
                           FRic.fruit.length = output.combined[[4]]$FRic,
                           FDis.all.traits   = output.combined[[1]]$FDis,
                           FDis.stem.height  = output.combined[[2]]$FDis,
                           FDis.blade.length = output.combined[[3]]$FDis,
                           FDis.fruit.length = output.combined[[4]]$FDis,
                           row.names = names(output.combined[[1]]$FRic)
                           )
  write.csv(fd.indices,
            file = paste0(fd.dir, "FD_observed_", id, ".csv"),
            eol = "\r\n",
            row.names = TRUE
            )

  return (0)
}

# Apply SampleFD function to all 100 gapfilled datasets
# -----------------------------------------------------
# Using parallel processing of course.
# Initiate cluster:
cluster <- makeCluster(num.cores)

# Initiate list of samples to run
sample.list <-
  seq_along(traits.gapfilled) %>%
  as.array() %>%
  t() %>%
  as.list()

# Export environment to cluster
clusterExport(cluster, ls())
clusterEvalQ(cluster,
             {
               library(magrittr)
               library(plyr)
               library(FD)
             }
             )

# Parallel apply SampleFD to each STANDARDIZED sample gapfilled trait dataset
parLapply(cluster,
          sample.list,
          function(x) {
            pcoa.traits <- 
              read.csv(file = paste0(fd.dir,
                                     "pcoa_traits_gapfilled_",
                                     x,
                                     ".csv"
                                     ),
                       header = TRUE,
                       row.names = 1,
                       check.names = FALSE
                       )
            pcoa.traits.single <-
              vector("list", length = ncol(traits.gapfilled.std[[x]]))
            names(pcoa.traits.single) <- colnames(traits.gapfilled.std[[x]])
            for (i in seq_along(pcoa.traits.single)) {
              pcoa.traits.single[[i]] <-
                read.csv(file = paste0(fd.dir,
                                       "pcoa_traits_gapfilled_",
                                       x,
                                       "_",
                                       names(pcoa.traits.single)[i],
                                       ".csv"
                                       ),
                          header = TRUE,
                          row.names = 1
                          )
            }
            SampleFD(trait.matrix = traits.gapfilled.std[[x]],
                     pres.abs.matrix = pres.abs.matrix,
                     id = x,
                     pcoa.traits = pcoa.traits,
                     pcoa.traits.single = pcoa.traits.single
                     )
          }
          )

# Gracefully end cluster
stopCluster(cluster)

# SampleFD for genus-mean filled data
# -----------------------------------
pcoa.traits.mean <-
  read.csv(file = paste0(fd.dir, "pcoa_traits_mean.csv"),
           header = TRUE,
           row.names = 1,
           check.names = FALSE
           )
pcoa.traits.single <-
  vector("list", length = ncol(traits.mean.std))
names(pcoa.traits.single) <- colnames(traits.mean.std)
for (i in seq_along(pcoa.traits.single)) {
  pcoa.traits.single[[i]] <-
    read.csv(file = paste0(fd.dir,
                           "pcoa_traits_mean_",
                           names(pcoa.traits.single)[i],
                           ".csv"
                           ),
              header = TRUE,
              row.names = 1
              )
}
SampleFD(trait.matrix = traits.mean.std,
         pres.abs.matrix = pres.abs.matrix,
         id = "genus_mean",
         pcoa.traits = pcoa.traits.mean,
         pcoa.traits.single = pcoa.traits.single
         )

# Summarize results of the 100 stochastic gapfilled datasets
# ----------------------------------------------------------
StochasticMeans(stat = "",
                samples = 100,
                output.dir = fd.dir,
                header = "FD_observed_"
                )


#######################
# Community trait means
#######################
cat("Calculating community trait means...\n")
# In addition to FD, it is useful to have community average trait values.
# Based on the log-transformed, non-standardized trait matrices.

# Calculate the trait means with a function
# -----------------------------------------
TraitMean <- function(trait.matrix, pres.abs.matrix) {
  mat <- matrix(nrow = nrow(pres.abs.matrix),
                ncol = ncol(trait.matrix),
                dimnames = list(row = rownames(pres.abs.matrix),
                                col = colnames(trait.matrix)
                                ),
                data = NA
                )
  for (community in seq_len(nrow(pres.abs.matrix))) {
    species <- colnames(pres.abs.matrix)[pres.abs.matrix[community, ] > 0]
    traits.subset <- trait.matrix[rownames(trait.matrix) %in% species, ]
    mat[community, ] <- apply(traits.subset, 2, mean)
  }
  mat
}

traitmeans.mean <- TraitMean(traits.mean, pres.abs.matrix)
traitmeans.stochastic <-
  llply(traits.gapfilled, TraitMean, pres.abs.matrix = pres.abs.matrix)

# Save results
# ------------
write.csv(traitmeans.mean,
          file = paste0(fd.dir, "community_trait_means_genus_mean.csv"),
          eol = "\r\n",
          row.names = TRUE
          )

for (id in seq_along(traitmeans.stochastic)) {
  write.csv(traitmeans.stochastic[[id]],
            file = paste0(fd.dir, "community_trait_means_", id, ".csv"),
            eol = "\r\n",
            row.names = TRUE
            )
}

# Summarize trait means of the 100 stochastic gapfilled datasets
# --------------------------------------------------------------
StochasticMeans(stat = "",
                samples = 100,
                output.dir = fd.dir,
                header = "community_trait_means_"
                )


cat("Done.\n")

