#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we explore how to properly implement calculation of FD indices.
# For each of the 100 gapfilled datasets, calculate functional richness (FRic)
# and functional dispersion (FRic).
# FD is calculated (1) using all three traits, and (2) for each trait individually.
# The community mean (CWM) trait values are also calculated.
#
# For comparison, FD is also calculated for the genus-mean gapfilled dataset.
#
#
# Input files:
#   dir 'output/test/stochastic_gapfilled/' with stochastically filled
#     trait matrices
#   data/palms_in_tdwg3.csv
#   output/test/traits_filled_genus_mean.csv
# Generated output files:
#   output/test/test_palm_tdwg3_pres_abs_gapfilled.csv
#   < output dir specified below, with 100 FD files and 100 CWM files >


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")

# Subset FD input for quick code testing?
subset <- TRUE
# Enable verbose reporting in the FD calculation?
verbose <- TRUE
# Directory for saved FD outputs (with trailing slash)
fd.dir <- "output/test/observed_FD/"
# Directory for saved trait matrices (with trailing slash)
trait.dir <- "output/test/trait_matrices/"
# Number of cores to use for parallel processing. Default is 80% of available cores.
num.cores <- 
  if (!is.na(detectCores())) {
    floor(detectCores() * 0.8)
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
    read.csv(file = paste0("output/test/stochastic_gapfilled/",
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
traits.mean <- read.csv(file = "output/test/traits_filled_genus_mean.csv")
traits.mean <- ToMatrix(traits.mean)
# This dataset has slightly more species. To keep the comparison honest,
# subset traits.mean to the species in traits.gapfilled
traits.mean %<>% .[rownames(.) %in% rownames(traits.gapfilled[[1]]), ]


# -----------------------
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

# Save pres/abs matrix and all trait matrices.
# We will need them later for the null models.
write.csv(pres.abs.matrix,
          file = "output/test/test_palm_tdwg3_pres_abs_gapfilled.csv",
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


##################################
# Calculating Functional Diversity
##################################
cat("Calculating Functional Diversity indices... (this may take a while)\n")

if (!dir.exists(fd.dir)) {
  cat("creating directory:", fd.dir, "\n")
  dir.create(fd.dir)
}

SampleFD <- function(trait.matrix, pres.abs.matrix, id) {
# Function to calculate all FD output for a single gapfilled trait matrix.
#
# Args:
#   trait.matrix: the trait matrix version to calculate FD for
#   pres.abs.matrix: the presence/absence matrix
#   id: a unique identifier, used for naming the produced files.
#
# Returns: nothing. The function saves one file with the FD indices and one file
#          with the community mean trait values.

  # Calculate FD indices
  output.all <- RunFD(trait.matrix,
                      pres.abs.matrix,
                      subset = subset,
                      verbose = verbose,
                      fast = FALSE
                      )
  output.single <- SingleFD(trait.matrix,
                            pres.abs.matrix,
                            subset = subset,
                            verbose = verbose,
                            fast = FALSE
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
            file = paste0(fd.dir, "test_FD_observed_", id, ".csv"),
            eol = "\r\n",
            row.names = TRUE
            )
  # Community weighted mean trait values
  # Not weighted in our case, so just the mean
  # KEEP IN MIND THESE ARE LOG10-TRANSFORMED VALUES
  community.means <- output.combined[[1]]$CWM
  write.csv(community.means,
            file = paste0(fd.dir, "test_community_trait_means_", id, ".csv"),
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

# Parallel apply SampleFD to each sample gapfilled trait dataset
parLapply(cluster,
          sample.list,
          function(x) {
            cat("sample 1\n")
            SampleFD(trait.matrix = traits.gapfilled[[x]],
                     pres.abs.matrix = pres.abs.matrix,
                     id = x
                     )
          }
          )

# Gracefully end cluster
stopCluster(cluster)


# SampleFD for genus-mean filled data
# -----------------------------------
SampleFD(trait.matrix = traits.mean,
         pres.abs.matrix = pres.abs.matrix,
         id = "genus_mean"
         )

# Summarize results of the 100 stochastic gapfilled datasets
# ----------------------------------------------------------
  StochasticMeans(stat = "",
                  samples = 100,
                  output.dir = "output/test/observed_FD/",
                  header = "test_FD_observed_"
                  )

  StochasticMeans(stat = "",
                  samples = 100,
                  output.dir = "output/test/observed_FD/",
                  header = "test_community_trait_means_"
                  )


cat("Done.\n")

