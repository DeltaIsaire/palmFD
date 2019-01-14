#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we generate FD null model one, based on 
# Random communities sampled from the same species pool.
# The definiton of species pool for this null model is is all species
# occurring in the same realm (New World, Old World East or Old World West).
# Which could be the called the "regional" null model.
#
# Input files:
#   data/TDWG_Environment_AllData_2014Dec.csv
#   data/palms_in_tdwg3.csv
#   output/test/test_palm_tdwg3_pres_abs_matrix.csv
#   output/test/test_palm_trait_matrix_transformed.csv
# Generated output files:
#   output/test/tdwg3_info.csv
#   output/test/test_fd.indices_nullmodel_one.csv
#   output/test/test_fd.indices_nullmodel_one_single_traits.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)

source(file = "functions/base_functions.R")

# Enable verbose reporting in the null model procedure?
verbose <- TRUE
# Number of cores to use for parallel processing.
# Default is 1 less than the number of available cores (so your computer
# doesn't lock up completely while running the code)
num.cores <- getOption("mc.cores", 2L) - 1


##################
# Data preparation
##################
cat("Preparing data...\n")

# Information on tdwg3 units
# --------------------------
env.data <- read.csv(file = "data/TDWG_Environment_AllData_2014Dec.csv",
                     header = TRUE)
tdwg3.info <- data.frame(tdwg3.code    = env.data$LEVEL_3_CO,
                         tdwg3.name    = env.data$LEVEL_NAME,
                         realm         = env.data$THREEREALM,
                         flora.region  = env.data$REALM_LONG,
                         lat           = env.data$LAT,
                         lon           = env.data$LONG,
                         is.island     = as.factor(env.data$ISISLAND),
                         has.palms     = as.factor(env.data$PalmPresAbs),
                         palm.richness = env.data$PALMSR
                         )
tdwg3.info %<>% .[order(.$tdwg3.code), ]

# Realm has one empty value, belonging to Antarctica, which isn't really in
# one of the three realms. So set the value to NA
tdwg3.info[which(tdwg3.info$realm == ""), "realm"] <- NA
tdwg3.info %<>% droplevels()
# Overwrite $palm.richness and $has.palms with our palm distribution data
palm.dist <- read.csv(file="data/palms_in_tdwg3.csv")
richness <- ddply(palm.dist, "Area_code_L3", nrow)
names(richness) <- c("tdwg3.code", "palm.richness")
missing.codes <- CrossCheck(x = tdwg3.info$tdwg3.code,
                            y = richness$tdwg3.code,
                            presence = FALSE,
                            value = TRUE)
richness %<>% rbind(.,
                    data.frame(tdwg3.code = missing.codes,
                               palm.richness = rep(0, length(missing.codes))
                               )
                    )
richness %<>% .[order(as.character(.$tdwg3.code)), ]
# We have no environmental data for tdwg3 unit 'VNA', so remove it
richness %<>% .[-which(.$tdwg3.code == "VNA"), ]
tdwg3.info$palm.richness <- richness$palm.richness
tdwg3.info$has.palms <-
  ifelse(tdwg3.info$palm.richness > 0, 1, 0) %>%
  as.factor()
write.csv(tdwg3.info,
          file = "output/test/tdwg3_info.csv",
          eol = "\r\n",
          row.names = FALSE
          )
tdwg3.info
# data frame with 368 observations of 9 variables:
# $tdwg3.code   - 3-letter code for each botanical country (tdwg3 unit).
# $tdwg3.name   - Name of the botanical country.
# $realm        - Factor with 3 levels "NewWorld", "OWEast", "OWWest".
# $flora.region - Floristic region to which the botanical country belongs
# $lat          - Latitude.
# $lon          - #Longitude
# is.island     - binary factor indicating whether botanical country is an island
# has.palms     - binary factor indicating whether botanical country has palms
# palm.richness - Palm species richness in the botanical country, including
#                 known absence (richness 0).

# Extract species list for each realm
# -----------------------------------
# First, extract list of tdwg3 units for each realm that have palms.
# Do not use tdwg3.info$has.palms, but use the presence/absence matrix,
# because trait-filling may have excluded some species.
pres.abs.matrix <- 
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_matrix.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
palm.richness <- data.frame(tdwg3.code    = rownames(pres.abs.matrix),
                            palm.richness = rowSums(pres.abs.matrix)
                            )
realm.tdwg3 <- list(new.world = tdwg3.info[which(tdwg3.info$realm == "NewWorld"),
                                           "tdwg3.code"],
                    old.world.west = tdwg3.info[which(tdwg3.info$realm == "OWWest"),
                                                "tdwg3.code"],
                    old.world.east = tdwg3.info[which(tdwg3.info$realm == "OWEast"),
                                                "tdwg3.code"]
                    )
realm.tdwg3 %<>% {
  llply(.,
        function(realm) {
          CrossCheck(x = rownames(pres.abs.matrix),
                     y = realm,
                     presence = TRUE,
                     value = TRUE
                     )
        }
        )
}
# Then, for each realm, extract the present species
realm.species <- llply(realm.tdwg3,
                       function(realm) {
                         indices <- CrossCheck(x = rownames(pres.abs.matrix),
                                               y = realm,
                                               presence = TRUE,
                                               value = FALSE)
                         pres.abs.subset <- pres.abs.matrix[indices, ]
                         species <-
                           pres.abs.subset %>%
                           { colnames(.)[which(colSums(.) > 0)] }
                         species
                       }
                       )

# Load trait matrix
# -----------------
trait.matrix <- 
  read.csv(file = "output/test/test_palm_trait_matrix_transformed.csv",
           row.names = 1) %>%
  as.matrix()
# This will be useful:
trait.names <- colnames(trait.matrix)


#############################
# Running the Null Model test
#############################
#
# This is a complicated multi-step procedure.
# The code is wrapped inside a function.

NullModelRegional <- function(verbose) {
#
# Args:
#   verbose: output extended info about progress.
#
# Returns:
#

  # ----------------------------------------
  # Randomly sampling null model communities
  # ----------------------------------------
  if (verbose) {
    cat("Randomly sampling null model communities from realm species pool...\n")
  }
  # For each tdwg3 unit in each realm, get the palm richness and sample that many
  # species from the realm species pool, without replacement,
  # and combine results into a list.

  # The FD requirement species > traits applies to functionally unique species.
  # That means it does NOT suffice to have richness > 3 (our number of traits).
  # Rather, we need > 3 species with unique trait combinations. 
  # This happens to be the case for our real dataset, but is not automatically
  # guaranteed for the null model samples.
  # So we have to enforce this explicitly, by resampling until this condition
  # is satisfied.
  do.sample <- TRUE
  sampling.iterations <- 0
  while (do.sample) {
    null.model.species <- list(NULL)
    i <- 1
    for (realm in seq_along(realm.tdwg3)) {
      for (country in seq_along(realm.tdwg3[[realm]])) {
        richness <-
          palm.richness %>% {
            .[which(.$tdwg3.code == realm.tdwg3[[realm]][country]),
              "palm.richness"
              ]
          }
        indices <-
          runif(n = richness * 3,
                min = 1,
                max = length(realm.species[[realm]])
                ) %>%
          round() %>%
          unique() %>%
          .[1:richness]
        null.model.species[[i]] <- realm.species[[realm]][indices]
        names(null.model.species)[i] <- realm.tdwg3[[realm]][country]
        i <- i + 1
      }
    }
    null.model.species %<>% .[order(names(.))]
    # Check if uniqueness condition is satisfied:
    # Step 1: assemble trait values for species in each community
    sample.traits <-
      llply(null.model.species,
            function(community) {
              indices <- CrossCheck(x = rownames(trait.matrix),
                                    y = community,
                                    presence = TRUE,
                                    value = FALSE
                                    )
              trait.matrix[indices, ]
            }
            )
    # Step 2: Identify communities with < 4 unique trait combinations
    sample.counts <-
      llply(sample.traits,
            function(community) {
              count(community) %>%
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
  }
  if (verbose) {
    cat("Null model sampling required", sampling.iterations, "sampling attempts\n")
  }
  # Now the fun part: transform null.model.species to a presence/absence matrix
  species <- unlist(null.model.species)
  community <- character(length = length(species))
  index <- 0
  for (i in seq_along(null.model.species)) {
    indices <- seq_along(null.model.species[[i]]) + index
    community[indices] <- rep(names(null.model.species)[i],
                              length(null.model.species[[i]])
                              )
    index <- max(indices)
  }
  nm.one.species <- data.frame(community = community,
                               species = species
                               )
  nm.one.pres.abs <-
    table(nm.one.species) %>%
    as.data.frame.matrix() %>%
    as.matrix()

  # --------------------------------
  # Calculating Functional Diversity
  # --------------------------------
  if (verbose) {
    cat("Calculating Functional Diversity indices... (this may take a while)\n")
  }
  # Subset trait matrix to species in the null model
  indices <- CrossCheck(x = rownames(trait.matrix),
                      y = colnames(nm.one.pres.abs),
                      presence = TRUE,
                      value = FALSE
                      )
  trait.matrix %<>% .[indices, ]

  # Code to subset data to expedite testing. This stuff is computationally
  # intensive.
  if (FALSE) {
    nm.one.pres.abs <- nm.one.pres.abs[1:10, ]
    orphaned.species <-
      which(colSums(nm.one.pres.abs) == 0) %>%
      colnames(nm.one.pres.abs)[.]
    nm.one.pres.abs %<>% .[, -which(colSums(.) == 0)]
    indices <- CrossCheck(x = rownames(trait.matrix),
                          y = orphaned.species,
                          presence = TRUE,
                          value = FALSE
                          )
    trait.matrix <- trait.matrix[-indices, ]
  }

# FD using all traits
# -------------------
cat("(1) FD with all traits...\n")
output.all <- suppressWarnings(dbFD(x = trait.matrix,
                                    a = nm.one.pres.abs,
                                    w.abun = FALSE,
                                    stand.x = TRUE,
                                    corr = "cailliez",
                                    calc.FRic = TRUE,
                                    m = "max",
                                    stand.FRic = FALSE,
                                    calc.CWM = FALSE,
                                    calc.FDiv = FALSE,
                                    messages = TRUE
                                    )
                               )
# FD for single traits
# --------------------
cat("(2) FD for single traits:\n")

Temp <- function(trait) {
  cat(trait, "...\n")
  subset.matrix <-
     trait.matrix[, trait] %>%
     as.matrix()
  colnames(subset.matrix) <- trait
  fd.trait <- dbFD(x = subset.matrix,
                   a = nm.one.pres.abs,
                   w.abun = FALSE,
                   stand.x = TRUE,
                   corr = "cailliez",  # Is this the best option?
                   calc.FRic = TRUE,
                   m = "max",
                   stand.FRic = FALSE,  # Would this be useful to do?
                   calc.CWM = FALSE,
                   calc.FDiv = FALSE,
                   messages = TRUE
                   )
  fd.trait
}

output.stem.height  <- suppressWarnings(Temp(trait.names[1]))
output.blade.length <- suppressWarnings(Temp(trait.names[2]))
output.fruit.length <- suppressWarnings(Temp(trait.names[3]))


# ----------------------
# Parse and save results
# ----------------------
cat("Parsing and saving output...\n")

# FD indices for all traits
fd.indices <- data.frame(TDWG3 = names(output.all$nbsp),
                         FRic  = output.all$FRic,
                         FDis  = output.all$FDis
                         )
write.csv(fd.indices,
          file = "output/test/test_fd.indices_nullmodel_one.csv",
          eol = "\r\n",
          row.names = FALSE
          )

# FD indices for single traits
single.fd <- data.frame(TDWG3       = names(output.all$nbsp),
                        height.FRic = output.stem.height$FRic,
                        height.FDis = output.stem.height$FDis,
                        blade.FRic  = output.blade.length$FRic,
                        blade.FDis  = output.blade.length$FDis,
                        fruit.FRic  = output.fruit.length$FRic,
                        fruit.FDis  = output.fruit.length$FDis
                        )
write.csv(single.fd,
          file = "output/test/test_fd.indices_nullmodel_one_single_traits.csv",
          eol = "\r\n",
          row.names = FALSE
          )

cat("Done.\n")

