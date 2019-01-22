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
#   output/test/test_palm_tdwg3_pres_abs_gapfilled.csv
#   output/test/test_palm_tdwg3_pres_abs_unfilled.csv
#   output/test/test_palm_traits_transformed_gapfilled.csv
#   output/test/test_palm_traits_transformed_unfilled.csv
#   output/test/test_fd_indices_gapfilled.csv
#   output/test/test_fd_indices_unfilled.csv
# Generated output files:
#   output/test/tdwg3_info.csv
#   Null model processing directory: output/test/nullmodel_regional/
#   output/test/test_fd_z_scores_regional_gapfilled.csv
#   output/test/test_fd_z_scores_regional_unfilled.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)
library(reshape2)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")


# Enable verbose reporting in the null model procedure?
verbose <- TRUE
# Number of cores to use for parallel processing (via parallel::mclapply).
# Default is 1 less than the number of available cores (so your computer
# doesn't lock up completely while running the code), and should work on most
# machines.
num.cores <- if (is.na(detectCores())) { 1 } else { max(1, detectCores() - 1) }
# Subset FD input for quick code testing?
subset <- FALSE


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
# any of the three realms. So set the value to NA
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

# Load datasets
# -------------
pres.abs.gapfilled <- 
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
pres.abs.unfilled <- 
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_unfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
traits.gapfilled <- 
  read.csv(file = "output/test/test_palm_traits_transformed_gapfilled.csv",
           row.names = 1
           ) %>%
  as.matrix()
traits.unfilled <- 
  read.csv(file = "output/test/test_palm_traits_transformed_unfilled.csv",
           row.names = 1
           ) %>%
  as.matrix()
# This will be useful:
trait.names <- colnames(traits.gapfilled)


#########################
# The Regional Null Model
#########################
cat("Creating Regional null model... (this will take a while)\n")

# Define grouping of tdwg3 units into realms
# ------------------------------------------
realm.tdwg3 <- list(new.world = tdwg3.info[tdwg3.info$realm == "NewWorld",
                                           "tdwg3.code"],
                    old.world.west = tdwg3.info[tdwg3.info$realm == "OWWest",
                                                "tdwg3.code"],
                    old.world.east = tdwg3.info[tdwg3.info$realm == "OWEast",
                                                "tdwg3.code"]
                    )

# Run the null model simulations
# ------------------------------
regional.gapfilled <-
  NullModel(trait.matrix = traits.gapfilled,
            pres.abs.matrix = pres.abs.gapfilled,
            groups = realm.tdwg3,
            process.dir = "output/test/nullmodel_regional/gapfilled/",
            iterations = 100,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = TRUE
            )

regional.unfilled <-
  NullModel(trait.matrix = traits.unfilled,
            pres.abs.matrix = pres.abs.unfilled,
            groups = realm.tdwg3,
            process.dir = "output/test/nullmodel_regional/unfilled/",
            iterations = 100,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = TRUE
            )


############################################################
# Use null model output to convert raw FD values to z-scores
############################################################
cat("Applying regional null model to raw FD indices...\n")

# Load raw FD indices
# -------------------
fd.gapfilled <- read.csv(file = "output/test/test_fd_indices_gapfilled.csv",
                      header = TRUE,
                      row.names = 1
                      )
fd.unfilled <- read.csv(file = "output/test/test_fd_indices_unfilled.csv",
                        header = TRUE,
                        row.names = 1
                        )

# Convert FD to z-cores and save output
# -------------------------------------
fd.regional.gapfilled <- NullTransform(fd.gapfilled, regional.gapfilled)
  write.csv(fd.regional.gapfilled,
            file = "output/test/test_fd_z_scores_regional_gapfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )
fd.regional.unfilled <- NullTransform(fd.unfilled, regional.unfilled)
  write.csv(fd.regional.unfilled,
            file = "output/test/test_fd_z_scores_regional_unfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )


############################################
# Determine optimal iterations for final run
############################################
cat("Assessing z-score convergence...\n")
# How many iterations is enough to get convergence on a z-score?
#
# We can test this quickly by running the null model n times for a few hand-picked
# communities.
# Which communities? The ones with highest/lowest raw FRic and FDis.
# That's MLY, WAU, SCZ and TCI, respectively.

# Subset data to these tdwg3 units
# --------------------------------
pres.abs.subset <- 
  pres.abs.gapfilled[rownames(pres.abs.gapfilled) %in% c("MLY",
                                                         "WAU",
                                                         "SCZ",
                                                         "TCI"
                                                          ),
                     ]
fd.test <- fd.gapfilled[rownames(fd.gapfilled) %in% rownames(pres.abs.subset),
                        c(1, 5)
                        ]
orphaned.species <-
  which(colSums(pres.abs.subset) == 0) %>%
  colnames(pres.abs.subset)[.]
pres.abs.subset %<>% .[, -which(colSums(.) == 0)]
indices <- CrossCheck(x = rownames(traits.gapfilled),
                      y = orphaned.species,
                      presence = TRUE,
                      value = FALSE
                      )
traits.subset <- traits.gapfilled[-indices, ]
realm.test <- realm.tdwg3[c(1, 3)]

# Run regional null model for n iterations and compute z-values
# -------------------------------------------------------------
# We should test for iterations in (100 * 1:5)
# To save processing time we can cheat: instead of running ten times for ten
# different iterations, we run it once for 500 iterations, and then
# calculate z-scores based on partial output.
# We do still need say 10 z-scores for each number of iterations, which means
# running 500 iterations ten times.
# To further save on processing time, exclude single-trait FD and focus on
# all-traits FD.
#
# Only run code if the output doesn't exist
if (!file.exists("output/test/nullmodel_regional_z_convergence_FRic.csv")) {
  test.results <- vector("list", length = 10)

  for (run in 1:10) {
    # init partial results list
    partial.results <- unlist(list(rep(list(fd.test), 5)), recursive = FALSE)
    names(partial.results) <- paste0("z.score.", 100 * 1:5)
    # clean processing dir
    if (dir.exists("output/test/nullmodel_test/")) {
      unlink("output/test/nullmodel_test/", recursive = TRUE)
    }
    # run nullmodel
    regional.test <-
      NullModel(trait.matrix = traits.subset,
                pres.abs.matrix = pres.abs.subset,
                groups = realm.test,
                process.dir = "output/test/nullmodel_test/",
                iterations = 500,
                mc.cores = num.cores,
                subset = subset,
                verbose = verbose,
                random.groups = TRUE,
                single.traits = FALSE
                )
    # compute z-scores
    for (i in 1:5) {
      # subset nullmodel output
      regional.test.subset <-
        llply(regional.test, function(x) { x[[1]] [1:(100 * i), ] } )
      names(regional.test.subset) <- c("FRic.all.traits", "FDis.all.traits")
      partial.results[[i]] <- NullTransform(fd.test, regional.test.subset)
    }
    test.results[[run]] <- partial.results
  }
  
  # Compute sd of z-values for each number of iterations
  # ----------------------------------------------------
  # Init result list
  mat <- matrix(ncol = 5,
                nrow = nrow(fd.test),
                dimnames = list(row = rownames(fd.test),
                                col = paste0("iterations.", 100 * 1:5)
                                )
                )
  sd.results <- list(FRic = mat, FDis = mat)
  mean.results <- list(FRic = mat, FDis = mat)
  # apply magic
  values <- numeric(length(test.results))
  for (run in seq_along(test.results)) {
    for (index in 1:2) {
      for (iterations in 1:5) {
        for (area in seq_len(nrow(fd.test))) {
          values[i] <- test.results[[run]] [[iterations]] [area, index]
          sd.results[[index]] [area, iterations] <- sd(values)
          mean.results[[index]] [area, iterations] <- mean(values)
        }
      }
    }
  }
  sd.results
  mean.results

  # Save results
  # ------------
  write.csv(sd.results[[1]],
            file = "output/test/nullmodel_regional_z_convergence_FRic_sd.csv",
            eol = "\r\n",
            row.names = TRUE
            )
  write.csv(sd.results[[2]],
            file = "output/test/nullmodel_regional_z_convergence_FDis_sd.csv",
            eol = "\r\n",
            row.names = TRUE
            )
  write.csv(mean.results[[1]],
            file = "output/test/nullmodel_regional_z_convergence_FRic_mean.csv",
            eol = "\r\n",
            row.names = TRUE
            )
  write.csv(mean.results[[2]],
            file = "output/test/nullmodel_regional_z_convergence_FDis_mean.csv",
            eol = "\r\n",
            row.names = TRUE
            )
}
# Interpretation of sd.results:
# FRic - WAU is stable after 200 iterations.
#        TCI is stable after 300 iterations.
#        SCZ is stable after 100 iterations.
#        MLY may or may not be stable after 500 iterations, but the sd is
#        very low - an order of magnitude lower than the next lowest (WAU) - so
#        it is probably fine.
# FDis - WAU is stable after 100 iterations.
#        TCI is stable after 100 iterations.
#        SCZ is stable after 100-200 iterations.
#        MLY is stable after 400 iterations.
# Overall, something like 300 should be quite sufficient.
# In fact the difference between 100 iterations sd or 500 iterations sd is very
# small (less than 10% on average). A few hundred iterations will be absolutely
# fine. 500 or more would be overkill for all practical purposes.
#
# Interpretation of mean.results:
# FRic - WAU is stable after 100-200 iterations
#        TCI is stable after 200 iterations
#        SCZ is stable after 100-300 iterations
#        MLY is stable after 300 iterations
# FDis - WAU is stable after 100 iterations
#        TCI is stable after 100 iterations
#        SCZ is stable after 100 iterations
#        MLY is stable after 400 iterations
# The means tell nearly the same story as the SDs.
#
# Final conclusion: You need at most 400 iterations. 
# 300 would also be fine in the majority of cases.

cat("Done.\n")

