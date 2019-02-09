###################################################
# Palm FD project: Functional Diversity calculation
###################################################
#
# In which we generate FD null model one, based on 
# Random communities sampled from the same species pool.
# The definiton of species pool for this null model is is all species
# occurring in the same realm (New World, Old World East or Old World West).
# Which could be the called the "regional" null model.
#
# Input files:
#   data/TDWG_Environment_AllData_2019Jan.csv
#   data/palms_in_tdwg3.csv
#   output/FD_palm_tdwg3_pres_abs_gapfilled.csv
#   output/FD_palm_tdwg3_pres_abs_unfilled.csv
#   output/FD_palm_traits_transformed_gapfilled.csv
#   output/FD_palm_traits_transformed_unfilled.csv
#   output/FD_fd_indices_gapfilled.csv
#   output/FD_fd_indices_unfilled.csv
# Generated output files:
#   output/tdwg3_info.csv
#   Null model processing directory: output/nullmodel_regional/
#   output/FD_z_scores_regional_gapfilled.csv
#   output/FD_z_scores_regional_unfilled.csv


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
env.data <- read.csv(file = "data/TDWG_Environment_AllData_2019Jan.csv",
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
          file = "output/tdwg3_info.csv",
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
  read.csv(file = "output/FD_palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
pres.abs.unfilled <- 
  read.csv(file = "output/FD_palm_tdwg3_pres_abs_unfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
traits.gapfilled <- 
  read.csv(file = "output/FD_palm_traits_transformed_gapfilled.csv",
           row.names = 1
           ) %>%
  as.matrix()
traits.unfilled <- 
  read.csv(file = "output/FD_palm_traits_transformed_unfilled.csv",
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
            process.dir = "output/nullmodel_regional/gapfilled/",
            iterations = 1000,
            mc.cores = num.cores,
            subset = subset,
            verbose = verbose,
            random.groups = TRUE
            )

regional.unfilled <-
  NullModel(trait.matrix = traits.unfilled,
            pres.abs.matrix = pres.abs.unfilled,
            groups = realm.tdwg3,
            process.dir = "output/nullmodel_regional/unfilled/",
            iterations = 1000,
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
fd.gapfilled <- read.csv(file = "output/FD_indices_gapfilled.csv",
                      header = TRUE,
                      row.names = 1
                      )
fd.unfilled <- read.csv(file = "output/FD_indices_unfilled.csv",
                        header = TRUE,
                        row.names = 1
                        )

# Convert FD to z-cores and save output
# -------------------------------------
fd.regional.gapfilled <- NullTransform(fd.gapfilled, regional.gapfilled)
  write.csv(fd.regional.gapfilled,
            file = "output/FD_z_scores_regional_gapfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )
fd.regional.unfilled <- NullTransform(fd.unfilled, regional.unfilled)
  write.csv(fd.regional.unfilled,
            file = "output/FD_z_scores_regional_unfilled.csv",
            eol = "\r\n",
            row.names = TRUE
            )

cat("Done.\n")

