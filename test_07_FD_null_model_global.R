#############################################################
# Palm FD project: Functional Diversity calculation test code
#############################################################
#
# In which we generate FD null model one, based on 
# Random communities sampled from the same species pool.
# The definiton of species pool for this null model is all palm species
# occurring worldwide.
# Which could be the called the "global" null model.
#
# Input files:
#   dir 'output/test/observed_FD/' with 100 observed FD files
#   output/test/observed_FD/test_community_trait_means_genus_mean.csv
#   output/test/test_palm_tdwg3_pres_abs_gapfilled.csv
# Generated output files:
#   "output/test/tdwg3_info.csv"
#   < output dir specified below, with 100 FD z-score files >
#   < nullmodel processing dir specified below, containing for each iteration
#     1 file with null FD values >

cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(FD)
library(parallel)
library(reshape2)

source(file = "functions/base_functions.R")
source(file = "functions/functional_diversity_functions.R")


# Subset FD input for quick code testing?
subset <- TRUE
# Enable verbose reporting in the FD calculation?
verbose <- TRUE
# Null model processing directory (with trailing slash!)
nm.dir <- "output/test/null_model_processing/"
# Null model output directory (with trailing slash!)
output.dir <- "output/test/null_model_global/"
# Number of cores to use for parallel processing. Default is 90% of available cores.
num.cores <- 
  if (!is.na(detectCores())) {
    floor(detectCores() * 0.9)
  } else {
    1
    warning("unable to detect cores. Parallel processing is NOT used!")
    cat("unable to detect cores. Parallel processing is NOT used!")
  }


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

# Load transformed gapfilled trait matrices
# -----------------------------------------
traits.mean <-
  read.csv(file = "output/test/trait_matrices/palm_trait_matrix_genus_mean.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
# This will be useful:
trait.names <- colnames(traits.mean)
traits.gapfilled <- vector("list", length = 100)
for (i in seq_along(traits.gapfilled)) {
  traits.gapfilled[[i]] <-
    read.csv(file = paste0("output/test/trait_matrices/palm_trait_matrix_filled_",
                           i,
                           ".csv"
                           ),
             row.names = 1,
             check.names = FALSE
             ) %>%
    as.matrix()
}

# Load presence/absence matrix
# ----------------------------
pres.abs.matrix <- 
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()


#######################
# The Global Null Model
#######################
cat("Creating Global null model... (this will take a while)\n")

if (!dir.exists(output.dir)) {
  cat("creating directory:", output.dir, "\n")
  dir.create(output.dir)
}

# Define grouping of tdwg3 units into realms
# ------------------------------------------
global.tdwg3 <- list(global = tdwg3.info[, "tdwg3.code"])

# Function to run the global null model
# -------------------------------------
RunGlobal <- function(trait.matrix, pres.abs.matrix, id) {
#   id: a unique identifier, used for naming the produced files.
#       should match the id used for the observed FD file

  # Run null model for the given dataset
  global.gapfilled <-
    NullModel(trait.matrix = trait.matrix,
              pres.abs.matrix = pres.abs.matrix,
              groups = global.tdwg3,
              process.dir = nm.dir,
              iterations = 1000,
              mc.cores = num.cores,
              subset = subset,
              verbose = verbose,
              random.groups = TRUE,
              fast = TRUE
              )

  # Using the null model output, find the z-scores (SES) for the observed FD
  # first, load observed FD indices
  fd.observed <- read.csv(file = paste0("output/test/observed_FD/",
                                        "test_FD_observed_",
                                        id,
                                        ".csv"
                                        ),
                          header = TRUE,
                          row.names = 1
                          )
  # Second, calculate the z-scores
  fd.global <- NullTransform(fd.observed, global.gapfilled)
  # Finally, save results
  cat("Saving output...\n")
  write.csv(fd.global,
            file = paste0(output.dir,
                          "test_FD_z_scores_global_",
                          id,
                          ".csv"
                          ),
            eol = "\r\n",
            row.names = TRUE
            )

  # When all is completed, delete the process.dir
  if (file.exists(paste0(output.dir,
                          "test_FD_z_scores_global_",
                          id,
                          ".csv"
                          )
                   )
      ) {
    if (identical(DetectOS(), "windows")) {
      remove.dir <- substr(nm.dir, (nchar(nm.dir) - 1))
    } else {
      remove.dir <- nm.dir
    }
    unlink(remove.dir, recursive = TRUE)
  }

  return (0)
}


# Global null model for genus-mean filled data
# --------------------------------------------
# This is time-consuming: run only if the output does not yet exist
cat("For genus-mean filled dataset:\n")
id <- "genus_mean"
  if (!file.exists(paste0(output.dir,
                          "test_FD_z_scores_global_",
                          id,
                          ".csv"
                          )
                   )
      ) {
    RunGlobal(trait.matrix = traits.mean,
              pres.abs.matrix = pres.abs.matrix,
              id = id
              )
  }


cat("Done.\n")

