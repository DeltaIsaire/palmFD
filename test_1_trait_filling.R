##########################################
# Palm FD project: trait_filling test code
##########################################
#
# This script generates a number of test scenarios of gap-filling,
# to help with finding the best way to gap-fill the palm traits matrix.
# Consider this the playground, while the main script is only the final
# implementation.
# WARNING: the runtime of this script is long, because it runs BHPMF
# multiple times. 
#
# TESTS: Gap-filling trait matrix
# -------------------------------
# We begin with the base unimputed trait matrix, which has a single 'oberved'
# trait value for each species, with gaps (NAs).

# Then, we consider the following gap-filling scenario's:

# - Gap-filling using genus means
# - Gap-filling with BHPMF, using different inputs:
# 1. Using our three traits of interest, extracting from the BHPMF output only
#    the estimates for the gaps in the original matrix.
# 2. Using our three traits of interest, extracting from the BHPMF output all
#    estimates, including those for 'observed' trait data.
# 3. Same as #1, but including a fourth dummy trait with all values = 1,
#    which circumvents the issue that BHPMF does not work for species
#    with missing values for all traits.
# 4. Same as #1, but including the following additional traits:
#    MaxStemDia_cm, MaxLeafNumber, Max_Rachis_Length_m, Max_Petiole_length_m,
#    AverageFruitWidth_cm
# 5. Same as #3, but with 10 such dummy traits instead of just one.
# 6. Same as #1, but with an extra hierarchy level between genus and species
#    corresponding to growth form: acaulescent / freestanding / climbing
# 7. Same as 4. but including the following additional binary traits, coded
#    with values 0 / 100: Climbing, Acaulescence, Errect, UnderstoryCanopy,
#    StemSolitary, StemArmed, LeavesArmed, FruitSizeBinary.
#    Actually we're going to do 7a: all these traits, and 7b: only the three
#    original traits + the three categorical growthform traits.
#
# Input files:
#   data/palms_in_tdwg3.csv
#   data/PalmTraits_10.csv
# Generated output files:
#   output/palm_traits.csv
#   output/palm_hierarchy.csv
#   output/traits_filled_genus_mean.csv
#   output/missing_genus_mean.csv
#   output/test_BHPMF_one_mean.txt
#   output/test_BHPMF_one_std.txt
#   output/traits_filled_BHPMF_one.csv
#   output/missing_BHPMF_one.csv
#   output/test_BHPMF_two_mean.txt
#   output/test_BHPMF_two_std.txt
#   output/traits_filled_BHPMF_two.csv
#   output/missing_BHPMF_two.csv
#   output/test_BHPMF_three_mean.txt
#   output/test_BHPMF_three_std.txt
#   output/traits_filled_BHPMF_three.csv
#   output/missing_BHPMF_three.csv
#   output/test_BHPMF_four_mean.txt
#   output/test_BHPMF_four_std.txt
#   output/traits_filled_BHPMF_four.csv
#   output/missing_BHPMF_four.csv
#   output/test_BHPMF_five_mean.txt
#   output/test_BHPMF_five_std.txt
#   output/traits_filled_BHPMF_five.csv
#   output/missing_BHPMF_five.csv

cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(BHPMF)

source(file = "functions/base_functions.R")

# Use verbose BHPMF?
verbose <- FALSE

# ------------------
# Prepare trait data
# ------------------
cat("Preparing trait data...\n")
palm.dist <- read.csv(file = "data/palms_in_tdwg3.csv")
# Long-list Data frame with two columns:
#	1. Area_code_L3 - 3-letter code for each botanical country
#	2. SpecName - Species name
trait.data <- 
  read.csv(file = "data/PalmTraits_10.csv") %>%
  .[order(.$SpecName), ]
trait.data$SpecName %<>% sub(pattern = " ",
                             replacement = "_",
                             x = .
                             )
# Data frame with 31 columns:
#	1. SpecName - binomial species name [sorted alphabetically]
#	2. accGenus - only the genus name
#	3. accSpecies - only the species name
#	4. PalmTribe
#	5. PalmSubfamily
#	6-30 - information on 25 traits
#	31. extra reference notes

# Do the datasets agree on which palm species exist?
species.list <- list(palm.dist  = levels(palm.dist$SpecName),
                     trait.data = trait.data$SpecName 
                     )
species.shared <-
  MultiCheck(species.list, presence=TRUE, value=TRUE, unique = TRUE) %>%
  sort()

# Subset palm traits dataset
palm.traits <- 
  data.frame(species      = trait.data$SpecName,
             genus        = trait.data$accGenus,
             stem.height  = trait.data$MaxStemHeight_m,
             blade.length = trait.data$Max_Blade_Length_m,
             fruit.length = trait.data$AverageFruitLength_cm
             ) %>%
  .[CrossCheck(x = .$species,
               y = species.shared,
               presence = TRUE,
               value = FALSE,
               unique = TRUE
               ),
    ]
# stem.height = maximum stem height in meters.
# blade.length = maximum blade length in meters.
# fruit.length = average fruit length in centimeters.
# Genus is included because we need it to calculate genus mean trait values
# Note: there should be no missing species, so the extra subset is for posterity.
write.csv(palm.traits, 
          file = "output/palm_traits.csv",
          eol = "\r\n",
          row.names = FALSE
          )
# This will be useful:
trait.names <- c("stem.height", "blade.length", "fruit.length")

# Extract hierarchical taxonomic data, which is needed for BHPMF
# # columns should be from low-level to high-level
palm.hierarchy <-   data.frame(species   = trait.data$SpecName,
                               genus     = trait.data$accGenus,
                               tribe     = trait.data$PalmTribe,
                               subfamily = trait.data$PalmSubfamily
                               )
write.csv(palm.hierarchy,
          file = "output/palm_hierarchy.csv",
          eol="\r\n",
          row.names = FALSE)


# Gap-filling with genus means
# ----------------------------
cat("Gap-filling with genus means...\n")
genus.means <- ddply(palm.traits, "genus", numcolwise(mean, na.rm=TRUE))
traits.filled.mean <- GapFill(palm.traits,
                              genus.means,
                              by = "genus",
                              fill = trait.names
                              )
# Beware NAs for genus means where all species of the genus have NA for the
# same trait.
# Solution: delete all cases where genus-mean filling is not possible.
genus.mean.missing <- 
  traits.filled.mean %>%
  .[!complete.cases(.), ]
traits.filled.mean %<>% 
  .[complete.cases(.), ]
write.csv(traits.filled.mean,
          file = "output/traits_filled_genus_mean.csv",
          eol="\r\n",
          row.names=FALSE
          )
write.csv(genus.mean.missing, 
          file = "output/missing_genus_mean.csv",
          eol="\r\n",
          row.names=FALSE
          )

# Prepare shared input for BHPMF runs
# -----------------------------------
cat("Preparing data for BHPMF...\n")
# Source dataframes are sorted by species name, so rows will correspond
# automatically.
trait.matrix <- 
  palm.traits[, trait.names] %>%
  as.matrix()
rownames(trait.matrix) <- trait.data$SpecName
# BHPMF is unable to handle observed trait values of exactly 0.
# These do occur in our dataset, for stem height.
# Solution: set these values to something close to 0
trait.matrix[which(trait.matrix == 0)] <- 0.0001

hierarchy.matrix <- as.matrix(palm.hierarchy)

TestBHPMF <- function(trial) {
# Standard code to run a BHPMF test trial.
# Trial: character vector of length 1 giving the name of the trial,
#        e.g. "one".
  # BHPMF wants a preprocessing directory, where it saves pre-processed files.
  # To avoid errors or erronous output when re-running the code, this directory
  # needs to be emptied.
  cat("Running BHPMF trial", trial, "\n")
  unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
  dir.create("output/BHPMF_preprocessing_test")
  GapFilling(X = test.matrix,
             hierarchy.info = test.hierarchy,
             prediction.level = 4,
             used.num.hierarchy.levels = 3,
             mean.gap.filled.output.path = paste0("output/test_BHPMF_",
                                                  trial,
                                                  "_mean.txt"
                                                  ),
             std.gap.filled.output.path = paste0("output/test_BHPMF_",
                                                 trial,
                                                 "_std.txt"
                                                 ),
             tmp.dir = "output/BHPMF_preprocessing_test", 
             rmse.plot.test.data=FALSE,
             verbose=verbose
             )
}

ParseBHPMF <- function(trial) {
# Standard code for parsing BHPMF output. Gapfills the trait matrix
# only for unobserved (missing) trait values.
# Trial: character vector of length 1 giving the name of the trial,
#        e.g. "one".
  traits.filled.BHPMF.test <-
    read.table(file = paste0("output/test_BHPMF_",
                             trial,
                             "_mean.txt"
                             ),
               header = TRUE,
               sep = "	"
               ) %>%
    data.frame(species = as.character(test.hierarchy[, 1]),
               genus   = as.character(test.hierarchy[, 2]),
               .
               ) %>%
    GapFill(palm.traits,
            .,
            by = "species",
            fill = trait.names
            ) %>%
    { .[complete.cases(.), ] }
  
  write.csv(traits.filled.BHPMF.test,
            file = paste0("output/traits_filled_BHPMF_",
                          trial,
                          ".csv"
                          ),
            eol="\r\n",
            row.names=FALSE
            )
  write.csv(missing.BHPMF.test,
            file = paste0("output/missing_BHPMF_",
                          trial,
                          ".csv"
                          ),
            eol="\r\n",
            row.names=FALSE
            )
}

# 1. Gap-filling with BHPMF: only estimates for missing values
# ------------------------------------------------------------
cat("Gap-filling with BHPMF... (this may take a while)\n")
cat("(1) Using only estimates for missing values\n")
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
to.remove <- 
  trait.matrix %>%
  { which(rowSums(is.na(.)) == ncol(.)) }
test.matrix <- trait.matrix[-to.remove, ]
test.hierarchy <- hierarchy.matrix[-to.remove, ]
missing.BHPMF.test <- palm.traits[to.remove, ]
rm(to.remove)

trial <- "one"
TestBHPMF(trial)
ParseBHPMF(trial)

# 2. Gap-filling with BHPMF: using all estimates
# ----------------------------------------------
cat("(2) Using estimates for all trait values\n")
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
to.remove <- 
  trait.matrix %>%
  { which(rowSums(is.na(.)) == ncol(.)) }
test.matrix <- trait.matrix[-to.remove, ]
test.hierarchy <- hierarchy.matrix[-to.remove, ]
missing.BHPMF.test <- palm.traits[to.remove, ]
rm(to.remove)

trial <- "two"
TestBHPMF(trial)
# Not using ParseBHPMF() because here we take all estimates.
traits.filled.BHPMF.test <-
  read.table(file = paste0("output/test_BHPMF_",
                           trial,
                           "_mean.txt"
                           ),
             header = TRUE,
             sep = "	"
             ) %>%
  data.frame(species = as.character(test.hierarchy[, 1]),
             genus   = as.character(test.hierarchy[, 2]),
             .
             ) %>%
  { .[complete.cases(.), ] }

write.csv(traits.filled.BHPMF.test,
          file = paste0("output/traits_filled_BHPMF_",
                        trial,
                        ".csv"
                        ),
          eol="\r\n",
          row.names=FALSE
          )
write.csv(missing.BHPMF.test,
          file = paste0("output/missing_BHPMF_",
                        trial,
                        ".csv"
                        ),
          eol="\r\n",
          row.names=FALSE
          )

# 3. Gap-filling with BHPMF: including 1 dummy trait
# --------------------------------------------------
cat("(3) including 1 dummy trait\n")
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
test.matrix <-
  trait.matrix %>%
  data.frame(.,
             dummy = rep(1, times = length(.[, 1]))
             ) %>%
  as.matrix()
test.hierarchy <- hierarchy.matrix
missing.BHPMF.test <- "none"

trial <- "three"
TestBHPMF(trial)
ParseBHPMF(trial)

# 4. Gap-filling with BHPMF: including more traits
# ------------------------------------------------
cat("(4) including additional palm traits\n")
test.matrix <- 
  cbind(trait.matrix,
        trait.data[, c("MaxStemDia_cm",
                       "MaxLeafNumber",
                       "Max_Rachis_Length_m",
                       "Max_Petiole_length_m",
                       "AverageFruitWidth_cm"
                       )
                   ]
        ) %>%
  as.matrix()

# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
to.remove <- 
  test.matrix %>%
  { which(rowSums(is.na(.)) == ncol(.)) }
test.matrix %<>% .[-to.remove, ]
test.hierarchy <- hierarchy.matrix[-to.remove, ]
missing.BHPMF.test <- palm.traits[to.remove, ]
rm(to.remove)

trial <- "four"
TestBHPMF(trial)
ParseBHPMF(trial)

# 5. Gap-filling with BHPMF: including 10 dummy traits
# ----------------------------------------------------
cat("(5) including 10 dummy traits\n")
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
test.matrix <-
  trait.matrix %>%
  data.frame(.,
             dummy.0 = rep(1, times = length(.[, 1])),
             dummy.1 = rep(1, times = length(.[, 1])),
             dummy.2 = rep(1, times = length(.[, 1])),
             dummy.3 = rep(1, times = length(.[, 1])),
             dummy.4 = rep(1, times = length(.[, 1])),
             dummy.5 = rep(1, times = length(.[, 1])),
             dummy.6 = rep(1, times = length(.[, 1])),
             dummy.7 = rep(1, times = length(.[, 1])),
             dummy.8 = rep(1, times = length(.[, 1])),
             dummy.9 = rep(1, times = length(.[, 1]))
             ) %>%
  as.matrix()

test.hierarchy <- hierarchy.matrix
missing.BHPMF.test <- "none"

trial <- "five"
TestBHPMF(trial)
ParseBHPMF(trial)

# 6. Gap-filling with BHPMF: Including growthform in hierarchy
# ------------------------------------------------------------
cat("(6) including growthform as hierarchy level\n")
# Construct categorical variable of growthform and add it to the hierarchy
growthform <- numeric(length = length(trait.data$SpecName))
growthform[which(trait.data$Climbing == 1)] <- "climbing"
growthform[which(trait.data$Acaulescence == 1)] <- "acaulescent"
growthform[which(trait.data$Errect == 1)] <- "freestanding"
growthform[which(growthform == 0)] <- NA
# Decisions decisions: can we code the '2's as "ambiguous"? Probably not.
# For now they are NA.
# Not done yet: BHPMF demands that lower hierarchy levels are unique to their
# parent. So 'growthform' must be recoded to 'growthform.genus' for each genus.
# Easiest to re-construct the matrix and recode at the same time:
hierarchy.matrix <- data.frame(species = palm.hierarchy$species,
                               growthform = growthform,
                               palm.hierarchy[, -1]
                               )
hierarchy.matrix$growthform %<>% as.character()  # friggin' factors :rolleyes:
hierarchy.matrix <- 
  ddply(hierarchy.matrix,
        "genus",
        function(subset) {
          genus <- as.character(subset$genus[1])
          for (i in 1:length(subset$growthform)) {
            if (!is.na(subset$growthform[i])) {
              subset$growthform[i] %<>% paste0(., ".", genus)
            }
          }
          subset
        }
        ) %>%
  as.matrix()

# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
# Additionally, remove remaining species with NA for growthform.
to.remove <- 
  trait.matrix %>% 
  { which(rowSums(is.na(.)) == ncol(.)) }
to.remove %<>% 
  c(., which(is.na(hierarchy.matrix[, 2]))) %>%
  unique()
test.matrix <- trait.matrix[-to.remove, ]
test.hierarchy <- hierarchy.matrix[-to.remove, ]
missing.BHPMF.test <- palm.traits[to.remove, ]
rm(to.remove)

trial <- "six"
# Not using TestBHPMF() because we use a different hierarchy
cat("Running BHPMF trial", trial, "\n")
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X = test.matrix,
           hierarchy.info = test.hierarchy,
           prediction.level = 5,
           used.num.hierarchy.levels = 4,
           mean.gap.filled.output.path = paste0("output/test_BHPMF_",
                                                trial,
                                                "_mean.txt"
                                                ),
           std.gap.filled.output.path = paste0("output/test_BHPMF_",
                                               trial,
                                               "_std.txt"
                                               ),
           tmp.dir = "output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE,
           verbose=verbose
           )
ParseBHPMF(trial)

# 7. Gap-filling with BHPMF: including even more traits
# ------------------------------------------------
cat("(7) including binary palm traits\n")
# Extract and clean trait data. Binary traits should get values 1 (for 0)
# and 100 (for 1).
extra.traits <- trait.data[, c("MaxStemDia_cm",
                               "MaxLeafNumber",
                               "Max_Rachis_Length_m",
                               "Max_Petiole_length_m",
                               "AverageFruitWidth_cm",
                               "Climbing",
                               "Acaulescence",
                               "Errect",
                               "UnderstoreyCanopy",
                               "StemSolitary",
                               "StemArmed",
                               "LeavesArmed",
                               "FruitSizeBinary"
                               )
                           ]
extra.traits$UnderstoreyCanopy %<>% as.character()
extra.traits$FruitSizeBinary %<>% as.character()

extra.traits$UnderstoreyCanopy %<>%
  { ifelse(. == "canopy", 1, ifelse(. == "understorey", 0, NA)) }
extra.traits$FruitSizeBinary %<>%
  { ifelse(. == "large", 1, ifelse(. == "small", 0, NA)) }

for (i in 6:13) {
  extra.traits[, i] %<>%
    { ifelse(. == 1, 100, ifelse(. == 0, 1, NA)) }
}

test.matrix.all <- cbind(trait.matrix, as.matrix(extra.traits))

test.matrix.growthform <- test.matrix.all[, c(1:3, 9:11)]

cat("(7a) three continuous traits + three growthform binary traits\n")
test.matrix <- test.matrix.growthform
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
to.remove <- 
  test.matrix %>%
  { which(rowSums(is.na(.)) == ncol(.)) }
if (!length(to.remove) == 0) {
  test.matrix %<>% .[-to.remove, ]
  test.hierarchy <- hierarchy.matrix[-to.remove, ]
  missing.BHPMF.test <- palm.traits[to.remove, ]
} else {
  test.hierarchy <- as.matrix(palm.hierarchy)
  missing.BHPMF.test <- 1
}
rm(to.remove)

trial <- "seven.a"
TestBHPMF(trial)
ParseBHPMF(trial)

cat("(7b) all usable traits\n")
test.matrix <- test.matrix.all
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
to.remove <- 
  test.matrix %>%
  { which(rowSums(is.na(.)) == ncol(.)) }
if (!length(to.remove) == 0) {
  test.matrix %<>% .[-to.remove, ]
  test.hierarchy <- hierarchy.matrix[-to.remove, ]
  missing.BHPMF.test <- palm.traits[to.remove, ]
} else {
  test.hierarchy <- as.matrix(palm.hierarchy)
  missing.BHPMF.test <- 1
}
rm(to.remove)

trial <- "seven.a"
TestBHPMF(trial)
ParseBHPMF(trial)


cat("Done.\n")
