###########################################
# Palm FD project: gap-filling trait matrix
###########################################
#
# In which the three palm trait variables we will use in this study
# (stem height, blade length, fruit length) are extracted from the trait
# dataset, and missing values in these three variables are filled using
# estimates based on (1) Genus-level means and (2) BHPMF.
#
# Input files:
#   data/palms_in_tdwg3.csv
#   data/PalmTraits_10.csv
# Generated output files:
#   output/palm_traits.csv
#   output/palm_hierarchy.csv
#   output/palm_traits_genus_mean.csv
#   output/genus_mean_missing.csv
#   output/BHPMF_preprocessing (directory containing many files)
#   output/BHPMF_mean.txt
#   output/BHPMF_std.txt
#   output/palm_traits_BHPMF.csv
#   output/BHPMF_missing.csv


library(magrittr)
library(plyr)
library(BHPMF)

source(file="functions/base_functions.R")


# ------------------
# Prepare trait data
# ------------------
cat("Preparing trait data...", "\n")
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
species <- list(palm.dist  = levels(palm.dist$SpecName),
                trait.data = trait.data$SpecName 
                )
present.list <- MultiCC(species, presence=TRUE, value=TRUE)
species.agreed <- sort(present.list$palm.dist$trait.data)

# Subset palm traits dataset
palm.traits <- data.frame(species      = trait.data$SpecName,
                          genus        = trait.data$accGenus,
                          stem.height  = trait.data$MaxStemHeight_m,
                          blade.length = trait.data$Max_Blade_Length_m,
                          fruit.length = trait.data$AverageFruitLength_cm
                          )
# stem.height = maximum stem height in meters.
# blade.length = maximum blade length in meters.
# fruit.length = average fruit length in centimeters.
# Genus is included because we need it to calculate genus mean trait values
palm.traits %<>% .[CrossCheck(x = palm.traits$species, 
                              y = species.agreed,
                              presence=TRUE,
                              value=FALSE
                              ),
                   ]
# Note: there should be no missing species, so this subset is for posterity.
write.csv(palm.traits, 
          file = "output/palm_traits.csv",
          eol = "\r\n",
          row.names = FALSE
          )

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




# ---------------------------------------
# Gap-fill trait matrix using genus means
# ---------------------------------------
cat("Gap-filling with genus means...", "\n")
genus.means <- ddply(palm.traits, "genus", numcolwise(mean, na.rm=TRUE))
traits.filled.mean <- GapFill(palm.traits,
                              genus.means,
                              by = "genus",
                              fill = c("stem.height",
                                       "blade.length",
                                       "fruit.length"
                                       )
                              )
# Beware NAs for genus means where all species of the genus have NA for the
# same trait.
# Solution: delete all cases where genus-mean filling is not possible.
genus.mean.missing <- 
  traits.filled.mean %>%
  .[!complete.cases(.), ]
traits.filled.mean %<>% .[complete.cases(.), ]
write.csv(traits.filled.mean,
          file = "output/palm_traits_genus_mean.csv",
          eol="\r\n",
          row.names=FALSE
          )
write.csv(genus.mean.missing, 
          file = "output/genus_mean_missing.csv",
          eol="\r\n",
          row.names=FALSE
          )


# ---------------------------------
# Gap-fill trait matrix using BHPMF
# ---------------------------------
cat("Gap-filling with BHPMF... (this may take a while)", "\n")
# First assemble the prerequisites data matrices.
# Source dataframes are sorted by species name, so rows will correspond
# automatically.
trait.matrix <- 
  palm.traits[, c("stem.height", "blade.length", "fruit.length")] %>%
  as.matrix()
rownames(trait.matrix) <- trait.data$SpecName
# columns of the hierarchy matrix should be from low-level to high-level:
hierarchy.matrix <-
  data.frame(species   = trait.data$SpecName,
             genus     = trait.data$accGenus,
             tribe     = trait.data$PalmTribe,
             subfamily = trait.data$PalmSubfamily
             ) %>%
  as.matrix()
# BHPMF does not want observations with NA for all trait values,
# So we have to remove those from both matrices.
to.remove <- 
  trait.matrix %>% { which(rowSums(is.na(.)) == ncol(.)) }
BHPMF.missing <- palm.traits[to.remove, ]
trait.matrix %<>% .[-to.remove, ]
hierarchy.matrix %<>% .[-to.remove, ]
rm(to.remove)
# BHPMF is unable to handle observed trait values of exactly 0.
# These do occur in our dataset, for stem height.
# Solution: set these values to something close to 0
trait.matrix[which(trait.matrix == 0)] <- 0.0001
# BHPMF wants a preprocessing directory, where it saves pre-processed files.
# To avoid errors or erronous output when re-running the code, this directory
# needs to be emptied.
unlink("output/BHPMF_preprocessing", recursive=TRUE)
dir.create("output/BHPMF_preprocessing")
GapFilling(X = trait.matrix,
           hierarchy.info = hierarchy.matrix,
           prediction.level = 4,
           used.num.hierarchy.levels = 3,
           mean.gap.filled.output.path = "output/BHPMF_mean.txt",
           std.gap.filled.output.path = "output/BHPMF_std.txt",
           tmp.dir = "output/BHPMF_preprocessing", 
           rmse.plot.test.data=FALSE,
           verbose=FALSE
           )
traits.filled.BHPMF <-
  read.table(file = "output/BHPMF_mean.txt",
             header = TRUE,
             sep = "	"
             ) %>%
  data.frame(species = as.character(hierarchy.matrix[, 1]),
             genus   = as.character(hierarchy.matrix[, 2]),
             .
             ) %>%
  GapFill(palm.traits,
          .,
          by = "species",
          fill = c("stem.height", "blade.length", "fruit.length")
          ) %>%
  { .[complete.cases(.), ] }

write.csv(traits.filled.BHPMF,
          file = "output/palm_traits_BHPMF.csv",
          eol="\r\n",
          row.names=FALSE
          )
write.csv(BHPMF.missing,
          file = "output/BHPMF_missing.csv",
          eol="\r\n",
          row.names=FALSE
          )

cat("Done.", "\n")

# TODO:we use the taxonomic data from the trait dataset.
# Will we do nothing with the phylogenetic tree?
