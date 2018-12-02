###################################################
# Palm FD Project: Assessing accuracy of gapfilling
###################################################
#
# In which the accuracy of trait gap-filling is assessed by predicting 
# known values from a subsetted palm trait matrix. We can create a test
# dataset by subsetting the palm traits matrix to only species with complete
# trait data for all our three traits of interest. Next we randomly delete
# some of the trait values, and let the gap-filling methods estimate them.
# Then we can compare the estimates with the 'real' values.
#
# Input files:
#   output/palm_traits.csv
# Generated output files:


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(BHPMF)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


# -----------------------------------
# Data preparation: Test trait matrix
# -----------------------------------
cat("Preparing data...\n")
# Load palm traits data and subset to complete cases:
palm.traits <- read.csv(file = "output/palm_traits.csv")
complete.traits <- palm.traits[complete.cases(palm.traits), ]
# calculate combined sparsity:
sparsity <- 1 - (dim(complete.traits)[1] / dim(palm.traits)[1])

# Create sparse test matrix
# -------------------------
# First, remove genera with only 1 member. A gap in the trait data of such
# genera means estimating a genus mean is not possible.
genera <-
  ddply(complete.traits, "genus", nrow) %>%
  .[!.[, 2] <= 1, ] %>%
  .[, 1] %>%
  as.character()
complete.traits <-
  CrossCheck(complete.traits$genus,
             genera,
             presence = TRUE,
             value = FALSE,
             unique = FALSE
             ) %>%
  complete.traits[., ]

# BHPMF does not like species with NA for all traits.
# So for each species, randomly select 1 trait value that is always present.
to.keep <-
  runif(n = dim(complete.traits)[1],
        min = 3,
        max = 5
        ) %>%
  round()
# The inverse selection is the set of values to create gaps in:
to.sparse <-
  as.list(to.keep) %>%
  llply(., function(x) {
             c(3, 4, 5) %>%
             .[!. %in% x]
           }
        ) %>%
  simplify2array() %>%
  t()
# to.sparse is, for each species, the column indices for trait values
# that *can* be set to NA.
# Now, to create gaps in complete.traits:
#   1. extract the values indicated by to.sparse into a vector
values <- numeric(length = length(to.sparse))
for (i in 1:dim(to.sparse)[1]) {
  val.ind <- c(2 * i - 1, 2 * i)
  values[val.ind] <- as.numeric(complete.traits[i, to.sparse[i, ]])
}
#   2. substitute some of the values in this vector with NA
gaps <-
  runif(n = round(sparsity * length(values)),
        min = 1,
        max = length(values)
        ) %>%
  round()
values[gaps] <- NA
#   3. put the resulting vector back in a copy of complete.traits
sparse.traits <- complete.traits
for (i in 1:dim(to.sparse)[1]) {
  val.ind <- c(2 * i - 1, 2 * i)
  sparse.traits[i, to.sparse[i, ]] <- values[val.ind]
}

# Almost done: verify that a genus mean can be calculated in all cases
means <- ddply(sparse.traits, "genus", numcolwise(mean, na.rm=TRUE))
success <- isTRUE(length(which(!complete.cases(means))) == 0)
# Solution if not successful: for these genera, put the known
# values back in the trait matrix
if (!success) {
  replace <- as.character(means[!complete.cases(means), "genus"])
  indices <- CrossCheck(sparse.traits$genus,
                        replace,
                        presence = TRUE,
                        value = FALSE,
                        unique = FALSE
                        )
  sparse.traits[indices, ] <- complete.traits[indices, ]
}

write.csv(sparse.traits,
          file = "output/test_sparse_traits.csv",
          eol="\r\n",
          row.names=FALSE
          )

