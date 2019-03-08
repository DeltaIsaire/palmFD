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
#   data/PalmTraits_10.csv
#   output/palm_hierarchy.csv
# Generated output files:
#   output/gapfilling_test_sparse_traits.csv
#   output/gapfilling_test_sparse_traits_complete.csv
#   output/gapfilling_test_sparse_filled_genus_mean.csv
#   directory: output/BHPMF_preprocessing_test
#   output/gapfilling_test_sparse_filled_BHPMF.csv
#   output/gapfilling_test_sparse_estimates_stats_table.csv
#   < a bunch of graphs starting with
#     'graphs/gapfilling_test_sparse_scatter_' >
#   output/gapfilling_test_sparse_est.vs.original.rsq.csv


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(BHPMF)
library(reshape2)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")

# Use verbose BHPMF?
verbose <- FALSE

# -----------------------------------
# Data preparation: Test trait matrix
# -----------------------------------
cat("Preparing data: creating artificially sparse trait matrix...\n")
# Load palm traits data and subset to complete cases:
palm.traits <- read.csv(file = "output/palm_traits.csv")
trait.names <- c("stem.height", "blade.length", "fruit.length")
complete.traits <- palm.traits[complete.cases(palm.traits), ]
# calculate combined sparsity:
sparsity <- 1 - (dim(complete.traits)[1] / dim(palm.traits)[1])
# calculate overall proportion of missing values
gap.count <- length(which(is.na(palm.traits[, trait.names])))
gap.proportion <- gap.count / length((as.matrix(palm.traits[, trait.names])))

# Create sparse test matrix
# -------------------------
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
# The amount of NAs to add is gap.proportion * the number of trait values
# in 'complete.traits'
new.gap.count <-
  complete.traits[, trait.names] %>%
  as.matrix() %>%
  length()
gaps <-
  runif(n = 1.3 * round(gap.proportion * new.gap.count),
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

# New sparsity and gaps:
new.sparsity <- 1 - (dim(sparse.traits)[1] / dim(complete.traits)[1])
# This will be 0 by construction.
new.gap.count <- length(which(is.na(sparse.traits[, trait.names])))
new.gap.proportion <- 
  new.gap.count / length((as.matrix(complete.traits[, trait.names])))

write.csv(sparse.traits,
          file = "output/gapfilling_test_sparse_traits.csv",
          eol = "\r\n",
          row.names = FALSE
          )
write.csv(complete.traits,
          file = "output/gapfilling_test_sparse_traits_complete.csv",
          eol = "\r\n",
          row.names = FALSE
          )


# -----------------------------------
# Gap-filling the sparse trait matrix
# -----------------------------------
cat("Gap-filling the sparse test matrix:\n")

# ----------------------------
# Gap-filling with genus means
# ----------------------------
cat("(1) with genus means...\n")
genus.means <- ddply(sparse.traits, "genus", numcolwise(mean, na.rm=TRUE))
traits.filled.mean <- GapFill(sparse.traits,
                              genus.means,
                              by = "genus",
                              fill = trait.names
                              )
write.csv(traits.filled.mean,
          file = "output/gapfilling_test_sparse_filled_genus_mean.csv",
          eol="\r\n",
          row.names=FALSE
          )


# ----------------------
# Gap-filling with BHPMF
# ----------------------
cat("(2) with BHPMF:\n")

# Prepare shared input for BHPMF runs
# -----------------------------------
cat("Preparing data for BHPMF...\n")
trait.data <- 
  read.csv(file = "data/PalmTraits_10.csv") %>%
  .[order(.$SpecName), ]
trait.data$SpecName %<>% sub(pattern = " ",
                             replacement = "_",
                             x = .
                             )

# Source dataframes are sorted by species name, so rows will correspond
# automatically.
trait.matrix <- 
  sparse.traits[, trait.names] %>%
  as.matrix()
rownames(trait.matrix) <- sparse.traits$species
# BHPMF is unable to handle observed trait values of exactly 0.
# These do occur in our dataset, for stem height.
# Solution: set these values to something close to 0
trait.matrix[which(trait.matrix == 0)] <- 0.0001
test.matrix <- trait.matrix

hierarchy.matrix <- 
  read.csv(file = "output/palm_hierarchy.csv") %>%
  .[CrossCheck(.$species,
               sparse.traits$species,
               presence = TRUE,
               value = FALSE
               ),
    ] %>%
  as.matrix()
test.hierarchy <- hierarchy.matrix

# 1. Gap-filling test data with BHPMF
# -----------------------------------
cat("Gap-filling with BHPMF... (this may take a while)\n")

# BHPMF wants a preprocessing directory, where it saves pre-processed files.
# To avoid errors or erronous output when re-running the code, this directory
# needs to be emptied.
unlink("output/BHPMF_preprocessing_test", recursive = TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X = test.matrix,
           hierarchy.info = test.hierarchy,
           prediction.level = 4,
           used.num.hierarchy.levels = 3,
           mean.gap.filled.output.path = paste0("output/gapfilling_test_BHPMF_",
                                                "mean.txt"
                                                ),
           std.gap.filled.output.path = paste0("output/gapfilling_test_BHPMF_",
                                               "std.txt"
                                               ),
           tmp.dir = "output/BHPMF_preprocessing_test", 
           rmse.plot.test.data = FALSE,
           verbose = verbose
           )
# Parsing results
traits.filled.BHPMF.test <-
  read.table(file = paste0("output/gapfilling_test_BHPMF_",
                                                "mean.txt"
                                                ),
             header = TRUE,
             sep = "	"
             ) %>%
  data.frame(species = as.character(test.hierarchy[, 1]),
             genus   = as.character(test.hierarchy[, 2]),
             .
             ) %>%
  GapFill(sparse.traits,
          .,
          by = "species",
          fill = trait.names
          ) %>%
  { .[complete.cases(.), ] }
write.csv(traits.filled.BHPMF.test,
          file = "output/gapfilling_test_sparse_filled_BHPMF.csv",
          eol = "\r\n",
          row.names = FALSE
          )


# -----------------------------
# Parsing gap-filled trait data
# -----------------------------
cat("Parsing gap-filled trait matrices...\n")

# Load filled trait matrices:
#   genus means
filled.mean <-
  read.csv(file = "output/gapfilling_test_sparse_filled_genus_mean.csv")
#   BHPMF
filled.BHPMF <- read.csv(file = "output/gapfilling_test_sparse_filled_BHPMF.csv")
# combine:
all.filled <- list(original = complete.traits,
                   mean     = filled.mean,
                   BHPMF    = filled.BHPMF
                   )
filled.names <- names(all.filled)

# For each trait, get the subset of values that were estimated
all.estimates <- llply(as.list(trait.names),
                       function(trait) { 
                         missing <- which(is.na(sparse.traits[, trait]))
                         species <- sparse.traits[missing, "species"]
                         estimates <- llply(all.filled,
                                            function(df) {
                                              indices <- CrossCheck(x = df$species,
                                                                    y = species,
                                                                    value = FALSE
                                                                    )
                                              df[indices, trait]
                                            }
                                            )
                         data.frame(species = species,
                                    estimates %>%
                                    simplify2array() %>%
                                    as.data.frame()
                         )
                       }
                       )
names(all.estimates) <- paste0(trait.names, ".estimates")

# For each trait, calculate the differences between estimates and original value
all.differences <- 
  llply(all.estimates,
        function(df) {
          cbind(df[, "species", drop = FALSE],
                df[, -c(1, 2)] - df[, "original"]
                )
        }
        )
names(all.differences) <- paste0(trait.names, ".differences")


# ------------------------------
# Analysing gap-filling accuracy
# ------------------------------
cat("Analysing accuracy of gap-filling...\n")

# First: for each trait, a table of mean and st.dev for all cases
all.means <- llply(all.estimates, numcolwise(mean))
all.sds <- llply(all.estimates, numcolwise(sd))
estimates.stats <- llply(as.list(seq_along(all.estimates)),
                         function(i) {
                           data.frame(id = colnames(all.means[[i]]),
                                      mean = as.numeric(all.means[[i]]),
                                      st.dev = as.numeric(all.sds[[i]])
                                      )
                         }
                         )
names(estimates.stats) <- paste0(trait.names, ".estimates")
# We should save these as a table.
est.stats.table <-
  data.frame("method"              = estimates.stats[[1]][, "id"],
             "mean (stem height)"  = estimates.stats[[1]][, "mean"],
             "sd (stem height)"    = estimates.stats[[1]][, "st.dev"],
             "mean (blade length)" = estimates.stats[[2]][, "mean"],
             "sd (blade length)"   = estimates.stats[[2]][, "st.dev"],
             "mean (fruit length)" = estimates.stats[[3]][, "mean"],
             "sd (fruit length)"   = estimates.stats[[3]][, "st.dev"],
             check.names = FALSE
             )
write.csv(est.stats.table,
          file = "output/gapfilling_test_sparse_estimates_stats_table.csv",
          eol = "\r\n",
          row.names = FALSE
          )

# As above, statistics for the differences:
diff.means <- llply(all.differences, numcolwise(mean))
diff.sds <- llply(all.differences, numcolwise(sd))
differences.stats <- llply(as.list(seq_along(all.differences)),
                         function(i) {
                           data.frame(id = colnames(diff.means[[i]]),
                                      mean = as.numeric(diff.means[[i]]),
                                      st.dev = as.numeric(diff.sds[[i]])
                                      )
                         }
                         )
names(differences.stats) <- paste0(trait.names, ".differences")


# ANOVA of trait estimates
# ------------------------
# We can use ANOVA to test if the estimate means differ.
# To do anova, we have to transform the data to long-list format.
# This is called stacking. We'll use melt() from package reshape2
stacked.estimates <-
  llply(all.estimates,
        function(df) {
          melt(df,
               id.vars = "species",
               variable.name = "filling.method",
               value.name = "trait.value"
               )
        }
        )
names(stacked.estimates) <- paste0(trait.names, ".estimates")

# Do ANOVA: trait.value ~ filling.method
anova.estimates <- llply(stacked.estimates,
                   function(df) {
                     aov(trait.value ~ filling.method, data = df)
                   }
                   )
names(anova.estimates) <- paste0(trait.names, ".anova")
# inspect results with summary(), e.g. summary(anova.estimates[[1]])

# ANOVA of estimate differences
# -----------------------------
# We can use ANOVA to test if the difference means differ.
# To do anova, we have to transform the data to long-list format.
stacked.differences <-
  llply(all.differences,
        function(df) {
          melt(df,
               id.vars = "species",
               variable.name = "filling.method",
               value.name = "trait.value"
               )
        }
        )
names(stacked.differences) <- paste0(trait.names, ".differences")

# Do ANOVA: estimate.diff ~ filling.method
anova.differences <- llply(stacked.differences,
                   function(df) {
                     aov(trait.value ~ filling.method, data = df)
                   }
                   )
names(anova.differences) <- paste0(trait.names, ".anova")
# inspect results with summary(), e.g. summary(anova.differences[[1]])
# Not exactly the same results as the estimates, but close enough.


# Comparing estimated with original values
# -----------------------------------------
# For each trait, a MultiScatter of estimated ~ original, for all methods
y.names <-
  names(all.estimates[[1]]) %>%
  .[!. %in% c("species", "original")]
for (df in seq_along(all.estimates)) {
  for (y in y.names) {
    GraphSVG(MultiScatter(all.estimates[[df]][, "original"],
                          all.estimates[[df]][, y],
                          x.name = paste("Original",
                                         trait.names[df]
                                         ),
                          y.name = paste0("Est. ",
                                          trait.names[df],
                                          " (",
                                          y,
                                          ")"
                                          )
                          ),
             file=paste0("graphs/gapfilling_test_sparse_scatter_",
                         trait.names[df],
                         "_original_vs_",
                         y,
                         ".svg"
                         ),
             width = 12,
             height = 4
             )
  }
}
# Subsequently, a table of the R^2 of estimates ~ original values, for all
# combinations of trait and method.
rsq.table <- matrix(ncol = 3,
                    nrow = length(y.names),
                    dimnames = list(y.names, trait.names)
                    )
for (df in seq_along(all.estimates)) {
  for (i in seq_along(y.names)) {
    mod <- lm(paste(y.names[i], "~", "original"), data = all.estimates[[df]])
    rsq.table[i, df] <- summary(mod)$r.squared
  }
}
write.csv(rsq.table,
          file = "output/gapfilling_test_sparse_est.vs.original.rsq.csv",
          eol = "\r\n",
          row.names = TRUE
          )

cat("Done.\n")

