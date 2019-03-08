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
#   output/test/palm_traits.csv
#   data/PalmTraits_10.csv
#   output/test/palm_hierarchy.csv
# Generated output files:
#   output/test/test_sparse_traits.csv
#   output/test/test_sparse_traits_complete.csv
#   output/test/test_sparse_filled_genus_mean.csv
#   output/test/test_BHPMF_one_mean.txt
#   output/test/SPARSE_test_BHPMF_one_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_one.csv
#   output/test/SPARSE_test_BHPMF_two_mean.txt
#   output/test/SPARSE_test_BHPMF_two_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_two.csv
#   output/test/SPARSE_test_BHPMF_three_mean.txt
#   output/test/SPARSE_test_BHPMF_three_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_three.csv
#   output/test/SPARSE_test_BHPMF_four_mean.txt
#   output/test/SPARSE_test_BHPMF_four_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_four.csv
#   output/test/SPARSE_test_BHPMF_five_mean.txt
#   output/test/SPARSE_test_BHPMF_five_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_five.csv
#   output/test/SPARSE_test_BHPMF_six_mean.txt
#   output/test/SPARSE_test_BHPMF_six_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_six.csv
#   output/test/SPARSE_test_BHPMF_seven_mean.txt
#   output/test/SPARSE_test_BHPMF_seven_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_seven.csv
#   output/test/SPARSE_test_BHPMF_eight_mean.txt
#   output/test/SPARSE_test_BHPMF_eight_std.txt
#   output/test/SPARSE_traits_filled_BHPMF_eight.csv
#
#   graphs/test/test_sparse_scatter_ < a bunch of graphs starting with this name >
#
#   output/test/test_sparse_estimates_stats_table.csv


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
palm.traits <- read.csv(file = "output/test/palm_traits.csv")
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
          file = "output/test/test_sparse_traits.csv",
          eol = "\r\n",
          row.names = FALSE
          )
write.csv(complete.traits,
          file = "output/test/test_sparse_traits_complete.csv",
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
          file = "output/test/test_sparse_filled_genus_mean.csv",
          eol="\r\n",
          row.names=FALSE
          )


# ----------------------
# Gap-filling with BHPMF
# ----------------------
cat("(2) with BHPMF methods:\n")

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
  read.csv(file = "output/test/palm_hierarchy.csv") %>%
  .[CrossCheck(.$species,
               sparse.traits$species,
               presence = TRUE,
               value = FALSE
               ),
    ] %>%
  as.matrix()
test.hierarchy <- hierarchy.matrix

TestBHPMF <- function(trial) {
# Standard code to run a BHPMF test trial.
# Trial: character vector of length 1 giving the name of the trial,
#        e.g. "one".
  # BHPMF wants a preprocessing directory, where it saves pre-processed files.
  # To avoid errors or erronous output when re-running the code, this directory
  # needs to be emptied.
  cat("Running BHPMF trial", trial, "\n")
  unlink("output/test/BHPMF_preprocessing_test", recursive=TRUE)
  dir.create("output/test/BHPMF_preprocessing_test")
  GapFilling(X = test.matrix,
             hierarchy.info = test.hierarchy,
             prediction.level = 4,
             used.num.hierarchy.levels = 3,
             mean.gap.filled.output.path = paste0("output/test/SPARSE_test_BHPMF_",
                                                  trial,
                                                  "_mean.txt"
                                                  ),
             std.gap.filled.output.path = paste0("output/test/SPARSE_test_BHPMF_",
                                                 trial,
                                                 "_std.txt"
                                                 ),
             tmp.dir = "output/test/BHPMF_preprocessing_test", 
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
    read.table(file = paste0("output/test/SPARSE_test_BHPMF_",
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
    GapFill(sparse.traits,
            .,
            by = "species",
            fill = trait.names
            ) %>%
    { .[complete.cases(.), ] }
  
  write.csv(traits.filled.BHPMF.test,
            file = paste0("output/test/SPARSE_traits_filled_BHPMF_",
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

trial <- "one"
TestBHPMF(trial)
ParseBHPMF(trial)

# 2. Gap-filling with BHPMF: using all estimates
# ----------------------------------------------
cat("(2) Using estimates for all trait values\n")

trial <- "two"
TestBHPMF(trial)
# Not using ParseBHPMF() because here we take all estimates.
traits.filled.BHPMF.test <-
  read.table(file = paste0("output/test/SPARSE_test_BHPMF_",
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
          file = paste0("output/test/SPARSE_traits_filled_BHPMF_",
                        trial,
                        ".csv"
                        ),
          eol="\r\n",
          row.names=FALSE
          )

# 3. Gap-filling with BHPMF: including 1 dummy trait
# --------------------------------------------------
cat("(3) including 1 dummy trait\n")

# Add dummy trait:
test.matrix <-
  trait.matrix %>%
  data.frame(.,
             dummy = rep(1, times = length(.[, 1]))
             ) %>%
  as.matrix()

trial <- "three"
TestBHPMF(trial)
ParseBHPMF(trial)
# Undo adding dummy trait:
test.matrix <- trait.matrix

# 4. Gap-filling with BHPMF: including more traits
# ------------------------------------------------
cat("(4) including additional palm traits\n")

indices <- CrossCheck(x = trait.data$SpecName,
                      y = rownames(trait.matrix),
                      presence = TRUE,
                      value = FALSE
                      )
test.matrix <- 
  cbind(trait.matrix,
        trait.data[indices,
                   c("MaxStemDia_cm",
                     "MaxLeafNumber",
                     "Max_Rachis_Length_m",
                     "Max_Petiole_length_m",
                     "AverageFruitWidth_cm"
                     )
                   ]
        ) %>%
  as.matrix()

trial <- "four"
TestBHPMF(trial)
ParseBHPMF(trial)
# Undo adding additional traits:
test.matrix <- trait.matrix

# 5. Gap-filling with BHPMF: including 10 dummy traits
# ----------------------------------------------------
cat("(5) including 10 dummy traits\n")

# Add ten dummy traits:
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
# Undo adding dummy traits:
test.matrix <- trait.matrix

# 6. Gap-filling with BHPMF: Including growthform in hierarchy
# ------------------------------------------------------------
cat("(6) including growthform as hierarchy level\n")

# Construct categorical variable of growthform and add it to the hierarchy
growthform <- numeric(length = length(trait.data$SpecName))
growthform[which(trait.data$Climbing == 1)] <- "climbing"
growthform[which(trait.data$Acaulescence == 1)] <- "acaulescent"
growthform[which(trait.data$Errect == 1)] <- "freestanding"
growthform[which(growthform == 0)] <- NA

test.hierarchy <- 
  read.csv(file = "output/test/palm_hierarchy.csv") %>%
  { data.frame(species = .$species,
             growthform = growthform,
             .[, -1]
             )
  }
test.hierarchy$growthform %<>% as.character()  # friggin' factors :rolleyes:
test.hierarchy <- 
  ddply(test.hierarchy,
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
  .[CrossCheck(.$species,
               sparse.traits$species,
               presence = TRUE,
               value = FALSE
               ),
    ] %>%
  as.matrix()
# Remove species with NA for growthform.
to.remove <- which(is.na(test.hierarchy[, 2]))
test.matrix %<>% .[-to.remove, ]
test.hierarchy %<>% . [-to.remove, ]
rm(to.remove)

trial <- "six"
# Not using TestBHPMF() because we use a different hierarchy
cat("Running BHPMF trial", trial, "\n")
unlink("output/test/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/test/BHPMF_preprocessing_test")
GapFilling(X = test.matrix,
           hierarchy.info = test.hierarchy,
           prediction.level = 5,
           used.num.hierarchy.levels = 4,
           mean.gap.filled.output.path = paste0("output/test/SPARSE_test_BHPMF_",
                                                trial,
                                                "_mean.txt"
                                                ),
           std.gap.filled.output.path = paste0("output/test/SPARSE_test_BHPMF_",
                                               trial,
                                               "_std.txt"
                                               ),
           tmp.dir = "output/test/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE,
           verbose=verbose
           )
ParseBHPMF(trial)
# Undo changing hierarchy matrix and trait matrix:
test.hierarchy <- hierarchy.matrix
test.matrix <- trait.matrix

# 7 + 8. Gap-filling with BHPMF: including even more traits
# ---------------------------------------------------------
# Extract and clean trait data. Binary traits should get values 1 (for 0)
# and 100 (for 1).
indices <- CrossCheck(x = trait.data$SpecName,
                      y = rownames(trait.matrix),
                      presence = TRUE,
                      value = FALSE
                      )
extra.traits <- trait.data[indices,
                           c("MaxStemDia_cm",
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

test.matrix.all <- cbind(test.matrix, as.matrix(extra.traits))

test.matrix.growthform <- test.matrix.all[, c(1:3, 9:11)]

cat("(7) three continuous traits + three growthform binary traits\n")
test.matrix <- test.matrix.growthform

trial <- "seven"
TestBHPMF(trial)
ParseBHPMF(trial)
# undo edit test matrix:
test.matrix <- trait.matrix

cat("(8) all usable traits\n")
test.matrix <- test.matrix.all

trial <- "eight"
TestBHPMF(trial)
ParseBHPMF(trial)
# undo edit test matrix:
test.matrix <- trait.matrix


# -----------------------------
# Parsing gap-filled trait data
# -----------------------------
cat("Parsing gap-filled trait matrices...\n")

# Load filled trait matrices:
#   genus means
filled.mean <- read.csv(file = "output/test/test_sparse_filled_genus_mean.csv")
#   BHPMF, in different ways
Load <- function(trial) {
  read.csv(file = paste0("output/test/SPARSE_traits_filled_BHPMF_",
                         trial,
                         ".csv"
                         )
           )
}
all.filled <- list(original = complete.traits,
                   mean     = filled.mean,
                   b.one    = Load("one"),
                   b.two    = Load("two"),
                   b.three  = Load("three"),
                   b.four   = Load("four"),
                   b.five   = Load("five"),
                   b.six    = Load("six"),
                   b.seven  = Load("seven"),
                   b.eight  = Load("eight")
                   )
filled.names <- names(all.filled)
# growthform.BHPMF (b.six) is incomplete (contains info for fewer species).
# For an honest comparison, we must subset the other dataframes to these
# same species
all.filled <- llply(all.filled,
                    function(df) {
                      indices <- CrossCheck(x = df$species,
                                            y = all.filled$b.six$species,
                                            presence = TRUE,
                                            value = FALSE
                                            )
                      df[indices, ]
                    }
                    )
# And don't forget to also subset sparse.traits:
indices <- CrossCheck(x = sparse.traits$species,
                      y = all.filled$b.six$species,
                      presence = TRUE,
                      value = FALSE
                      )
sparse.traits %<>% .[indices, ]

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
          file = "output/test/test_sparse_estimates_stats_table.csv",
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
             file=paste0("graphs/test/test_sparse_scatter_",
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
          file = "output/test/test_sparse_est.vs.original.rsq.csv",
          eol = "\r\n",
          row.names = TRUE
          )



cat("Done.\n")

