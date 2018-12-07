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
# Generated output files:
#   output/test_sparse_traits.csv
#   output/test_sparse_filled_genus_mean.csv
#   output/test_sparse_std_BHPMF_mean.txt
#   output/test_sparse_std_BHPMF_std.txt
#   output/test_sparse_filled_std_BHPMF.csv
#   output/test_sparse_dummy_BHPMF_mean.txt
#   output/test_sparse_dummy_BHPMF_std.txt
#   output/test_sparse_filled_dummy_BHPMF.csv
#   graphs/test_sparse_scatter_ <a bunch of graphs starting with this name>


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(BHPMF)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


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
          file = "output/test_sparse_traits.csv",
          eol="\r\n",
          row.names=FALSE
          )


# -----------------------------------
# Gap-filling the sparse trait matrix
# -----------------------------------
cat("Gap-filling the sparse test matrix:\n")

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
          file = "output/test_sparse_filled_genus_mean.csv",
          eol="\r\n",
          row.names=FALSE
          )

# Gap-filling with standard BHPMF
# -------------------------------
cat("(2) with standard BHPMF:\n")
cat("Preparing data...\n")
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

hierarchy.matrix <- 
  read.csv(file = "output/palm_hierarchy.csv") %>%
  .[CrossCheck(.$species,
               sparse.traits$species,
               presence = TRUE,
               value = FALSE
               ),
    ] %>%
  as.matrix()
cat("Running BHPMF... (this may take a moment)\n")
# BHPMF wants a preprocessing directory, where it saves pre-processed files.
# To avoid errors or erronous output when re-running the code, this directory
# needs to be emptied.
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X = trait.matrix,
           hierarchy.info = hierarchy.matrix,
           prediction.level = 4,
           used.num.hierarchy.levels = 3,
           mean.gap.filled.output.path = "output/test_sparse_std_BHPMF_mean.txt",
           std.gap.filled.output.path = "output/test_sparse_std_BHPMF_std.txt",
           tmp.dir = "output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE,
           verbose=FALSE
           )
traits.filled.std.BHPMF <-
  read.table(file = "output/test_sparse_std_BHPMF_mean.txt",
             header = TRUE,
             sep = "	"
             ) %>%
  data.frame(species = as.character(hierarchy.matrix[, 1]),
             genus   = as.character(hierarchy.matrix[, 2]),
             .
             ) %>%
  GapFill(sparse.traits,
          .,
          by = "species",
          fill = trait.names
          ) %>%
  { .[complete.cases(.), ] }

write.csv(traits.filled.std.BHPMF,
          file = "output/test_sparse_filled_std_BHPMF.csv",
          eol="\r\n",
          row.names=FALSE
          )

# Gap-filling with dummy BHPMF
# -------------------------------
cat("(3) with dummy BHPMF:\n")
cat("Preparing data...\n")
# Source dataframes are sorted by species name, so rows will correspond
# automatically.
trait.matrix <- 
  sparse.traits[, trait.names] %>%
  as.matrix()
rownames(trait.matrix) <- sparse.traits$species
trait.matrix %<>%
  data.frame(.,
             dummy = rep(1, times = length(.[, 1]))
             ) %>%
  as.matrix()
# BHPMF is unable to handle observed trait values of exactly 0.
# These do occur in our dataset, for stem height.
# Solution: set these values to something close to 0
trait.matrix[which(trait.matrix == 0)] <- 0.0001

hierarchy.matrix <- 
  read.csv(file = "output/palm_hierarchy.csv") %>%
  .[CrossCheck(.$species,
               sparse.traits$species,
               presence = TRUE,
               value = FALSE
               ),
    ] %>%
  as.matrix()
cat("Running BHPMF... (this may take a moment)\n")
# BHPMF wants a preprocessing directory, where it saves pre-processed files.
# To avoid errors or erronous output when re-running the code, this directory
# needs to be emptied.
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X = trait.matrix,
           hierarchy.info = hierarchy.matrix,
           prediction.level = 4,
           used.num.hierarchy.levels = 3,
           mean.gap.filled.output.path = "output/test_sparse_dummy_BHPMF_mean.txt",
           std.gap.filled.output.path = "output/test_sparse_dummy_BHPMF_std.txt",
           tmp.dir = "output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE,
           verbose=FALSE
           )
traits.filled.dummy.BHPMF <-
  read.table(file = "output/test_sparse_dummy_BHPMF_mean.txt",
             header = TRUE,
             sep = "	"
             ) %>%
  data.frame(species = as.character(hierarchy.matrix[, 1]),
             genus   = as.character(hierarchy.matrix[, 2]),
             .
             ) %>%
  GapFill(sparse.traits,
          .,
          by = "species",
          fill = trait.names
          ) %>%
  { .[complete.cases(.), ] }

write.csv(traits.filled.dummy.BHPMF,
          file = "output/test_sparse_filled_dummy_BHPMF.csv",
          eol="\r\n",
          row.names=FALSE
          )

# Gap-filling with growthform BHPMF
# ---------------------------------
cat("(4) with growthform BHPMF:\n")
cat("Preparing data...\n")
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

# Construct hierarchy matrix with growthform:
trait.data <- 
  read.csv(file = "data/PalmTraits_10.csv") %>%
  .[order(.$SpecName), ]
trait.data$SpecName %<>% sub(pattern = " ",
                             replacement = "_",
                             x = .
                             )
growthform <- numeric(length = length(trait.data$SpecName))
growthform[which(trait.data$Climbing == 1)] <- "climbing"
growthform[which(trait.data$Acaulescence == 1)] <- "acaulescent"
growthform[which(trait.data$Errect == 1)] <- "freestanding"
growthform[which(growthform == 0)] <- NA

hierarchy.matrix <- 
  read.csv(file = "output/palm_hierarchy.csv") %>%
  { data.frame(species = .$species,
             growthform = growthform,
             .[, -1]
             )
  }
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
  .[CrossCheck(.$species,
               sparse.traits$species,
               presence = TRUE,
               value = FALSE
               ),
    ] %>%
  as.matrix()
# Remove species with NA for growthform.
to.remove <- which(is.na(hierarchy.matrix[, 2]))
trait.matrix %<>% .[-to.remove, ]
hierarchy.matrix %<>% . [-to.remove, ]
rm(to.remove)

cat("Running BHPMF... (this may take a moment)\n")
# BHPMF wants a preprocessing directory, where it saves pre-processed files.
# To avoid errors or erronous output when re-running the code, this directory
# needs to be emptied.
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
GapFilling(X = trait.matrix,
           hierarchy.info = hierarchy.matrix,
           prediction.level = 5,
           used.num.hierarchy.levels = 4,
           mean.gap.filled.output.path = 
             "output/test_sparse_growthform_BHPMF_mean.txt",
           std.gap.filled.output.path =
             "output/test_sparse_growthform_BHPMF_std.txt",
           tmp.dir = "output/BHPMF_preprocessing_test", 
           rmse.plot.test.data=FALSE,
           verbose=FALSE
           )
traits.filled.growthform.BHPMF <-
  read.table(file = "output/test_sparse_growthform_BHPMF_mean.txt",
             header = TRUE,
             sep = "	"
             ) %>%
  data.frame(species = as.character(hierarchy.matrix[, 1]),
             genus   = as.character(hierarchy.matrix[, 2]),
             .
             ) %>%
  GapFill(sparse.traits,
          .,
          by = "species",
          fill = trait.names
          ) %>%
  { .[complete.cases(.), ] }

write.csv(traits.filled.growthform.BHPMF,
          file = "output/test_sparse_filled_growthform_BHPMF.csv",
          eol="\r\n",
          row.names=FALSE
          )


# -----------------------------
# Parsing gap-filled trait data
# -----------------------------
cat("Parsing gap-filled trait matrices...\n")
all.filled <- list(original = complete.traits,
                   mean = traits.filled.mean,
                   std.BHPMF = traits.filled.std.BHPMF,
                   dummy.BHPMF = traits.filled.dummy.BHPMF,
                   growthform.BHPMF = traits.filled.growthform.BHPMF
                   )
filled.names <- names(all.filled)
# growthform.BHPMF is incomplete (contains info for fewer species).
# For an honest comparison, we must subset the other dataframes to these
# same species
all.filled <- llply(all.filled,
                    function(df) {
                      indices <- CrossCheck(x = df$species,
                                            y = all.filled[[5]]$species,
                                            presence = TRUE,
                                            value = FALSE
                                            )
                      df[indices, ]
                    }
                    )
# And don't forget to also subset sparse.traits:
indices <- CrossCheck(x = sparse.traits$species,
                      y = all.filled[[5]]$species,
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
          data.frame(species = df$species,
                     mean = df$mean - df$original,
                     std.BHPMF = df$std.BHPMF - df$original,
                     dummy.BHPMF = df$dummy.BHPMF - df$original,
                     growthform.BHPMF = df$growthform.BHPMF - df$original
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
transformed.estimates <-
  llply(all.estimates,
        function(df) {
          data.frame(trait.value = c(df$original,
                                     df$mean,
                                     df$std.BHPMF,
                                     df$dummy.BHPMF,
                                     df$growthform.BHPMF
                                     ),
                     filling.method = c(rep("original", dim(df)[1]),
                                        rep("genus.mean", dim(df)[1]),
                                        rep("std.BHPMF", dim(df)[1]),
                                        rep("dummy.BHPMF", dim(df)[1]),
                                        rep("growthform.BHPMF", dim(df)[1])
                                        )
                     )
        }
        )
names(transformed.estimates) <- paste0(trait.names, ".estimates")

# Do ANOVA: trait.value ~ filling.method
anova.estimates <- llply(transformed.estimates,
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
transformed.differences <-
  llply(all.differences,
        function(df) {
          data.frame(estimate.diff = c(df$mean,
                                     df$std.BHPMF,
                                     df$dummy.BHPMF,
                                     df$growthform.BHPMF
                                     ),
                     filling.method = c(rep("genus.mean", dim(df)[1]),
                                        rep("std.BHPMF", dim(df)[1]),
                                        rep("dummy.BHPMF", dim(df)[1]),
                                        rep("growthform.BHPMF", dim(df)[1])
                                        )
                     )
        }
        )
names(transformed.differences) <- paste0(trait.names, ".differences")

# Do ANOVA: estimate.diff ~ filling.method
anova.differences <- llply(transformed.estimates,
                   function(df) {
                     aov(trait.value ~ filling.method, data = df)
                   }
                   )
names(anova.differences) <- paste0(trait.names, ".anova")
# inspect results with summary(), e.g. summary(anova.differences[[1]])
# OBSERVATION: THESE ANOVA RESULTS ARE *EXACTLY* THE SAME AS FOR THE ESTIMATES.

# Scatterplots comparing estimated with original values
# -----------------------------------------------------
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
             file=paste0("graphs/test_sparse_scatter_",
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

# Scatterplots of standard BHPMF vs dummy BHPMF estimates
# -------------------------------------------------------
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "std.BHPMF"],
                        all.estimates[[df]][, "dummy.BHPMF"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "std.BHPMF"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "dummy.BHPMF"
                                       )
                        ),
             file=paste0("graphs/test_sparse_scatter_",
                         trait.names[df],
                         "_std.BHPMF_vs_dummy.BHPMF.svg"
                         ),
             width = 12,
             height = 4
             )
}

# Scatterplots of genus mean vs standard BHPMF estimates
# ------------------------------------------------------
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "mean"],
                        all.estimates[[df]][, "std.BHPMF"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "mean"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "std.BHPMF"
                                       )
                        ),
             file=paste0("graphs/test_sparse_scatter_",
                         trait.names[df],
                         "_genus_mean_vs_std.BHPMF.svg"
                         ),
             width = 12,
             height = 4
             )
}

# Scatterplots of standard BHPMF vs growthform BHPMF estimates
# ------------------------------------------------------------
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "std.BHPMF"],
                        all.estimates[[df]][, "growthform.BHPMF"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "std.BHPMF"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "growthform.BHPMF"
                                       )
                        ),
             file=paste0("graphs/test_sparse_scatter_",
                         trait.names[df],
                         "_std.BHPMF_vs_growthform.BHPMF.svg"
                         ),
             width = 12,
             height = 4
             )
}

# Scatterplots of original vs estimates for BHPMF: standard + growthform
# ----------------------------------------------------------------------
# For each trait, plot the original values vs the estimates from standard BHPMF
# and growthform BHPMF. 
OneScatter <- function(index) {
  Scatterplot(x = all.estimates[[index]][, "original"],
              y = all.estimates[[index]][, "std.BHPMF"],
              xlab = paste("Original",
                           trait.names[index]
                           ),
              ylab = paste("Estimated",
                           trait.names[index]
                           ),
              pch = 17
              )
    points(x = all.estimates[[index]][, "original"],
           y = all.estimates[[index]][, "growthform.BHPMF"],
           col = "red",
           pch = 19
           )
  # Add 1:1 line for reference:
  lines(x = c(par("usr")[1], par("usr")[2]), 
        y = c(par("usr")[1], par("usr")[2])
        )
  # legend
  legend("bottomright",
         legend = c("Standard BHPMF", "Growthform BHPMF"),
         col = c("black", "red"),
         pch = c(17, 19),
         bg = "white"
         )
}

for (i in seq_along(all.estimates)) {
  GraphSVG(OneScatter(i),
           file = paste0("graphs/test_sparse_scatter_",
                         trait.names[i],
                         "_original_vs_std_and_growthform_BHPMF.svg"
                         ),
           width = 6,
           height = 4
           )
}


cat("Done.\n")

# TODO: compare uncertainties of BHPMF, i.e. the st.dev outputs.
# TODO: compare accuracy of estimates for species with NA for all real traits.
