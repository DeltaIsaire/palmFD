####################################
# Test gapfilling code for coherence
####################################
#
# Having established that the genus-mean gap-filling code works,
# we need to track down the observed negative correlation between
# 'observed' trait value and 'estimated' trait value.
#


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

# Normally I do:
# source(file = "functions/base_functions.R")
# For the sake of portability, I'll just paste my GapFill() function here:

GapFill <- function(x, y, by, fill) {
# Fills missing values in the selected columns of x using the data provided
# in y.
#
# Args:
#   x: Dataframe with columns to be gap-filled
#   y: Dataframe with data to use for gap filling. Must contain columns with
#      names matching the columns in x which are to be gap-filled.
#   by: Character vector of length 1 giving the column name of the column
#       used to match x to y. Both dataframes must contain a column with this name.
#   fill: Character vector with the column names in x that should be gap-filled.
#
# Returns:
#   Dataframe x where the selected columns are gap-filled
  if (!is.data.frame(x)) {
    stop("argument 'x' must be a dataframe")
  }
  if (!is.data.frame(y)) {
    stop("argument 'y' must be a dataframe")
  }
  if (!(is.character(by) && length(by) == 1)) {
    stop("argument 'by' must be a character vector of length 1")
  }
  if (!is.character(fill)) {
    stop("argument 'fill' must be a character vector")
  }
  x.by <- which(names(x) == by)
  y.by <- which(names(y) == by)
  # Robustness against cases where the 'by' column is a factor with differing
  # levels between dataframes:
  if (is.factor(x[, x.by])) {
    x[, x.by] %<>% as.character()
    convert <- TRUE
  }
  if (is.factor(y[, y.by])) {
    y[, y.by] %<>% as.character()
  }
  filled <- 
    lapply(fill, 
           function(name) {
           missing <- x[, x.by][which(is.na(x[, name]))]
           indices <- as.numeric(lapply(missing, 
                                        function(x) {
                                          which(y[, y.by] == x)
                                        }
                                        )      
                                 )
           x[, name][which(is.na(x[, name]))] <- y[, name][indices]
           x[, name]
         }
         ) %>%
    simplify2array() %>%
    data.frame()
  names(filled) <- fill
  fill.indices <- 
    lapply(fill, function(a) { which(names(x) == a) } ) %>%
    simplify2array()
  x[, fill.indices] <- filled
  if (convert) {
    x[, x.by] %<>% as.factor()
  }
  x
}


# ---------------------
# differing genus means
# ---------------------
# We make a large artificial test dataset.
# 100 genera with 20 species each.
# one trait, with values from a normal distribution with mean = 1:100
# and sd = 1.
# Notice that each genus has a different mean.
trait.vals <- numeric(length = 100 * 20)
for (i in 1:100) {
  trait.vals[(20 * i - 19):(20 * i)] <- rnorm(n = 20, mean = i, sd = 1)
}
test.matrix <- data.frame(species = as.character(1:(100 * 20)),
                          genus = as.character(sort(rep(1:100, 20))),
                          trait = trait.vals
                          )

# Create artificial gaps for say 30% of values.
# 1.3 is the magic number, to compensate for unique() excluding some indices.
indices <- 
  runif(n = 0.3 * length(test.matrix$trait) * 1.3,
        min = 1,
        max = length(test.matrix$trait)
        ) %>%
  round() %>%
  unique()
test.matrix.sparse <- test.matrix
test.matrix.sparse[indices, "trait"] <- NA

# gap-fill the sparse test dataset with genus means
means.sparse <- ddply(test.matrix.sparse, "genus", numcolwise(mean, na.rm = TRUE))
trait.names <- "trait"
filled <- GapFill(test.matrix.sparse,
                  means.sparse,
                  by = "genus",
                  fill = trait.names
                  )

# Scatterplot of original vs estimated values
# Which indices were NA?
missing.indices <- which(is.na(test.matrix.sparse$trait))
# Extract original values
originals <- test.matrix[missing.indices, "trait"]
# extract estimates
estimates <- filled[missing.indices, "trait"]

plot(x = originals,
     y = estimates,
     type = "p",
     pch = 21,
     col = "black",
     xlab = "Original trait value",
     ylab = "Estimate (genus mean)"
     )
# Add 1:1 line for reference:
lines(x = c(par("usr")[1], par("usr")[2]), 
      y = c(par("usr")[1], par("usr")[2])
      )
# This looks fine.

# Could the problem be my plotting function?
#source(file = "functions/plotting_functions.R")
#MultiScatter(x = originals,
#             y = estimates,
#             x.name = "Original",
#             y.name = "Estimate"
#             )
# Nope, that's the exact same thing.
# (Code commented out for portability)


# -----------------
# common genus mean
# -----------------
# Here, we will use trait data from a single normal distribution.
trait.vals <- rnorm(n = 100 * 20, mean = 10, sd = 1)

test.matrix <- data.frame(species = as.character(1:(100 * 20)),
                          genus = as.character(sort(rep(1:100, 20))),
                          trait = trait.vals
                          )

# Create artificial gaps for say 30% of values
# Remember 1.3 is the magic number.
indices <- 
  runif(n = 0.3 * length(test.matrix$trait) * 1.3,
        min = 1,
        max = length(test.matrix$trait)
        ) %>%
  round() %>%
  unique()
test.matrix.sparse <- test.matrix
test.matrix.sparse[indices, "trait"] <- NA

# gap-fill the sparse test dataset with genus means
means.sparse <- ddply(test.matrix.sparse, "genus", numcolwise(mean, na.rm = TRUE))
trait.names <- "trait"
filled <- GapFill(test.matrix.sparse,
                  means.sparse,
                  by = "genus",
                  fill = trait.names
                  )

# Scatterplot of original vs estimated values
# Which indices were NA?
missing.indices <- which(is.na(test.matrix.sparse$trait))
# Extract original values
originals <- test.matrix[missing.indices, "trait"]
# extract estimates
estimates <- filled[missing.indices, "trait"]

plot(x = originals,
     y = estimates,
     type = "p",
     pch = 21,
     col = "black",
     xlab = "Original trait value",
     ylab = "Estimate (genus mean)",
     ylim = c(0, 1.1 * max(estimates))
     )
# Add 1:1 line for reference:
lines(x = c(par("usr")[1], par("usr")[2]), 
      y = c(par("usr")[1], par("usr")[2])
      )

# Gotcha. Notice how the error is horizontal, while the x:y line is not.
# My multiscatter plot shows this well:
#MultiScatter(x = originals,
#             y = estimates,
#             x.name = "Original",
#             y.name = "Estimate"
#             )

# EXPLANATION:
# I was right all along, only my supervisors didn't notice :P
# The estimate is the mean. If a 'real' value is above the mean, then
# the estimate is too low. Similarly, if the 'real' value is below
# the mean, then the estimate is too high.
# In general, the error is (mean - x) where x is the original value.
# So we can expect the negative correlation between observation and estimate
# to follow a line with slope -1, which crosses the x-axis
# for x = genus mean.
#
# Let's make a fancy graph to visualize this.


# -----------------------------
# Visualizing gap-filling error
# -----------------------------
# Prepare a test dataset with 2 genera, 100 species each.
# We should vary both the mean and sd.
trait.vals <- numeric(length = 100 * 2)
trait.vals[1:100] <- rnorm(n = 100, mean = 2, sd = 0.5)
trait.vals[101:200] <- rnorm(n = 100, mean = 4, sd = 2)

test.matrix <- data.frame(species = as.character(1:(100 * 2)),
                          genus = rep(LETTERS[1:2], 100) %>%
                                  sort() %>%
                                  as.character(),
                          trait = trait.vals
                          )

# Create artificial gaps for say 30% of values
# Remember 1.3 is the magic number.
indices <- 
  runif(n = 0.3 * length(test.matrix$trait) * 1.3,
        min = 1,
        max = length(test.matrix$trait)
        ) %>%
  round() %>%
  unique()
test.matrix.sparse <- test.matrix
test.matrix.sparse[indices, "trait"] <- NA

# gap-fill the sparse test dataset with genus means
means.sparse <- ddply(test.matrix.sparse, "genus", numcolwise(mean, na.rm = TRUE))
trait.names <- "trait"
filled <- GapFill(test.matrix.sparse,
                  means.sparse,
                  by = "genus",
                  fill = trait.names
                  )

# Extract the interesting bits
missing.indices <- which(is.na(test.matrix.sparse$trait))
originals <- test.matrix[missing.indices, "trait"]
estimates <- filled[missing.indices, "trait"]

# plot the data.
# For portability, I'll keep it simple.
ErrorPlot <- function() {
  par(mfrow = c(1, 2))
  plot(x = originals,
       y = estimates,
       type = "p",
       pch = 21,
       xlab = "Original trait value",
       ylab = "Estimate (genus mean)",
       ylim = c(0, 1.1 * max(estimates))
     )
lines(x = c(par("usr")[1], par("usr")[2]), 
      y = c(par("usr")[1], par("usr")[2])
      )
  plot(x = originals,
       y = (estimates - originals),
       type = "p",
       pch = 21,
       xlab = "Original trait value",
       ylab = "(Mean estimate - original)"
       )
abline(h = 0)
}
ErrorPlot()
# Speaks for itself, really.

# I'll save a good-quality graph:
#source(file = "functions/plotting_functions.R")
#GraphSVG(MultiScatter(x = originals,
#                      y = estimates,
#                      x.name = "Original",
#                      y.name = "Mean estimate"
#                      ),
#         file = "graphs/test_gapfilling_error_genus_mean.svg",
#         width = 12,
#         height = 4
#         )

