################################################
# Palm FD Project: identifying gap-filling error
################################################
#
# In the gap-filling accuracy test we have observed a negative correlation between
# the 'observed' trait value and the 'estimated' trait value.
# Here, I explore the origin of this pattern by applying genus-mean gapfilling to
# a completely artificial dataset.
#
# Input files:
#   none
#
# Generated output files:
#   graphs/gapfilling_error_one.svg
#   graphs/gapfilling_error_two.svg
#   graphs/gapfilling_error_three.svg


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")

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

FirstPlot <- function() {
  Scatterplot(x = originals,
              y = estimates,
              xlab = "Original trait value",
              ylab = "Estimate (genus mean)"
              )
  # Add 1:1 line for reference:
  lines(x = c(par("usr")[1], par("usr")[2]), 
        y = c(par("usr")[1], par("usr")[2])
        )
}
GraphSVG(FirstPlot(),
         file = "graphs/gapfilling_error_one.svg",
         width = 6,
         height = 6
         )
# This looks fine.


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

SecondPlot <- function() {
  Scatterplot(x = originals,
              y = estimates,
              xlab = "Original trait value",
              ylab = "Estimate (genus mean)",
              )
  # Add 1:1 line for reference:
  lines(x = c(par("usr")[1], par("usr")[2]), 
        y = c(par("usr")[1], par("usr")[2])
        )
}
GraphSVG(SecondPlot(),
         file = "graphs/gapfilling_error_two.svg",
         width = 6,
         height = 6
         )
# Gotcha. Notice how the error is horizontal, while the x:y line is not.
#
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
GraphSVG(MultiScatter(x = originals,
                      y = estimates,
                      x.name = "Original",
                      y.name = "Mean estimate"
                      ),
         file = "graphs/gapfilling_error_three.svg",
         width = 12,
         height = 4
         )

cat("Done.\n")

