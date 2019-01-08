####################################
# Test gapfilling code for coherence
####################################
#
# In which I test if my gapfilling code is working as intended.
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


# -----------------------------------------
# Data preparation: artificial trait matrix
# -----------------------------------------
# Ten genera (A-J) with 3 species each.
# One trait, with values of increasing integers, i.e.
# genus A has species with values c(1, 2, 3), genus B with c(4, 5, 6), etc.
test.matrix <- data.frame(species = as.character(1:30),
                          genus = sort(rep(LETTERS[1:10], 3)),
                          trait = 1:30
                          )
# Create artificial gaps for say 10 out of 30 values.
# 1.3 is the magic number, to compensate for unique() excluding some indices.
indices <- 
  runif(n = 10 * 1.3, min = 1, max = length(test.matrix$trait)) %>%
  round() %>%
  unique()
test.matrix.sparse <- test.matrix
test.matrix.sparse[indices, "trait"] <- NA

# ---------------------------------
# Test gap-filling with genus means
# ---------------------------------
# All code here is a direct copy of the code I've been using.

# Step 1. Calculate genus-means for both complete and sparse test matrix
means.complete <- ddply(test.matrix, "genus", numcolwise(mean, na.rm = TRUE))
# The means should be
2 + 3 * 0:9
# and they are.
means.sparse <- ddply(test.matrix.sparse, "genus", numcolwise(mean, na.rm = TRUE))
# The same code should work again. I checked: it does.

# Step 2. Verify gap-filling.
# The code I have been using:
trait.names <- "trait"
filled <- GapFill(test.matrix.sparse,
                  means.sparse,
                  by = "genus",
                  fill = trait.names
                  )
# Now to verify. Which indices were NA?
missing.indices <- which(is.na(test.matrix.sparse$trait))
# To which genera do these species belong?
genera <- as.character(test.matrix.sparse[missing.indices, "genus"])
# What were the estimated means for these genera?
est.means <- numeric(length = length(genera))
for (i in 1:length(genera)) {
  est.means[i] <- means.sparse[which(means.sparse$genus == genera[i]), "trait"]
}
# Which values were used to fill the trait matrix?
filled.values <- filled[missing.indices, "trait"]
# And the big test:
test <- identical(est.means, filled.values)
cat(test, "\n")
# [1] TRUE
# Now you can rerun this script as many times as you want
# source("test_gapfilling_code.R")
# The outcome doesn't change. 
# MY CODE IS FINE.

