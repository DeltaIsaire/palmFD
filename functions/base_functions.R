###################################
# Palm FD project: custom functions
###################################
#
# In which custom general-purpose functions are defined.
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
# IsOneDimensional
# CrossCheck
# MultiCheck
# GapFill
# CountObserved
# IntervalRnorm
# StochasticFill
# DetectOS
# WhichMax
# ParseData
# AutoVIF
# GetFD
# RealmSubset


library(plyr)
library(magrittr)
library(car)


IsOneDimensional <- function(x) {
# Checks whether x has more than one element or column.
#
# Args:
#   x: object to be evaluated
#
# Returns: TRUE if x has only one element or column, FALSE if x has more than
#          one element or column
  if (is.list(x)) {  # is.list() returns TRUE for data frames
    if (!length(x) == 1) {  # for data frames, length() returns # of columns
      return (FALSE)
    }
  } else {
    if (!length(as.data.frame(x)) == 1) {
      return (FALSE)
    }
  }
  TRUE
}


CrossCheck <- function(x, y, presence = TRUE, value = TRUE, unique = TRUE) {
# For each value in x, check if it is present or absent in y, then return
# the subset of x values that were present/absent, or their indices.
# 
# Args:
#   x: object whose values should be cross-referenced with the values in y
#   y: object against which values in x are compared
#   presence: Logical indicating whether to check for presence (TRUE) or
#             absence (FALSE) of the values of x in y. Default is TRUE.
#   value: Logical indicating whether to return the values of x (TRUE) which
#          are present/absent in y, or their index in x (FALSE). Default is TRUE.
#   unique: Logical. If and only if 'value = TRUE', subset the returned
#           vector of present/absent values to unique values (TRUE), or 
#           return duplicates as well (FALSE).
#
# Returns:
#  Vector of the same type as x containing the subset of x present/absent in y,
#  or an integer vector with the indices of x whose values are present/absent
#  in y.
  if (is.data.frame(x)) {
    x %<>% as.matrix(.)
  }
  if (is.data.frame(y)) {
    y %<>% as.matrix(.)
  }
  if (isTRUE(presence)) {
    if (isTRUE(value)) {
      if (isTRUE(unique)) {
        return (unique(x[x %in% y]))
      } else {
        return (x[x %in% y])
      }
    } else {
      return (which(x %in% y))
    }
  } else {
    if (isTRUE(value)) {
      if (isTRUE(unique)) {
        return (unique(x[!x %in% y]))
      } else {
        return (x[!x %in% y])
      }
    } else {
      return (which(!x%in% y))
    }
  }
}


MultiCheck <- function(x = NULL, list = NULL, ...) {
# Wrapper function for CrossCheck. Takes a list and cross-references all
# elements of the list with each other.
#
# Args:
#   x: Object whose values are compared against all elements of 'list'.
#      Can be missing, in which case the first element of 'list' is
#      compared to the other elements.
#   list: A list with elements to compare x against.
#   ... : Additional arguments to pass to function 'CrossCheck'
#
# Returns:
#   A vector with the values that are present/absent in all elements of of
#   the provided list.
  if (!missing(x) && !missing(list)) {
    result <- x
  }
  if (missing(x) && !missing(list)) {
    result <- list[[1]]
  }
  if (!missing(x) && missing(list)) {
    if (!is.list(x)) {
      stop("Input is not a list")
    } else {
    list <- x
    result <- list[[1]]
    }
  }
  if (missing(x) && missing(list)) {
    stop("This function cannot run without input")
  }
  iterations <- 
    seq_along(list) %>%
    .[-max(.)]
  for (i in iterations) {
  result <- CrossCheck(x = result, y = list[[i + 1]], ...)
  }
  result
}


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


CountObserved <- function(x) {
# A simple function: find the number of values in x that is NOT missing
#
# Args:
#   x: any object whose elements will be checked for being NA
#
# Returns:
#   A single integer value with the number of non-missing values in x.
  sum(!is.na(x))
}


IntervalRnorm <- function(n, mean = 0, sd = 1, min, max) {
# From a normal distribution with given mean and sd, draw n values that fall within
# the interval given by the min and max.
#
# Args:
#   n: Number of random values sto draw
#   mean: mean of the normal distribution to draw from
#   sd: standard deviation of the normal distribution to draw from
#   min: lower limit for randomly drawn values
#   max: upper limit for randomly drawn values
#
# Returns:
#   A numeric vector of length n with randomly drawn values.
  output <- numeric(length = n)
  for (i in seq_along(output)) {
    do.sample <- TRUE
    while (do.sample) {
      value <- rnorm(1, mean = mean, sd = sd)
      if (value >= min & value <= max) {
        do.sample <- FALSE
      }
    }
    output[i] <- value
  }
  output
}


StochasticFill <- function(trait, genus, distribution) {
# For every missing value in 'trait', get the genus-level trait distribution from
# 'distribution' using 'genus' as index, and sample a random value
# via IntervalRnorm().
#
# Args:
#   trait: numeric vector giving trait data, with missing values
#   genus: character vector of the same length as trait, giving for each trait
#          value the genus to which that trait value belongs.
#   distribution: dataframe giving trait distribution data for each genus.
#                 must have columns with names in 
#                 c("genus", "mean", "sd, "min", "max")
#
# Returns:
#   The trait vector with all the NA values filled in.
  if (!identical(length(trait), length(genus))) {
    stop("lengths of 'trait' and 'genus' differ")
  }
  output <- trait
  for (i in seq_along(trait)) {
    if (is.na(trait[i])) {
      dist <- distribution[match(genus[i], distribution[, "genus"]), ]
      output[i] <- IntervalRnorm(1,
                                 mean = dist[, "mean"],
                                 sd = dist[, "sd"],
                                 min = dist[, "min"],
                                 max = dist[, "max"]
                                 )
    }
  }
  output
}


DetectOS <- function() {
# Detect the operating system for this R session, which can be useful for choosing 
# OS-specific methods. 
# Function code based on an R-bloggers post by "will":
# https://www.r-bloggers.com/identifying-the-os-from-r/
#
# Returns:
#   character vector of length 1 with value in c("osx", "linux", "windows")
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)[[1]]
}


WhichMax <- function(x, n = 1) {
# Function to find the indices of the n highest values in a vector.
#
# Args:
#   x: numeric vector in which to find the n highest values
#   n: number of highest values to find
# Returns:
#   The n indices in x corresponding to the highest values in x, in order of
#   highest to lowest value in x.
  values <- sort(x, decreasing = TRUE)[seq_len(n)]
  CrossCheck(x = x, y = values, presence = TRUE, value = FALSE)
}



ParseData <- function(x, response, standardize = TRUE, numeric.only = FALSE,
                      double.std = FALSE) {
# Function to subset a dataframe to complete cases, and optionally standardize
# all numeric variables EXCEPT the response variable.
# Unless 'double.std' is TRUE, in which case the response variable is also
# standardized. 
  if (numeric.only) {
  x.num <- numcolwise(function(x) { x } ) (x)
  x.complete <- x.num[complete.cases(x.num), ]
  } else {
    x.complete <- x[complete.cases(x), ]
  }
  if (standardize) {
    x.std <- numcolwise(scale, center = TRUE, scale = TRUE) (x.complete)
    # Re-add non-numeric columns:
    missing.cols <- names(x.complete)[!names(x.complete) %in% names(x.std)]
    x.std[, missing.cols] <- x.complete[, missing.cols]
    if (!isTRUE(double.std)) {
      # un-standardize the response variable:
      x.std[, response] <- x.complete[, response]
      x.complete <- x.std
    }
  }
  # move response variable to last column
  resp.ind <- which(names(x.complete) == response)
  if (! resp.ind == ncol(x.complete)) {
    x.complete <- x.complete[, c(1:(resp.ind - 1), (resp.ind + 1):ncol(x), resp.ind)]
  }
  x.complete
}


AutoVIF <- function(x, response, standardize = TRUE, threshold = 5,
                    numeric.only = FALSE) {
# Function to automatically exclude predictors based on VIF. An OLS model is
# fitted with all predictors, then the predictor with highest VIF is removed. This
# is repeated until the highest remaining VIF is below the threshold.
#
# In the case of categorical predictor variables, instead of VIF we compute
# GVIF^[1/(2*df)], square it, and then interpret it the same as regular VIF scores.
# see ?vif and Fox & Monette 1990 (https://www.tandfonline.com/doi/abs/10.1080/01621459.1992.10475190#.U2jkTFdMzTo)
#
# Args:
#   x: data.frame with all variables (response + predictors)
#   response: character vector giving column name of the response variable.
#   standardize: Should predictors be standardized to mean 0 and unit variance?
#   threshold: VIF value deemed to be the upper acceptable limit. Traditionally,
#     people have used a threshold of 10, but Zuur et al (2010) reccommends a
#     threshold of 3, or even lower when signal is weak.
#   numeric.only: Subset dataset to only continuous numerical predictors?
# Returns:
#   A character vector with the names of remaining PREDICTOR variables

  # Validate data: subset to complete, numeric data and standardize
  x.complete <- ParseData(x = x,
                          response = response,
                          standardize = standardize,
                          numeric.only = numeric.only
                          )

  # Iteratively exclude most collinear predictor.
  # First construct initial mod formula
  predictors <- colnames(x.complete)[-match(response, colnames(x.complete))]
  form <- paste(response, "~", predictors[1])
  for (i in seq_len(length(predictors) - 1)) {
    form <- paste(form, "+", predictors[i + 1])
  }
  form %<>% as.formula()
  # Then apply loop for the selection
  do.exclude <- TRUE
  while (do.exclude) {
    # Identify most collinear predictor
    if (length(predictors) > 1) {
      mod <- lm(form, data = x.complete)
      scores <- vif(mod)
      if (!IsOneDimensional(scores)) {
        scores <- scores[, "GVIF^(1/(2*Df))"] ^ 2
      }
    } else {
      scores <- 0
    }
    # IF highest VIF is above threshold, construct new formula
    if (max(scores) > threshold) {
      predictors <- predictors[-match(names(scores)[which.max(scores)], predictors)]
      form <- paste(response, "~", predictors[1])
      for (i in seq_len(length(predictors) - 1)) {
        form <- paste(form, "+", predictors[i + 1])
      }
      form %<>% as.formula()
    } else {
      do.exclude <- FALSE
    }
  }

  # Return names of predictors not excluded
  predictors
}

# Example code
if (FALSE) {
  fric <- read.csv(file = "output/FD_summary_FRic.csv", row.names = 1)
  env <- read.csv(file = "output/tdwg3_environmental_predictors.csv", row.names = 1)
  # These are all approximately normally distributed numerical predictors.
  x <- env
  x$FRic <- fric$FRic.global.SES
  response <- "FRic"
  x <- x[, c(1:9, 12:24, 27:28, 29)]
  # Using bio7, instead of bio5 + bio6
  # using LGM climate anomaly, instead of LGM climate

  a <- AutoVIF(x, response, threshold = 5)
  b <- AutoVIF(x, response, threshold = 3)
}



# Function to extract data of each null model from fd.indices
# -----------------------------------------------------------
GetFD <- function(fd.indices, name) {
#  fd.indices: list of FD index dataframes (generated above)
#  name: name of column to extract from each df in 'fd.indices'
  fd.list <- llply(fd.indices, function(x) { x[, name] } )
  df <- as.data.frame(fd.list)
  rownames(df) <- rownames(fd.indices[[1]])
  df
}


# Function to subset the input data to realm, combining results into a list
# -------------------------------------------------------------------------
RealmSubset <- function(x) {
# Where x is a dataframe, such as the output of GetFD(), with tdwg3 codes as
# rownames
  output <- vector("list", length = length(realm.tdwg3))
  names(output) <- names(realm.tdwg3)
  for (i in seq_along(output)) {
    output[[i]] <- x[rownames(x) %in% realm.tdwg3[[i]], ]
  }
  output
}

