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


library(plyr)
library(magrittr)


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

