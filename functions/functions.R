###################################
# Palm FD project: custom functions
###################################
#
# In which custom functions are defined.
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
# isOneDimensional
# crossCheck
# multiCC
# gapFill

library(plyr)


isOneDimensional <- function(x) {
# Checks whether x has more than one element or column.
#
# Args:
#   x: object to be evaluated
#
# Returns: TRUE if x has only one element or column, FALSE if x has more than
#          one element or column
  if (is.list(x)) {  # note: is.list returns TRUE for data frames
    if (!length(x) == 1) {
      return(FALSE)
    }
  } else {
    if (!length(as.data.frame(x)) == 1) {
      return(FALSE)
    }
  }
  return(TRUE)
}


crossCheck <- function(x, y, presence = TRUE, value = TRUE) {
# For each value in x, check if it is present or absent in y, then return
# the subset of x values that were present/absent, or their indices.
# 
# Args:
#   x: Vector whose values should be cross-referenced with the values in y
#   y: Vector against which values in x are compared
#   presence: Logical indicating whether to check for presence (TRUE) or
#             absence (FALSE) of the values of x in y. Default is TRUE.
#   value: Logical indicating whether to return the values of x (TRUE) which
#          are present/absent in y, or their index (FALSE). Default is TRUE.
#
# Returns:
#  Vector of the same type as x containing the subset of x present/absent in y,
#  or an integer vector with the indices of x whose values are present/absent
#  in y.
  if (!isOneDimensional(x)) {
    stop("argument x is not one-dimensional")
  }
  if (!isOneDimensional(y)) {
    stop("argument y is not one-dimensional")
  }
  checklist <- lapply(x, function(x) {
                 which(x == y)
                 }
               )
  if (isTRUE(presence)) {
    present <- lapply(checklist, function(x) {
                 length(which(!is.na(x)))
                 }
               )
    if (isTRUE(value)) {
      return(x[which(present >= 1)])
    } else {
      return(which(present >= 1))
    }
  } else {
    absent <- lapply(checklist, function(x) {
                which(length(x) == 0)
                }
              )
    if (isTRUE(value)) {
      return(x[which(absent >= 1)])
    } else {
      return(which(absent >= 1))
    }
  }
}


multiCC <- function(x, presence = TRUE, value = TRUE) {
# Wrapper function for crossCheck. Takes a list and cross-references all
# elements of the list with each other.
#
# Args:
#   x: List with elements to cross-reference.
#   presence: Logical indicating whether to check for presence (TRUE) or
#             absence (FALSE) of the values being cross-referenced.
#              Default is TRUE.
#   value: Logical (see crossCheck)
#
# Returns:
#   When value = TRUE, returns a list of the values present/absent in each
#   pairwise comparison. When value = FALSE, returns a matrix listing the
#   number of values missing in each pairwise comparison.
  if (isOneDimensional(x)) {
    stop("Argument x is one-dimensional")
  }
  if (isTRUE(value)) {
    a <- plyr::llply(x, function(a) {
           plyr::llply(x, crossCheck, y=a, presence=presence, value=TRUE)
         }
         )
    return(a)
  } else {
    a <- plyr::llply(x, function(a) { 
           plyr::llply(x, crossCheck, y=a, presence=presence, value=FALSE)
         }
         )
    b <- plyr::laply(simplify2array(a), length)
    c <- matrix(data=b, nrow=length(x), ncol=length(x),byrow=FALSE,
      dimnames=list(x=names(x), y=names(x)))
    return(c)
  }
}


gapFill <- function(x, y, by, fill) {
# description
#
# Args:
#   x: Dataframe with columns to be gap-filled
#   y: Dataframe with data to use for gap filling. Must contain columns with
#      names matching the columns in x which are to be gap-filled.
#   by: Character vector of length 1 giving the column name used to match x
#       to y. Both dataframes must contain a column with this name.
#   fill: Character vector with the column names in x that should be gap-filled
#
# Returns:
#   Dataframe x where the selected columns are gap-filled
  if (!is.data.frame(x)) {
    stop("argument 'x' must be a dataframe")
  }
  if (!is.data.frame(y)) {
    stop("argument 'y' must be a dataframe")
  }
  if (! (is.character(by) && length(by) == 1)) {
    stop("argument 'by' must be a character vector of length 1")
  }
  if (!is.character(fill)) {
    stop("argument 'fill' must be a character vector")
  }
  x.by <- which(names(x) == by)
  y.by <- which(names(y) == by)
  name = fill[1]  # temporary for development
  missing <- x[, x.by][which(is.na(x[, name]))]
  indices <- unlist(lapply(missing, function(x) {
               which(y[, y.by] == x)
             }
             ))
  x[, name][which(is.na(x[, name]))] <- y[, name][indices]
  return(x)
}
