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
# crossReference
# multiCrossRef

library(plyr)


isOneDimensional <- function(x) {
# Checks whether x has more than one element or column.
#
# Args:
#   x: object to be evaluated
#
# Returns: TRUE if x has only one element or column, FALSE if x has more than
#          one element or column
  if(is.list(x)) {  # note: is.list returns TRUE for data frames
    if(!length(x) == 1) {
      return(FALSE)
    }
  } else {
    if(!length(as.data.frame(x)) == 1) {
      return(FALSE)
    }
  }
  return(TRUE)
}


crossReference <- function(x, y, presence = TRUE, value = TRUE) {
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
  if(!isOneDimensional(x)) {
    stop("argument x is not one-dimensional")
  }
  if(!isOneDimensional(y)) {
    stop("argument y is not one-dimensional")
  }
  a <- as.numeric(
         lapply(x, function(x) {
           which(x == y)
          }
       ) )
  if (isTRUE(presence)) {
    present <- which(!is.na(a))
    if (isTRUE(value)) {
      return(x[present])
    } else {
      return(present)
    }
  } else {
    absent <- which(is.na(a))
    if (isTRUE(value)) {
      return(x[absent])
    } else {
      return(absent)
    }
  }
}


multiCrossRef <- function(x, presence = TRUE) {
# Wrapper function for crossReference. Takes a list and cross-references all
# elements of the list with each other.
#
# Args:
#   x: List with elements to cross-reference.
#   presence: Logical indicating whether to check for presence (TRUE) or
#             absence (FALSE) of the values being cross-refernced.
#              Default is TRUE.
#
# Returns:
#   A dataframe listing the number of values missing in each pairwise comparison.
  if(isOneDimensional(x)) {
    stop("Argument x is one-dimensional")
  }
  a <- plyr::llply(x, function(a) { 
       plyr::llply(x, crossReference, y=a, presence=presence, value=FALSE)
       }
       )
  b <- plyr::laply(simplify2array(a), length)
  c <- matrix(data=b, nrow=length(x), ncol=length(x),byrow=FALSE,
    dimnames=list(names(x), names(x)))
  return(data.frame(c))
}

