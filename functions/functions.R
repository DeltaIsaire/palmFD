###################################
# Palm FD project: custom functions
###################################
#
# In which custom functions are defined.
# Input files:
#   none
# Generated output files:
#   none

# F01: crossReference
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
  if(is.list(y)) {  # note: is.list returns TRUE for data frames
    if(!length(y) == 1) {
      stop("argument y has more than one element")
    }
  } else {
    if(!length(as.data.frame(y)) == 1) {
      stop("argument y is not one-dimensional")
    }
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

