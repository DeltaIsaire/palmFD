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
crossReference <- function(x,y,value = TRUE) {
# description
# 
# Args:
#   x: Vector whose values should be cross-referenced with the values in y
#   y: Vector against which values in x are compared
#   value: if TRUE, return the values of x missing from y. If FALSE, then index
#          numbers of x are returned instead. Default is TRUE.
#
# Returns:
#  Vector of values (or indices) of x that are missing from y.
  a <- which(
    is.na(
      as.numeric(
        lapply(x, function(x) {
          which(x == y)
          }
    ))))
  if (value == TRUE) {
    print(x[a])
  } else {
    print(a)
  }
}

