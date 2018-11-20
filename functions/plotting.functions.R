############################################
# Palm FD project: custom plotting functions
############################################
#
# In which custom plotting functions are defined.
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
# Scatterplot


source(file="functions/base.functions.R")


Scatterplot <- function(x, y, grid = TRUE, pch = 21, col = "black") {
# Wrapper function of plot() for making scatterplots.
# Helps creating scatterplots by automatically implementing good default
# parameters and plotting options. Makes it simple to quickly generate
# scatterplots of decent quality for testing purposes.
#
# Args:
#   x: vector with x-axis values. Can be a list or dataframe of length 1.
#   y: vector with y-axis values. Can be a list or dataframe of length 1.
#   grid: Logical indicating whether to draw grid lines in the plotting area.
#   pch: plotting symbol for points (see ?par for details). Default is 21
#        (hollow circles)
#   col: Color for plotted points. Default is black.
#
# Returns:
#
  if (!IsOneDimensional(x)) {
    stop("argument x is not one-dimensional")
  }
  if (!IsOneDimensional(y)) {
    stop("argument y is not one-dimensional") 
  }
  data <- data.frame(x=x, y=y)
  if (length(which(complete.cases(data))) != length(data[, 1])) {
    data <- data[complete.cases(data), ]
    warning("Omitted observations with missing values")
  }
 # set margins to be nicely narrow
 par(mar=c(4.1, 4.1, .5, .5))
 # Initiate plot, but do not plot values yet
 plot(y ~ x, data=data, type="n", ylim=c(0, 1.02 * max(y)),
      xlab= if (length(names(x)) > 0) {
              names(x)
            } else {
              "x"
            }
      , ylab= if (length(names(y)) > 0) {
                names(y)
              } else {
                "y"
              }
      )
  # Set background color for plotting region. Color is very light grey.
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col=rgb(220, 220, 220, max=255)
       )
  # Add gridlines and zero lines. Grid color is grey.
  if (isTRUE(grid)) {
    grid(nx=NULL, ny=NULL, lty=1, lwd=1, col=rgb(150, 150, 150, max=255))
    abline(h=0, lty=2, col="black")
    abline(v=0, lty=2, col="black")
  }
  # Plot points
  points(y ~ x, data=data, pch=pch, cex=1.0, col=col, bg="transparent")
  # Draw a clean box around the plotting region
  box(lty=1, col="black")
  return (0)
}

