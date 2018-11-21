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
# GraphPDF
# MultiScatter


source(file="functions/base.functions.R")


Scatterplot <- function(x, y, grid = TRUE, pch = 21, col = "black",
                        xlab = NULL, ylab = NULL) {
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
#   A new plot frame with the plot. In addition, Scatterplot sets custom
#   values for par(mar) without reverting them afterwards.
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
 plot(y ~ x, data=data, type="n", 
      ylim= if (min(y) < 0) {
              c(1.02 * min(y), 1.02 * max(y))
            } else {      
              c(0, 1.02 * max(y))
            }
      , xlab= if (!is.null(xlab)) {
                xlab
              } else {
                if (length(names(x)) > 0) {
                  names(x)
                } else {
                  "x"
                }
              }
      , ylab= if (!is.null(ylab)) {
                ylab
              } else {
                if (length(names(y)) > 0) {
                  names(y)
                } else {
                  "y"
                }
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
}
# TODO: add formula notation option.
# to check for formula, use inherits(object, "formula") which returns TRUE or
# FALSE. 


GraphPDF <- function(expr, file, ...) {
# Wrapper function for 'pdf', to simplify saving graphs as pdf files.
#
# Args:
#   expr: Function which creates a graph
#   file: A character string giving the name of the file to create.
#   ...: additional arguments to pass to the 'pdf' function
#
# Returns:
#   Creates a PDF files with the given filename, containing the graph.
  pdf(file=file, ...)
  eval.parent(substitute(expr))
  on.exit(dev.off())
}


MultiScatter <- function(x, y, x.name="x", y.name="y", ...) {
# Produce a neatly formatted multi-frame graph. IN DEVELOPMENT
#
# Args:
#   x: vector with x-axis values. Can be a list or dataframe of length 1.
#   y: vector with y-axis values. Can be a list or dataframe of length 1.
#   x.name: character vector of length 1 giving x-axis label. Will be used in
#            all frames.
#   y.name: character vector of length 1 giving a name for the y-values.
#            Will be used in the y-axis labels.
#   ...: additional arguments to pass to Scatterplot function. These arguments
#        will apply to all frames.

# Output:
#   A graph with 2 frames
  par(mfrow=c(1, 2))
  # Genus mean vs BHPMF:
  Scatterplot(x=x, y=y, xlab=x.name, ylab=y.name, ...)
  # Add frame label in top left:
  text(x=par("usr")[1], y=par("usr")[4], 
       labels="A", cex=2, adj=c(-0.5, 1.5))
  # Add 1:1 line for reference:
  lines(x=c(par("usr")[1], par("usr")[2]), y=c(par("usr")[3], par("usr")[4]))
  # Second frame: a simple kind of residuals plot, which is visually easier
  # to interpret.
  Scatterplot(x=x, y=(y-x), xlab=x.name, ylab=paste0("(", x.name, " - ",
                                                     y.name, ")"), ...)
  # Add frame label in top right:
  text(x=par("usr")[2], y=par("usr")[4], 
       labels="B", cex=2, adj=c(1.2, 1.5))
  # Add solid horizontal line for comparison:
  abline(h=0, lty=1, col="black")
  # return par to default:
  par(mfrow=c(1, 1))
}

