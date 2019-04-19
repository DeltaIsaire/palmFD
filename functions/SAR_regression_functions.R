###############################################################
# Palm FD project: custom functions for SAR regression modeling
###############################################################
#
# In which functions to help with SAR regression model construction and selection
# are defined. It is assumed that SAR error models are used.
#
# Input files:
#   none
# Generated output files:
#   none
#
# FUNCTION LIST:
# PseudoRsq
# SingleSAR
# MultiSingleSAR
# RunMSSAR
# SoiNB
# SWMat
# Nb2Segments

library(magrittr)
library(plyr)
library(spdep)

source(file = "functions/base_functions.R")



PseudoRsq <- function(sarlm) {
  # No such thing as a real Rsq for SAR models, because SAR is likelihood-based.
  # However, Kissling et al (2012a) use Pseudo-Rsq values, based on the method of
  # Kissling & Carl (2008), defined as "squared Pearson correlation between
  # predicted and observed values."
  # HOWEVER: if you simply use fitted() on a SAR model object, this pseudo-rsq
  # will include the influence of the spatial autocorrelation term,
  # not only the effect of the tested predictor. If the spatial structure is
  # strong, pseudo-Rsq can be high even if the included predictor is not
  # significant.
  # To see this effect, you could include a 'null' predictor with random values.
  # Here, we instead calculate fitted values directly from the estimated intercept
  # and slope(s).
  if (!identical(class(sarlm), "sarlm")) {
    stop("mod object must be of class 'sarlm")
  }
  coeff.mat <- t(t(sarlm[["X"]]) * sarlm[["coefficients"]])
  cor(x = rowSums(coeff.mat), y = sarlm[["y"]], method = "pearson") ^ 2
}



SingleSAR <- function(x, listw, response, standardize = TRUE, digits = "all",
                      numeric.only = FALSE) {
# Fit every single-predictor SAR-error model and report the slope, predictor p-value,
# p-value for spatial autocorrelation of the model residuals via moran.test(),
# and the pseudo-Rsq following the methods of Kissling & Carl (2008).
#
# Args:
#   x: data.frame with all variables (response + predictors)
#   listw: a spatial weights matrix (of class 'listw') for the SAR model
#   response: character vector giving column name of the response variable.
#   standardize: Should predictors be standardized to mean 0 and unit variance?
#     NOTE: if this is FALSE, then the reported slopes are not directly comparable!
#   digits: integer giving number of decimal places to report, or "all" to skip
#     rounding.
#   numeric.only: Subset dataset to only continuous numerical predictors?
# Returns:
#  A matrix with 4 columns (slope, p-value, moran.p, pseudo-Rsq) reporting results
#  for each single predictor.

  # Validate data: subset to complete, numeric data and standardize
  x.complete <- ParseData(x = x,
                          response = response,
                          standardize = standardize,
                          numeric.only = numeric.only
                          )

  # Init output matrix
  y.index <- match(response, colnames(x.complete))
  mat <- matrix(data = NA,
                ncol = 4,
                nrow = ncol(x.complete) - 1,
                dimnames = list(colnames(x.complete)[-y.index],
                                c("slope", "p.value", "moran.p", "pseudo.Rsq")
                                )
                )

  # Fit each single-predictor model and extract results
  for (index in seq_len(ncol(x.complete))[-y.index]) {
    cat("\t\tPredictor:", colnames(x.complete)[index], "\n")
    sar.mod <-
      errorsarlm(as.formula(paste(response, "~", colnames(x.complete)[index])),
                 data = x.complete,
                 listw = listw
                 )
    # slope:
    mat[match(names(x.complete)[index], rownames(mat)), 1] <-
      summary(sar.mod)[["Coef"]] [2, 1]
    # predictor p.value. Sadly unreliable for categorical predictors,
    # so replace those with NA
    mat[match(names(x.complete)[index], rownames(mat)), 2] <-
      summary(sar.mod)[["Coef"]] [2, 4]
    mat[which(as.logical(colwise(is.discrete) (x.complete[, -y.index]))),
        2
        ] <- NA
    # moran.p
    mat[match(names(x.complete)[index], rownames(mat)), 3] <-
      moran.test(resid(sar.mod), listw = listw)[["p.value"]]
    # Pseudo.Rsq
    mat[match(names(x.complete)[index], rownames(mat)), 4] <-
      PseudoRsq(sar.mod)
  }

  # Return results, rounding as specified
  if (!identical(digits, "all")) {
    mat <- round(mat, digits = digits)
  }
  mat
}



MultiSingleSAR <- function(responses, predictors, listw, ...) {
# Wrapper function to run SingleSAR() multiple times. For each response variable,
# merge it to 'predictors' then apply SingleSAR() to the merged dataframe.
#
# Args:
#   responses: data frame with response variables
#   Predictors: data frame with predictor variables
#   ...: additional arguments passed to SingleSAR()

  # init output list
  mat <- matrix(data = NA,
                nrow = ncol(responses),
                ncol = ncol(predictors),
                dimnames = list(colnames(responses), colnames(predictors))
                )
  output <- list(mat, mat, mat, mat)
  names(output) <- c("slope", "p.value", "moran.p", "pseudo.Rsq")

  # Generate single-predictor models for each response variable
  for (i in seq_along(responses)) {
    response <- colnames(responses)[i]
        cat("\tResponse:", response, "\n")
    model.data <- predictors
    model.data[, response] <- responses[, response]
    mat <- SingleSAR(model.data,
                     listw = listw,
                     response = response,
                     ...
                     )
    indices <- match(rownames(mat), colnames(output[[1]]))
    for (j in 1:4) {
      output[[j]] [i, indices] <- mat[, j]
    }
  }
  output
}



RunMSSAR <- function(name, responses, predictors, listw, ...) {
# Wrapper function to run MultiSingleSAR() and save the results to disk.
#
# Args:
#   Name: character string as unique data identifier. Used as prefix for file names.
#   ...: Additional arguments to pass to MultiSingleSAR()

  output <- MultiSingleSAR(responses = responses,
                           predictors = predictors,
                           listw = listw,
                           ...
                           )

  for (i in seq_along(output)) {
    write.csv(output[[i]],
              file = paste0(name, "_", names(output)[i], ".csv"),
              eol = "\r\n"
              )
  }
  output
}



SoiNB <- function(tdwg.map) {
# Function to create the Sphere of Influence neighbourhood. It is defined as follows:
#   1. For each polygon, draw a circle around the centroid with radius = distance
#      to nearest neighbour
#   2. Polygons whose "influence" circles overlap are neighbours.
# Args:
#   tdwg.map: the map object as created by read_sf(), or a subset thereof
  tdwg.spatial <- as(tdwg.map, "Spatial")
  tdwg.coords <- coordinates(tdwg.spatial)
  nb.dlny <- tri2nb(tdwg.coords)
  nb.soi <-
    soi.graph(nb.dlny, tdwg.coords) %>%
    graph2nb()
  nb.soi
}



SWMat <- function(nb, tdwg.map = NULL, longlat = TRUE, dist.weight = FALSE, ...) {
# Function to create spatial weights matrix, optionally with neighbour weighting
# by distance.
# Args:
#   nb: a neighbourhood object
#   tdwg.map: the map object as created by read_sf() that was used to create 'nb'
#   longlat: argument passed to nbdists()
#   dist.weight: whether to use weighting by distance
#   ...: further arguments to pass to nb2listw()

  if (dist.weight) {
    # 1. For each polygon, the (great circle) distances to neighbours:
    tdwg.spatial <- as(tdwg.map, "Spatial")
    nb.dist <- nbdists(nb, coordinates(tdwg.spatial), longlat = longlat)
    # 2. The greatest (great circle) distance between any two neighbours:
    nb.maxdist <- max(unlist(nb.dist))
    # 3. For each polygon, weight neighbour distances by the greatest distance
    nb.wdist <- lapply(nb.dist, function(x) { 1 - x / nb.maxdist } )
    # 4. Create distance-weighted weights matrix:
    swmat <- nb2listw(nb, glist = nb.wdist, ...)
  } else {
    swmat <- nb2listw(nb, ...)
  }
  swmat
}




Nb2Segments <- function(nb, tdwg.map) {
# Function to create segments dataframe from a neighbourhood object
# Args:
#   nb: a neighbourhood object
#   tdwg.map: the map object as created by read_sf() that was used to create 'nb'
  nb.num <- laply(nb, length)
  tdwg3.names <- tdwg.map$LEVEL_3_CO
  segments <-
    matrix(data = NA,
           nrow = sum(nb.num),
           ncol = 4,
           dimnames = list(rep(tdwg3.names, times = nb.num),
                           c("x", "y", "xend", "yend")
                           )
                     )
  segments[, 1] <- tdwg.map$Long[match(rownames(segments), tdwg.map$LEVEL_3_CO)]
  segments[, 2] <- tdwg.map$Lat[match(rownames(segments), tdwg.map$LEVEL_3_CO)]
  indices <- match(tdwg3.names[unlist(nb)], tdwg.map$LEVEL_3_CO)
  segments[, 3] <- tdwg.map$Long[indices]
  segments[, 4] <- tdwg.map$Lat[indices]
  as.data.frame(segments)
}



