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
# AllSARSingles
# ParseSARSingle
# FitGlobalSAR
# NbPlot

library(magrittr)
library(plyr)
library(spdep)
library(ncf)

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
    stop("mod object must be of class 'sarlm'")
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




AllSARSingles <- function(fd.indices, colname, predictors, tdwg.map,
                          dist.weight = FALSE, name.all, ...) {
# Function to run all models via RunMSSAR
# Args:
#   fd.indices: the list object 'fd.indices' or a subset thereof
#   colname: name of column to extract from each df in 'fd.indices' via GetFD()
#   predictors: dataframe with predictor variables
#   tdwg.map: the map object as created by read_sf(). Is automatically
#             subsetted as required. Used for creating the spatial weights matrix.
#   dist.weight: whether to apply neighbour distance weighting in the spatial
#                weights matrix. See SWMat()
#   name.all: 'name' argument passed to RunMSSAR(). Should be of length 1; suffixes
#         will be added for the cases 'full' and each realm subset
#   ...: Additional arguments to pass to RunMSSAR()

  fd.all <- GetFD(fd.indices, colname)

  # full dataset:
  cat(name.all, "full...\n")
  ind.complete <- complete.cases(fd.all) & complete.cases(predictors)
  tdwg.map.subset <- tdwg.map[ind.complete, ]
  nb <- SoiNB(tdwg.map.subset)
  swmat <- SWMat(nb, tdwg.map.subset, dist.weight = dist.weight, style = "W")
  RunMSSAR(name = paste0(name.all, "_full"),
           responses = fd.all,
           predictors = predictors,
           listw = swmat,
           ...
           )

  # realm subsets:
  fd.all.realms <- RealmSubset(fd.all)
  predictors.realms <- RealmSubset(predictors[, !colnames(predictors) %in% "realm"])
  for (i in seq_along(fd.all.realms)) {
    cat(name.all, names(fd.all.realms)[i], "...\n")
    ind.complete <-
      complete.cases(fd.all.realms[[i]]) & complete.cases(predictors.realms[[i]])
    tdwg.ind <- tdwg.map$LEVEL_3_CO %in% rownames(fd.all.realms[[i]])[ind.complete]
    tdwg.map.subset <- tdwg.map[tdwg.ind, ]
    nb <- SoiNB(tdwg.map.subset)
    swmat.realm <- SWMat(nb, tdwg.map.subset, dist.weight = dist.weight, style = "W")
    RunMSSAR(name = paste0(name.all, "_", names(fd.all.realms)[i]),
             responses = fd.all.realms[[i]],
             predictors = predictors.realms[[i]],
             listw = swmat.realm,
             ...
             )
  }
  return (0)
}



ParseSARSingle <- function(name.all, cases, statistics = c("slope", "p.value",
                                                          "moran.p", "pseudo.Rsq"),
                           filter.p = TRUE) {
# Automatically read the files created by AllSARSingles() and combine the data
# into a dataframe.
#
# filenames generated by ALLSARSingles() have the following structure:
#   paste0(name.all, "_", case, ".csv")
#   Where 'name.all' is the name argument passed to ALLSarSingles() and 'case' is
#   either "full" or one of the realm names.
#
# Args:
#   name.all: the 'name.all' argument previously passed to ALLSarSingles()
#   cases: cases to return, vector with values "full" and the realm names, or a
#          subset thereof
#   statistics: model statistics to parse, in c("slope", "p.value", "moran.p",
#                                             "pseudo.Rsq")
#              See SingleSAR(). If more than one statistic is provided, the output
#              will be a list with the dataframe for each statistic
#   filter.p: whether to replace statistics from models where p > 0.05 with NA
  stat.list <- vector("list", length = length(statistics))
  names(stat.list) <- statistics
  for (i in seq_along(stat.list)) {
    output <- vector("list", length = length(cases))
    names(output) <- cases
    for (j in seq_along(output)) {
      # load requested statistic
      stat <- read.csv(file = paste0(name.all,
                                     "_",
                                     cases[j],
                                     "_",
                                     statistics[i],
                                     ".csv"
                                     )
                       )
        if (filter.p) {
          # load corresponding p-values
          p.vals <- read.csv(file = paste0(name.all, "_", cases[j], "_p.value.csv"))
          # filter data for p < 0.05
          for (row in seq_len(nrow(stat))) {
            nonsignificant <- suppressWarnings(which(p.vals[row, ] > 0.05))
           # warnings suppressed because the 1st column is not numeric (originally
           # rownames) which is not in fact a problem
            stat[row, nonsignificant] <- "n.s."
          }
        }
      output[[j]] <- stat
    }
    # Collapse output into a dataframe. 
    output <- ldply(output, function(x) { x } , .id = "realm.subset")
    stat.list[[i]] <- output
  }
  if (length(stat.list) == 1) {
    stat.list <- stat.list[[1]]
  }
  stat.list
}



FitGlobalSAR <- function(predictors, response, tdwg.map, dist.weight = TRUE) {
# Function to fit a SAR error model to the provided data, using a Sphere of
# Influence neighbourhood based on the provided spatial data.
# Provided data is subsetted to complete cases via ParseData().
# NOTE: this function uses global assignments, for compatibility with dredge() from
# package MuMIn.
# Args:
#   response: vector of response data
#   predictors: dataframe with predictor variables, in the same order as the response
#               data.
#   tdwg.map: simple features polygon data
#   dist.weight: whether to apply neighbour distance weighting in the spatial
#                weights matrix. See SWMat()
  sar.mod.data <<- 
    data.frame(predictors, response = response) %>%
    ParseData(., response = "response", standardize = TRUE, numeric.only = TRUE)

  # Generate spatial weights matrix from SOI neighbourhood
  ind.complete <- complete.cases(predictors) & complete.cases(response)
  tdwg.map.subset <- tdwg.map[which(ind.complete), ]
  nb <- SoiNB(tdwg.map.subset)
  sar.swmat <<- SWMat(nb, tdwg.map.subset, dist.weight = dist.weight, style = "W")

  # Fit global model
  sar.mod.formula <<-
    paste0("response ~ ", paste(colnames(predictors), collapse = " + ")) %>%
    as.formula()
  errorsarlm(formula = sar.mod.formula,
             data = sar.mod.data,
             listw = sar.swmat,
             na.action = na.fail
             )
}



NbPlot <- function(nb, tdwg.map, tdwg.map.subset, title = NULL, subtitle = NULL,
                   filename = NULL) {
# Function to plot the Sphere of Influence neighbourhood for a particular dataset.
#
# Args:
#   nb: A neighbourhood object, such as the output of SoiNB()
#   tdwg.map: the tdwg.map loaded with geom_sf()
#   tdwg.map.subset: the (subsetted) tdwg.map object used to create the neighbourhood
#                    (e.g. the object used in the SoiNB() call)
#   title: title for plot.
#   subtitle: subtitle for plot
#   filename: filename for saving the graph (file type by extension, see ggsave()).
#             Can be NULL, in which case the graph is NOT saved.
#
#  Returns:
#    The neighbourhood graph as a ggplot object
  segments <- Nb2Segments(nb, tdwg.map.subset)
  presence <- tdwg.map$LEVEL_3_CO %in% tdwg.map.subset$LEVEL_3_CO
  presence[!presence] <- NA
  soi.plot <- SpatialPlotNB(tdwg.map[!tdwg.map$LEVEL_3_CO == "ANT", ],
                            presence[!tdwg.map$LEVEL_3_CO == "ANT"],
                            segments,
                            title = title,
                            subtitle = subtitle
                          )
  if (!is.null(filename)) {
    ggsave(soi.plot, filename = filename, width = 8, height = 4)
  }
  invisible(soi.plot)
}

