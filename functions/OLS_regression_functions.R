###############################################################
# Palm FD project: custom functions for OLS regression modeling
###############################################################
#
# In which functions to help with OLS regression model construction and selection
# are defined.
#
#
# FUNCTION LIST:
# SelectOLS
# SingleOLS
# SingleModels
# MultiSelect
# FitGlobalOLS
# PlotResids
# FDResidPlot

library(magrittr)
library(plyr)
library(leaps)
library(car)

source(file = "functions/base_functions.R")



SelectOLS <- function(x, response, k = 10, standardize = TRUE, method = "forward",
                      numeric.only = FALSE, vif.threshold = 3, double.std = FALSE) {
# Function to perform automated OLS regression model selection for a given dataset.
# The provided data is subsetted to complete cases of all numeric variables. It is 
# then further subsetted to only non-collinear predictors, using function AutoVIF().
# The best OLS model of each length is then found via 'regsubsets' from the 'leaps'
# package. Finally, the optimal length is estimated using cross-validation, adjusted
# Rsq, AIC and BIC.
#
# NOTE: beware of linear dependencies. With bioclim data for instance, bio7 equals
# (bio5 - bio6), which doesn't provide any new information. You must use either 5 + 6
# OR 7, but not all three at once.
#
# Code based on https://uc-r.github.io/model_selection
#
# Args:
#   x: data.frame with all variables (response + predictors)
#   response: character vector giving column name of the response variable.
#   k: number of folds to use for k-fold cross-validation.
#   standardize: Should predictors be standardized to mean 0 and unit variance?
#   method: model selection method to use, see regsubsets()
#   numeric.only: Subset dataset to only continuous numerical predictors?
#   vif.threshold: Threshold of VIF score that is deemed acceptable.
# Returns:
#   A list of three:
#     formulae: A list with the formula for the model selected by each method
#       (cross.validation, adj.rsq, cp, bic)
#     mods: model objects of OLS models fitted with these formulae
#     performance: a matrix reporting the length, Rsq, AIC, BIC and highest VIF of
#       the four selected models

  # Validate data: subset to complete, numeric data and standardize
  x.complete <- ParseData(x = x,
                          response = response,
                          standardize = standardize,
                          numeric.only = numeric.only,
                          double.std = double.std
                          )

  # If categorical predictors are included, ensure appropriate parsing.
  # Convert characters to factors:
  char.inds <- which(sapply(x.complete, is.character))
  x.complete[char.inds] <- lapply(x.complete[char.inds], as.factor)

  # Filter predictor variables using AutoVIF()
  to.keep <- AutoVIF(x.complete, response = response, threshold = vif.threshold)
  x.complete %<>% .[, colnames(.) %in% c(to.keep, response)]

  # If we have factors, argument nvmax to regsubsets() must be adjusted:
  factor.inds <- which(sapply(x.complete, is.factor))
  factor.levels <- lapply(x.complete[factor.inds], levels)
  total.levels <- length(unlist(factor.levels))
  adj.nvmax <- total.levels - 2 * length(factor.inds)

  # Find best models of each length
  fitted.all <- regsubsets(as.formula(paste(response, "~", ".")),
                            data = x.complete,
                            nvmax = ncol(x.complete) + adj.nvmax,
                            method = method
                            )

  # Find optimal length using k-fold cross-validation:
  # Randomly assign each pairwise observation to a fold, in a balanced fashion
  folds <- sample(x = rep(seq_len(k), length.out = nrow(x.complete)),
                  size = nrow(x.complete),
                  replace = FALSE
                  )
  # init matrix for prediction errors
  cv.errors <-
    matrix(data = NA,
           nrow = k,
           ncol = ncol(x.complete) - 1,
           dimnames = list(fold = seq_len(k),
                           mod.length = seq_len(ncol(x.complete) - 1)
                           )
           )
  # Forward selection for each of the k training datasets
  for (i in seq_len(k)) {
    fitted <- regsubsets(as.formula(paste(response, "~", ".")),
                          data = x.complete[!folds == i, ],
                          nvmax = ncol(x.complete) + adj.nvmax,
                          method = method
                          )
    # Predict test values using the best model of each size
    for (j in seq_len(ncol(x.complete) - 1)) {
      my.formula <- as.formula(paste(response, "~", "."))
      mat <- model.matrix(my.formula, x.complete[folds == i, ])
      coeff <- coef(fitted, id = j)
      xvars <- names(coeff)
      predictions <- mat[, xvars] %*% coeff
    # and find the mean square difference with the test data
      cv.errors[i, j] <-
        mean((x.complete[, response][folds == i] - predictions) ^ 2)    
    }
  }
  # Average prediction error over the folds, and find optimal model length
  avg.error <- colMeans(cv.errors)
  tolerance <- min(avg.error)
  subset <- avg.error[avg.error <= tolerance]
  cv.optimal <- as.numeric(names(subset)[1])

  # Init list of best model formulas
  best.formulae <- vector("list", length = 0)
  # Find optimal length based on cros-validation
  GetFormula <- function(length) {
    best <- summary(fitted.all) [["outmat"]] [length, , drop = FALSE]
    best.predictors <- colnames(best)[best == "*"]
    # If any best predictors are factor levels, use the factor name instead.
    # Some serious magic going on here...
    indices <- which(best.predictors %in% paste0(names(factor.levels),
                                                 unlist(factor.levels)
                                                 )
                     )
    for (fac in seq_along(factor.levels)) {
      factor.levels[[fac]] <-
        paste0(names(factor.levels)[fac], factor.levels[[fac]])
    }
    for (i in indices) {
      best.predictors[i] <-
        llply(factor.levels, function(x) { best.predictors[i] %in% x } ) %>%
        { which(. == TRUE) } %>%
        names()
    }
    best.predictors %<>% unique()
    best.formula <- paste(response, "~", best.predictors[1])
    for (i in seq_len(length(best.predictors) - 1)) {
      best.formula <- paste(best.formula, "+", best.predictors[i + 1])
    }
    as.formula(best.formula)
  }
  best.formulae$cross.validation <- GetFormula(cv.optimal)
  # Find optimal length based on Rsq
  adj.rsq <- summary(fitted.all) [["adjr2"]]
  rsq.optimal <- which.max(adj.rsq)
  best.formulae$adj.rsq <- GetFormula(rsq.optimal)
  # Find optimal length based on Cp (AIC)
  cp <- summary(fitted.all) [["cp"]]
  cp.optimal <- which.min(cp)
  best.formulae$cp <- GetFormula(cp.optimal)
  # Find optimal length based on BIC
  bic <- summary(fitted.all)[["bic"]]
  bic.optimal <- which.min(bic)
  best.formulae$bic <- GetFormula(bic.optimal)

  # Compute the optimal models
  best.mods <- vector("list", length = 0)
  best.mods$cross.validation <- lm(best.formulae$cross.validation, data = x.complete)
  best.mods$adj.rsq <- lm(best.formulae$adj.rsq, data = x.complete)
  best.mods$cp <- lm(best.formulae$cp, data = x.complete)
  best.mods$bic <- lm(best.formulae$bic, data = x.complete)

  # Create matrix with model performance indicators
  mat <- matrix(data = NA,
                nrow = 4,
                ncol = 5,
                dimnames = list(c("cross.validation", "adj.rsq", "cp", "bic"),
                                c("mod.length", "Rsq", "AIC", "BIC", "max.VIF")
                                )
                )
  mat[, 1] <- c(cv.optimal, rsq.optimal, cp.optimal, bic.optimal)
  for (i in seq_len(nrow(mat))) {
    mat[i, 2] <- summary(best.mods[[i]]) [["r.squared"]]
    mat[i, 3] <- AIC(best.mods[[i]])
    mat[i, 4] <- BIC(best.mods[[i]])
    if (length(attr(summary(best.mods[[i]])$terms, "term.labels")) > 1) {
      mat[i, 5] <- max(vif(best.mods[[i]]))
    } else {
      mat[i, 5] <- 0
    }
  }

  # Return combined output
  list(formulae = best.formulae,
       mods = best.mods,
       performance = mat
       )
}

# Example code
if (FALSE) {
  fric <- read.csv(file = "output/FD_summary_FRic.csv", row.names = 1)
  env <- read.csv(file = "output/tdwg3_environmental_predictors.csv", row.names = 1)
  # These are all approximately normally distributed numerical predictors.
  x <- env
  x$FRic <- fric$FRic.global.SES
  response <- "FRic"
  k <- 10
  x <- x[, c(1:9, 12:24, 27:28, 29)]
  # Using bio7, instead of bio5 + bio6
  # using LGM climate anomaly, instead of LGM climate

  # Example fittings with different methods:
  a <- SelectOLS(x, response, k, method = "forward")
  b <- SelectOLS(x, response, k, method = "backward")
  c <- SelectOLS(x, response, k, method = "seqrep")
  d <- SelectOLS(x, response, k, method = "exhaustive")

  compare <- list(forward = a[[3]],
                  backward = b[[3]],
                  seqrep = c[[3]],
                  exhaustive = d[[3]]
                  )
  # Interpreting results remains non-trivial...
}



SingleOLS <- function(x, response, standardize = TRUE, digits = "all",
                      numeric.only = FALSE, double.std = FALSE) {
# Fit every single-predictor OLS model and report the Rsq, slope and P-value.
#
# Args:
#   x: data.frame with all variables (response + predictors)
#   response: character vector giving column name of the response variable.
#   standardize: Should predictors be standardized to mean 0 and unit variance?
#     NOTE: if this is FALSE, then the reported slopes are not directly comparable!
#   digits: integer giving number of decimal places to report, or "all" to skip
#     rounding.
#   numeric.only: Subset dataset to only continuous numerical predictors?
# Returns:
#  A matrix with 3 columns (Rsq, slope, p-value) reporting results for each single
#  predictor.


  # Validate data: subset to complete, numeric data and standardize
  x.complete <- ParseData(x = x,
                          response = response,
                          standardize = standardize,
                          numeric.only = numeric.only,
                          double.std = double.std
                          )

  # Init output matrix
  y.index <- match(response, colnames(x.complete))
  mat <- matrix(data = NA,
                ncol = 3,
                nrow = ncol(x.complete) - 1,
                dimnames = list(colnames(x.complete)[-y.index],
                                c("Rsq", "slope", "p-value")
                                )
                )

  # Fit each single-predictor model and extract results
  for (index in seq_len(ncol(x.complete))[-y.index]) {
    mod.summary <-
      lm(as.formula(paste(response, "~", colnames(x.complete)[index])),
         data = x.complete
         ) %>%
      summary()
    mat[match(names(x.complete)[index], rownames(mat)), 1] <-
      mod.summary[["r.squared"]]
    mat[match(names(x.complete)[index], rownames(mat)), 2] <-
      mod.summary[["coefficients"]] [2, 1]  # slope (of first predictor)
    # p-value: return the MODEL p-value. This is a bit more complex to obtain.
    f <- mod.summary[["fstatistic"]]
    mat[match(names(x.complete)[index], rownames(mat)), 3] <-
      pf(f[1], f[2], f[3], lower.tail = FALSE)  # p-value
  }

  # Return results, rounding as specified
  if (!identical (digits, "all")) {
    mat <- round(mat, digits = digits)
  }
  mat
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

  a <- SingleOLS(x, response)
  b <- SingleOLS(x, response, digits = 4)
}


# Function to generate single-predictor model data
SingleModels <- function(responses, predictors, double.std = FALSE) {
# responses: data frame with response variables
# Predictors: data frame with predictor variables

  # init output list
  mat <- matrix(data = NA,
                nrow = ncol(responses),
                ncol = ncol(predictors),
                dimnames = list(colnames(responses), colnames(predictors))
                )
  output <- list(mat, mat, mat)
  names(output) <- c("Rsq", "slope", "p-value")

  # Generate single-predictor models for each response variable
  for (i in seq_along(responses)) {
    response <- colnames(responses)[i]
    model.data <- predictors
    model.data[, response] <- responses[, response]
    mat <- SingleOLS(model.data,
                     response,
                     standardize = TRUE,
                     double.std = double.std
                     )
    indices <- match(rownames(mat), colnames(output[[1]]))
    output[[1]] [i, indices] <- mat[, 1]
    output[[2]] [i, indices] <- mat[, 2]
    output[[3]] [i, indices] <- mat[, 3]
  }
  output
}


# Wrapper function to automate the auto-selection process
# -------------------------------------------------------
MultiSelect <- function(response.var, response.name, predictors, k = 10,
                        double.std = FALSE) {
  model.data <- predictors
  model.data[, response.name] <- response.var
  all.mods <- vector("list", length = 4)
  names(all.mods) <- c("exhaustive", "backward", "forward", "seqrep")
  for (i in seq_along(all.mods)) {
    all.mods[[i]] <- SelectOLS(x = model.data,
                               response = response.name,
                               method = names(all.mods)[i],
                               k = k,
                               double.std = double.std
                               )
  }
  all.mods
}



FitGlobalOLS <- function(response, predictors, double.std = FALSE,
                         numeric.only = TRUE) {
# Function to generate a global OLS model object. Provided data is subsetted to
# complete cases via ParseData().
# NOTE: this function uses global assignments, for compatibility with dredge() from
# package MuMIn.
# Args:
#   response: vector of response data
#   predictors: dataframe with predictor variables, in the same order as the response
#               data.
  ols.mod.data <<-
    data.frame(predictors, response = response) %>%
    ParseData(.,
              response = "response",
              standardize = TRUE,
              numeric.only = numeric.only,
              double.std = double.std
              )
  ols.mod.formula <<-
    paste0("response ~ ", paste(colnames(predictors), collapse = " + ")) %>%
    as.formula()
  lm(formula = ols.mod.formula, data = ols.mod.data, na.action = na.fail)
}



PlotResids <- function(df, filename) {
# A Function to generate a multi-graph
# Assumes the first column in df is the response variable
  plot.rows <- ceiling((ncol(df) - 1) / 3)
  DoPlot <- function() {
    par(mfrow = c(plot.rows, 3))
    for (name in colnames(df)[-1]) {
      mod <- 
        paste(colnames(df)[1], "~", name) %>%
        as.formula() %>%
        lm(., data = df)
      Histogram(resid(mod), xlab = paste0("Residuals (", name, ")"))
    }
  }
  GraphSVG(DoPlot(),
           file = filename,
           width = 9,
           height = plot.rows * 3
           )
  return (0)
}



FDResidPlot <- function(index, null.model, predictors, filename) {
# A wrappper function to apply PlotResids()
# filename without extension!
  df <- data.frame(GetFD(fd.indices[index], null.model),
                   predictors
                   )
  PlotResids(df, filename = paste0(filename, "_full.svg"))
  df.realms <- RealmSubset(df)
  for (i in seq_along(df.realms)) {
    PlotResids(df.realms[[i]][, !colnames(df.realms[[i]]) %in% "realm"],
               filename = paste0(filename, "_", names(df.realms)[i], ".svg")
               )
  }
  return (0)
}

