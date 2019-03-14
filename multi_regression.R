library(magrittr)
library(plyr)
library(corrplot)
library(leaps)
library(car)

# Load basic data
# ---------------
fric <- read.csv(file = "output/FD_summary_FRic.csv", row.names = 1)


env <- read.csv(file = "output/tdwg3_environmental_predictors.csv", row.names = 1)
# These are all approximately normally distributed numerical predictors.


###########################################
# Functions for (automated) model selection
###########################################

# Function to create corrgrams
# ----------------------------
Correlogram <- function(x) {
  subset <- x[, unlist(lapply(x, is.numeric))]
  corr.matrix <- 
    cor(subset[complete.cases(subset), ]) %>%
    round(., digits = 2)
  corrplot(corr.matrix, method = "square", order = "FPC", addCoef.col = "red")
  return (0)
}


# Compare best models using cross-validation
# ------------------------------------------
SelectOLS <- function(x, response, k = 10, standardize = TRUE, method = "forward") {
# Function to perform automated OLS regression model selection for a given dataset.
# The provided data is subsetted to complete cases of all numeric variables. The best
# OLS model of each length is found via 'regsubsets' from the 'leaps' package.
# Subsequently, the optimal length is estimated using cross-validation, adjusted Rsq,
# AIC and BIC.
#
# NOTE: beware of linear dependencies. With bioclim data for instance, bio7 equals
# bio5 - bio6, which doesn't provide any new information. You must use either 5 + 6
# OR 7, but not all three at once.
#
# based on https://uc-r.github.io/model_selection
#
# Args:
#   x: data.frame with all variables (response + predictors)
#   response: character vector giving column name of the response variable.
#   k: number of folds to use for k-fold cross-validation.
#   standardize: Should predictors be standardized to mean 0 and unit variance?
#   method: model selection method to use, see regsubsets()
# Returns:
#   dd

  # Validate data: subset to complete, numeric data and standardize
  x.num <- numcolwise(function(x) { x } ) (x)
  x.complete <- x.num[complete.cases(x.num), ]
  if (standardize) {
    x.std <- 
      apply(x.complete, 2, scale, center = TRUE, scale = TRUE) %>%
      as.data.frame()
    colnames(x.std) <- colnames(x.complete)
    x.std[, response] <- x.complete[, response]
    x.complete <- x.std
  }

  # Find best models of each length
  fitted.all <- regsubsets(as.formula(paste(response, "~", ".")),
                            data = x.complete,
                            nvmax = ncol(x.complete),
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
           dimnames = list(fold = seq_len(k), mod.length = seq_len(ncol(x.complete) - 1))
           )
  # Forward selection for each of the k training datasets
  for (i in seq_len(k)) {
    fitted <- regsubsets(as.formula(paste(response, "~", ".")),
                          data = x.complete[!folds == i, ],
                          nvmax = ncol(x.complete) - 1,
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
  # Average prediction error over the folds, and find optimal model length, defined
  # as the shortest model
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
    if (mat[i, 1] > 1) {
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













