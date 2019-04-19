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
# SingleSAR

library(magrittr)
library(plyr)
library(spdep)

source(file = "functions/base_functions.R")




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
      moran.test(resid(sar.mod), listw = nb.soi.swmat)[["p.value"]]
    # Pseudo.Rsq
    # No such thing as a real Rsq, because SAR is likelihood-based.
    # However, Kissling et al (2012a) use Pseudo-Rsq values, based on the method of
    # Kissling & Carl (2008), defined as "squared Pearson correlation between
    # predicted and observed values."
    mat[match(names(x.complete)[index], rownames(mat)), 4] <-
      cor(x = fitted(sar.mod), y = x.complete[, response]) ^ 2
    # BEWARE: this rsq includes the influence of the spatial autocorrelation term,
    # not only the effect of the tested predictor. If the spatial structure is
    # strong, pseudo-Rsq can be high even if the included predictor is not
    # significant.
    # To see this effect, you could include a 'null' predictor with random values
  }

  # Return results, rounding as specified
  if (!identical(digits, "all")) {
    mat <- round(mat, digits = digits)
  }
  mat
}

