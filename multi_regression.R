library(corrplot)
library(leaps)
library(car)

##############################
# Correlation plots of FD data
##############################
# exploratory correlations, to see how single-trait FD relates to all-trait FD,
# and how FRic relates to FDis.
# Here we focus on the stochastic data.

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

# Observed FD
# -----------
Correlogram(fd.observed.stochastic)
# FDis and FRic are not strongly correlated, except for blade length and
#   stem height.
# all.traits FRic is moderately correlated with each single trait FRic.
# all.traits FDis is correlated with FDis blade length and fruit length,
#   but not stem height.

# We can explore the predictors of FRic and FDis more closely.

mod.FRic <- lm(FRic.all.traits ~ FRic.stem.height + FRic.blade.length +
               FRic.fruit.length, data = fd.observed.stochastic
               )
# RSQ 0.685, all VIF < 2. Removing any predictors does not decrease AIC.
# Residuals are quite normally distributed. This model is optimal.
# All predictors highly significant of course, all slopes positive

mod.FDis <- lm(FDis.all.traits ~ FDis.stem.height + FDis.blade.length +
               FDis.fruit.length, data = fd.observed.stochastic
               )
# RSQ 0.971, all VIF < 2. Removing any predictors does not decrease AIC.
# Residuals are a bit skewed. This model is optimal.
# All predictors highly significant, all slopes positive.



# Cross-validated model selection using forward and backward selection
# based on https://uc-r.github.io/model_selection


# Forward selection
# -----------------
mod.forward <- 
  regsubsets(FRic.all.traits ~ .,
             data = fd.observed.stochastic[, -1],
             nvmax = ncol(fd.observed.stochastic[, -1]),
             method = "forward"
             )
forward <- summary(mod.forward)
# tells you which variables are included in the best model of a given size.
# Subsequently, you need to know which size model is best.
# The summary provides three statistics: adj. rsq, BIC and Cp (kind of like AIC).
select <- data.frame(adj.rsq = forward$adjr2,
                     BIC     = forward$bic,
                     cp      = forward$cp,
                     row.names = paste0("mod.size.", seq_len(length(forward$cp)))
                     )
best.size.forward <- 
  c(which.max(select$adj.rsq), which.min(select$BIC), which.min(select$cp))
names(best.size.forward) <- names(select)

# Backward selection
# ------------------
mod.backward <- 
  regsubsets(FRic.all.traits ~ .,
             data = fd.observed.stochastic[, -1],
             nvmax = ncol(fd.observed.stochastic[, -1]),
             method = "backward"
             )
backward <- summary(mod.backward)
# tells you which variables are included in the best model of a given size.
# Subsequently, you need to know which size model is best.
# The summary provides three statistics: adj. rsq, BIC and Cp (kind of like AIC).
select <- data.frame(adj.rsq = backward$adjr2,
                     BIC     = backward$bic,
                     cp      = backward$cp,
                     row.names = paste0("mod.size.", seq_len(length(backward$cp)))
                     )
best.size.backward <- 
  c(which.max(select$adj.rsq), which.min(select$BIC), which.min(select$cp))
names(best.size.backward) <- names(select)

# Compare best models using cross-validation
# ------------------------------------------
CrossValidate <- function(x, response, k = 10) {
  x.complete <- x[complete.cases(x), ]
  folds <- sample(sample.int(k), nrow(x.complete), replace = TRUE)  ## FALSE!
  cv.errors <-
    matrix(data = NA,
           nrow = k,
           ncol = ncol(x.complete) - 1,
           dimnames = list(NULL, seq_len(ncol(x.complete) - 1))
           )
  # Forward selection for each of the k training datasets
  for (i in seq_len(k)) {
    forward <- regsubsets(as.formula(paste(response, "~", ".")),
                          data = x.complete[!folds == i, ],
                          nvmax = ncol(x.complete) - 1,
                          method = "forward"
                          )
    # Predict test values using the best model of each size,
    # and find the mean square difference with the test data
    for ( j in seq_len(ncol(x.complete) - 1)) {
      my.formula <- as.formula(paste(response, "~", "."))
      mat <- model.matrix(my.formula, x.complete[folds == i, ])
      coeff <- coef(forward, id = j)
      xvars <- names(coeff)
      predictions <- mat[, xvars] %*% coeff
      cv.errors[i, j] <-
        mean((x.complete[, response][folds == i] - predictions) ^ 2)    
    }
  }
  cv.errors
}

x <- fd.observed.stochastic[, -1]
response <- "FRic.all.traits"

CrossValidate(x, response)




# Global null model SES
# ---------------------
Correlogram(fd.global.stochastic)
# FDis and FRic are moderately correlated.
# FDis is correlated more strongly with FDis blade length and fruit length than 
#   stem height.
# FRic is only modestly correlated with single-trait FRic indices.

# Realm null model SES
# --------------------
Correlogram(fd.regional.stochastic)
# More or less the same thing

# dispersion field null model SES
# -------------------------------
Correlogram(fd.local.stochastic)
# Again, qualitatively the same thing.



# Overall, it looks like z-scores for FRic and FDis are quite correlated,
# and this is driven mostly by fruit length and blade length, and not stem height.





