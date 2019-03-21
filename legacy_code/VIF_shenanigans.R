################################################
# Filtering predictors by collinearity using VIF
################################################
cat("Filtering predictors by collinearity using VIF:\n")

# Function to do the hard work
# ----------------------------
FilterPred <- function(fd, predictors, name) {
# Name: character string as unique data identifier. Used as suffix for file names.

  # init output list
  output <- vector("list", length = ncol(fd))
  names(output) <- colnames(fd)
  # Select non-collinear predictors for each FD metric
  for (i in seq_along(fd)) {
    response <- colnames(fd)[i]
    model.data <- predictors
    model.data[, response] <- fd[, response]
    output[[i]] <- AutoVIF(model.data, response, standardize = TRUE, threshold = 3)
  }
  # parse and save output
  if (length(unique(output)) == 1) {
    cat("Selected predictors are identical for each response variable.\n")
    write.csv(output[[1]],
              file = paste0(output.dir,
                            "OLS_noncollinear_predictors_", 
                            name,
                            ".csv"
                            ),
              eol = "\r\n"
              )
  } else {
    cat("Selected predictors differ between response variables.",
        "Output is NOT saved. We apologize for the inconvenience :(\n"
        )
  }
  output
}

# Apply the function
# ------------------
# The input data is already generated above, but we have to exclude predictors
# with linear dependencies. Hence,
# we use bio7 instead of bio5 + bio6,
# and LGM climate anomalies instead of absolute LGM climate.
exclude <- c("bio5_mean", "bio6_mean", "lgm_ens_Tmean", "lgm_ens_Pmean")

cat("(1) for null model Global...\n")
env.global %<>% .[, -which(names(.) %in% exclude)]
env.global.noCH %<>% .[, -which(names(.) %in% exclude)]

a <- FilterPred(fd.global, env.global, "global")
FilterPred(fd.global, env.global.noCH, "global_noCH")

cat("(2) for null model Realm...\n")
env.realm %<>% .[, -which(names(.) %in% exclude)]
env.realm.noCH %<>% .[, -which(names(.) %in% exclude)]

a <- FilterPred(fd.realm, env.realm, "realm")
FilterPred(fd.realm, env.realm.noCH, "realm_noCH")

cat("(2) for null model ADF...\n")
env.adf %<>% .[, -which(names(.) %in% exclude)]
env.adf.noCH %<>% .[, -which(names(.) %in% exclude)]

a <- FilterPred(fd.adf, env.adf, "adf")
FilterPred(fd.adf, env.adf.noCH, "adf_noCH")

