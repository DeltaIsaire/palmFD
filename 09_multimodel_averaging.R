###########################################################
# Palm FD project: Multimodel Averaging on OLS + SAR models
###########################################################
#
# In which automated model selection and multimodel averaging is applied for OLS
# and SAR-error models correlating palm functional diversity indices to
# environmental variation.


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(sf)
library(ncf)
library(spdep)
library(MuMIn)
library(parallel)
library(spatialreg)
library(car)

source(file = "functions/base_functions.R")
source(file = "functions/SAR_regression_functions.R")
source(file = "functions/plotting_functions.R")

# Graphics output directory (with trailing slash!)
plot.dir <- "graphs/multimodels/"

# Model data output directory (with trailing slash!)
output.dir <- "output/multimodel_averaging/"

# Number of cores to use for parallel processing. Default is 95% of available cores.
num.cores <- 
  if (!is.na(detectCores())) {
    floor(detectCores() * 0.95)
  } else {
    1
    warning("unable to detect cores. Parallel processing is NOT used!")
    cat("unable to detect cores. Parallel processing is NOT used!")
  }


if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
}
if (!dir.exists(output.dir)) {
  cat("creating directory:", output.dir, "\n")
  dir.create(output.dir)
}



############################################
# Load functional diversity and spatial data
############################################
cat("loading data...\n")

# TDWG3 data
# ----------
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")

tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")

# List with the tdwg3 units in each realm
realm.tdwg3 <- vector("list", length = 3)
names(realm.tdwg3) <- levels(tdwg3.info[, "realm"])
for (i in 1:3) {
  realm.tdwg3[[i]] <-
    tdwg3.info[tdwg3.info[, "realm"] == names(realm.tdwg3)[i], "tdwg3.code"]
}


# Functional diversity indices
# ----------------------------
# A list with all the FD index data.
# Note these dataframes include richness and endemism variables.
fd.indices <- vector("list", length = 0)
fric <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
          "FRic.fruit.length")
fdis <- c("FDis.all.traits", "FDis.stem.height", "FDis.blade.length",
          "FDis.fruit.length")
fd.names <- c(fric, fdis)
for (index in fd.names) {
  fd.indices [[index]] <-
    read.csv(file = paste0("output/FD_summary_", index, ".csv"), row.names = 1)
}


cat("Preparing predictor variables...\n")
# ---------------------------------------
env <- read.csv(file = "output/tdwg3_predictors_complete_noch.csv",
                row.names = 1
                )

# The following are the predictors of primary interest:
# UPDATE: we no longer include CH_SD or bio1_mean or bio12_mean
pred.names <- c("alt_range", "soilcount", "bio1_sd", "bio12_sd", "bio4_mean",
                "bio15_mean", "lgm_Tano", "lgm_Pano", "plio_Tano", "plio_Pano",
                "mio_Tano", "mio_Pano"
                )
# Collinearity must be accounted for, so filter the predictors with AutoVIF()
GraphSVG(Correlogram(env[, pred.names]),
         file = paste0(plot.dir, "Correlogram_of_predictors.svg"),
         width = 7,
         height = 7
         )
test.data <- env[, pred.names]
test.data[, "response"] <- fd.indices[["FRic.all.traits"]] [, "global.SES"]
pred.names <- AutoVIF(test.data,
                      response = "response",
                      standardize = TRUE,
                      threshold = 2.5,
                      numeric.only = TRUE
                      )
# alt_range was excluded at threshold 3.
# The correlogram shows it is strongly correlated with bio1_sd, which makes sense.

env.complete <- env[, pred.names]



################################
# Automated best model selection
################################
cat("Automated best model selection:\n")

# Function to fit global OLS model
# --------------------------------
FitGlobalOLS <- function(response, predictors, double.std = FALSE,
                         numeric.only = TRUE) {
# Function to generate a global OLS model object. Provided data is subsetted to
# complete cases via ParseData().
# NOTE: this function uses global assignments, for compatibility with dredge() from
# package MuMIn.
# Args:
#   response: vector of response data
#   predictors: dataframe with predictor variables, in the same order as the
#   response data.
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


# Function for model averaging
# ----------------------------
DredgeBestMod <- function(global.model, beta = "none", cluster,
                          long.output = FALSE, ...) {
# Function to apply model averaging based on a global model object.
# The global model object should be called with 'na.action = na.fail', and all
# objects referenced in the model call should exist in the global environment. This
# is a condition required by dredge().
#
# Args:
#   global.model: the fitted global model object
#   beta: type of coefficient standardization to use, see dredge()
#   cluster: cluster to use for parallel computation in pdredge() and
#             model.avg()
#   long.output: return all model average objects in a list, instead of only
#                the coefficient average dataframe. The full objects are rather
#                heavy on RAM, so use long output only for testing and on a machine
#                with a lot of RAM
#   ...: additional arguments to pass to dredge()

  # Fit all models
  cat("\t\tFitting all models...\n")
#  fits <- pdredge(global.model = global.model,
#                  cluster = cluster,
#                  beta = beta,
#                  ...
#                  )
# PARALLEL FUNCTION DISABLED BECAUSE IT BREAKS WITH spatialreg::errorsarlm()
# IN R 3.6.0
  fits <- dredge(global.model = global.model, beta = beta, ...)
  # Extract model coefficients
  cat("\t\tGet model coefficients...\n")
  fits.table <-
     model.sel(fits, beta = beta) %>%
     as.data.frame()
  # model averaging
  cat("\t\tmodel averaging...\n")
#  fits.avg <- model.avg(fits, beta = beta, cluster = cluster)
  fits.avg <- model.avg(fits, beta = beta)
  cat("\t\tCalculating statistics of results...\n")
  fits.coefs <- coefTable(fits.avg, full = TRUE, adjust.se = FALSE)
  CI95 <- confint(fits.avg, level = 0.95, full = TRUE)
  stats <-
    cbind(fits.coefs, CI95)[, -3] %>%
    as.data.frame()
  best.coefs <- unlist(fits.table[1, ])
  stats[, "best.model"] <- best.coefs[match(rownames(stats), names(best.coefs))]

  # Return output as list
  if (long.output) {
    return (list(global.model = global.model,
                 fits = fits,
                 coefs = fits.table,
                 mod.avg = fits.avg,
                 coef.avg = stats
                 )
            )
  } else {
    return (list(coef.avg = stats))
  }
}


# Apply model averaging for each FD response variable
# ---------------------------------------------------
# prepare globals
null.models <- c("global.SES", "realm.SES", "realm.SES.noMDG", "adf.SES",
                 "observed")
env.complete.realms <- RealmSubset(env.complete)
n <- 1

# Prepare cluster for parallel processing
cluster <- makeCluster(num.cores)
clusterEvalQ(cluster, library(spdep))

for (index in fd.names) {
  for (null.model in null.models) {
    # Prepare data for cases full + realm subsets
    index.data <- GetFD(fd.indices[index], null.model)
    index.data.realms <- RealmSubset(index.data)
    response.list <- c(list(full = index.data[, index]), index.data.realms)
    predictors.list <- c(list(full = env.complete), env.complete.realms)
    cases <- names(response.list)

    for (case in cases) {
      # verbose output
      nmax <- length(fd.names) * length(null.models) * length(cases)
      cat("\t",
          paste0("(", n, " of ", nmax, ")"),
          paste(index, null.model, "for", case, "dataset"),
          "\n"
          )

      # Handle process flow: skip if output already exists
      header <-
        paste0(output.dir, "multimod_avg_", index, "_", null.model, "_", case, "_")
      if (file.exists(paste0(header, "SAR.csv"))) {
        cat("\tSKIP - output files already exist\n")
        n <- n + 1
        next
      }

      # Fit global models
      response <- response.list[[case]]
      predictors <- predictors.list[[case]]
      global.ols <- FitGlobalOLS(response = response,
                                 predictors = predictors,
                                 double.std = TRUE
                                 )
      global.sar <- suppressMessages(FitGlobalSAR(response = response,
                                                  predictors = predictors,
                                                  tdwg.map = tdwg.map,
                                                  dist.weight = TRUE,
                                                  double.std = TRUE
                                                  )
                                     )

      cat("Calculating global model explained variance and testing residual SAC..\n")
      # Calculate explained variance of global models and residual SAC
      # First init dataframe
      rsq.file <- paste0(output.dir, "00_multimod_avg_global_rsq.csv")
      if (file.exists(rsq.file)) {
        global.rsq <- read.csv(file = rsq.file, as.is = c(1:3))
      }
      if (!exists("global.rsq")) {
        global.rsq <-
          matrix(data = NA, nrow = nmax, ncol = 10,
                 dimnames = list(NULL,
                                 c("index", "null.model", "case", "OLS.Rsq",
                                   "SAR.PsRsq.Nagel", "SAR.PsRsq.KC",
                                   "SAR.PsRsq.pred", "OLS.resid.SAC.P",
                                   "SAR.resid.SAC.P", "SAR.bptest.P"
                                   )
                                 )
                 ) %>%
         as.data.frame()
      }
      # Second, calculate the statistics. 
      global.rsq[n, "index"] <- index
      global.rsq[n, "null.model"] <- null.model
      global.rsq[n, "case"] <- case
      # For SAR models we use 3 different methods of calculating pseudo-Rsq.
      global.rsq[n, "OLS.Rsq"] <- summary(global.ols)[["r.squared"]]
      # SAR 1: Nagelkerke pseudo-Rsq (Nagelkerke 1991, see ?summary.sarlm
      global.rsq[n, "SAR.PsRsq.Nagel"] <-
        summary(global.sar, Nagelkerke = TRUE)[["NK"]]
      # SAR 2: PseudoRsq following Kissling & Carl (2008)
      global.rsq[n, "SAR.PsRsq.KC"] <- cor(x = predict(global.sar),
                                           y = global.sar[["y"]],
                                           method = "pearson"
                                           ) ^ 2
      # SAR 3: PseudoRsq following Kissling & Carl (2008) but excluding lambda
      global.rsq[n, "SAR.PsRsq.pred"] <- PseudoRsq(global.sar)
      # Test for spatial autocorrelation (SAC) in OLS model residuals
      global.rsq[n, "OLS.resid.SAC.P"] <-
        moran.test(resid(global.ols), listw = sar.swmat)[["p.value"]]
      # Test for SAC in SAR model residuals
      global.rsq[n, "SAR.resid.SAC.P"] <-
        moran.test(resid(global.sar), listw = sar.swmat)[["p.value"]]
      # bruch pagan test for heteroskedasticity in SAR model
      global.rsq[n, "SAR.bptest.P"] <- bptest.sarlm(global.sar)[["p.value"]]
      # Third, save (partial) dataframe to disk
      write.csv(global.rsq, file = rsq.file, row.names = FALSE)

      cat("Plotting partial residuals of SAR model...\n")
      # Assemble spatially-corrected data
      sar.y <- global.sar[["tary"]]
      sar.x <- as.data.frame(global.sar[["tarX"]])
      names(sar.x) <- names(global.sar[["coefficients"]])
      sar.data <- cbind(sar.y, sar.x[, -1])
      # Fit linear model to this data. The result is not the same as a real SAR
      # model, but the coefficient estimates and therefore the residuals are
      # (approximately) the same, which we can then plot
      tarlm <- lm(sar.y ~ ., data = sar.data)
      # Make and save crPlot
      GraphPNG(crPlots(tarlm, layout = c(4, 3), main = paste("Partial Residuals")),
               file = paste0(plot.dir,
                             "partial_resids_",
                             index, "_",
                             null.model, "_",
                             case, ".png"
                             ),
               width = 900,
               height = 1200,
               pointsize = 18
               )

      # Export GLOBAL environment to cluster
      clusterExport(cluster,
                    varlist = ls(name = .GlobalEnv),
                    envir = .GlobalEnv
                    )

      # Perform multimodel averaging
      # For Realm subsets, to prevent overfitting we limit the models to having
      # at most 5 predictors.
      if (!case == "full") {
        m.max <- 5
        cat("Limiting model length to", m.max, "predictors\n")
      } else {
        m.max <- NA
      }

      cat("\tOLS:\n")
      mod.avg.ols <- DredgeBestMod(global.ols,
                                   beta = "none",
                                   cluster = cluster,
                                   m.max = m.max
                                   )
      cat("\tOLS with standardization:\n")
      mod.avg.ols.std <- DredgeBestMod(global.ols,
                                       beta = "partial.sd",
                                       cluster = cluster,
                                       m.max = m.max
                                       )
      cat("\tSAR error:\n")
      mod.avg.sar <- DredgeBestMod(global.sar,
                                   beta = "none",
                                   cluster = cluster,
                                   m.max = m.max
                                   )

      # Save results summary to disk
      write.csv(mod.avg.ols[["coef.avg"]],
                file = paste0(header, "OLS.csv"),
                eol = "\r\n"
                )
      write.csv(mod.avg.ols.std[["coef.avg"]],
                file = paste0(header, "OLS_std.csv"),
                eol = "\r\n"
                )
      write.csv(mod.avg.sar[["coef.avg"]],
                file = paste0(header, "SAR.csv"),
                eol = "\r\n"
                )

      n <- n + 1
    }
  }
}


# Gracefully end cluster
stopCluster(cluster)



##############################################
cat("Processing model averaging results...\n")
##############################################

# Function to load output files into array
# ----------------------------------------
ParseDredge <- function(header, fd.names, null.models, cases, model) {
  # Init array of results, getting basic dimensions from a template
  template <- read.csv(file = paste0(header,
                                     fd.names[1], "_",
                                     null.models[1], "_",
                                     cases[1], "_",
                                     model, ".csv"
                                     ),
                       row.names = 1
                       )
  output <- array(dim = c(nrow(template),
                          ncol(template),
                          length(fd.names),
                          length(cases),
                          length(null.models)
                          ),
                  dimnames = list(predictors = rownames(template),
                                  statistics = colnames(template),
                                  fd.indices = fd.names,
                                  cases      = cases,
                                  null.model = null.models
                                  )
                  )

  # Load data from files into array, matching row and column names
  for (null.model in null.models) {
    for (case in cases) {
      for (index in fd.names) {
        # load data file
        data <- read.csv(file = paste0(header,
                                       index, "_",
                                       null.model, "_",
                                       case, "_",
                                       model, ".csv"
                                       ),
                         row.names = 1
                         )
        # Match row and column names
        data <- data[match(rownames(output), rownames(data)), ]
        data <- data[, match(colnames(output), colnames(data))]
        # Insert data into array
        output[, , index, case, null.model] <- as.matrix(data)
      }
    }
  }

  output
}


# Apply ParseDredge() to load data
# --------------------------------
avg.FDis.OLS <- ParseDredge(header = paste0(output.dir, "multimod_avg_"),
                            fd.names = fdis,
                            null.models = null.models,
                            cases = cases,
                            model = "OLS"
                            )

avg.FDis.OLS.std <- ParseDredge(header = paste0(output.dir, "multimod_avg_"),
                                fd.names = fdis,
                                null.models = null.models,
                                cases = cases,
                                model = "OLS_std"
                                )

avg.FDis.SAR <- ParseDredge(header = paste0(output.dir, "multimod_avg_"),
                            fd.names = fdis,
                            null.models = null.models,
                            cases = cases,
                            model = "SAR"
                            )

avg.FRic.OLS <- ParseDredge(header = paste0(output.dir, "multimod_avg_"),
                            fd.names = fric,
                            null.models = null.models,
                            cases = cases,
                            model = "OLS"
                            )

avg.FRic.OLS.std <- ParseDredge(header = paste0(output.dir, "multimod_avg_"),
                                fd.names = fric,
                                null.models = null.models,
                                cases = cases,
                                model = "OLS_std"
                                )

avg.FRic.SAR <- ParseDredge(header = paste0(output.dir, "multimod_avg_"),
                            fd.names = fric,
                            null.models = null.models,
                            cases = cases,
                            model = "SAR"
                            )


# (Manual) comparison of results
# ------------------------------
# Array tables make the data accessible, with the incantation
#   object[, , index, case, null.model]
# e.g.
avg.FDis.SAR[, , "FDis.all.traits", "full", "observed"]
#
# Fun stuff can be done with apply functions. For instance, find if the 95% CI
# includes 0 (i.e. the average coefficient is not significantly different from 0):
CheckCI <- function(x, value = 0) {
  output <- value >= x[, "X2.5.."] & value <= x[, "X97.5.."]
  sign.ind <- which(output == FALSE)
  output[sign.ind] <- x[sign.ind, "Estimate"]
  output[output == TRUE] <- "n.s."
  output
}
adply(avg.FDis.SAR, .margins = c(3, 4, 5), CheckCI)

# Or, we could use a function to report the estimate +/- 95% CI as a string
ReportCI <- function(x, digits = 3) {
  x <- round(x, digits = digits)
  output <- paste0(x[, "Estimate"], " (", x[, "X2.5.."], " - ", x[, "X97.5.."], ")")
  names(output) <- rownames(x)
  output
  }


write.csv(adply(avg.FDis.OLS, .margins = c(3, 4, 5), ReportCI),
          file = paste0(output.dir, "00_95CI_FDis_OLS.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )

write.csv(adply(avg.FDis.OLS.std, .margins = c(3, 4, 5), ReportCI),
          file = paste0(output.dir, "00_95CI_FDis_OLS_std.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )

write.csv(adply(avg.FDis.SAR, .margins = c(3, 4, 5), ReportCI),
          file = paste0(output.dir, "00_95CI_FDis_SAR.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )

write.csv(adply(avg.FRic.OLS, .margins = c(3, 4, 5), ReportCI),
          file = paste0(output.dir, "00_95CI_FRic_OLS.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )

write.csv(adply(avg.FRic.OLS.std, .margins = c(3, 4, 5), ReportCI),
          file = paste0(output.dir, "00_95CI_FRic_OLS_std.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )

write.csv(adply(avg.FRic.SAR, .margins = c(3, 4, 5), ReportCI),
          file = paste0(output.dir, "00_95CI_FRic_SAR.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )

cat("Done.\n")

