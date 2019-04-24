###################################################
# Palm FD project: Preparation for SAR error models
###################################################
#
# In which automated model selection and  multimodel averaging is applied for SAR
# error models correlating palm functional diversity indices to environmental data.
#
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(sf)
library(ncf)
library(spdep)
library(MuMIn)
library(parallel)

source(file = "functions/base_functions.R")
source(file = "functions/SAR_regression_functions.R")
source(file = "functions/OLS_regression_functions.R")
source(file = "functions/plotting_functions.R")

# Graphics output directory (with trailing slash!)
plot.dir <- "graphs/SAR_multimod/"

# Single-predictor model output directory (with trailing slash!)
output.dir <- "output/SAR_multimod/"

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

set.seed(357)


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
fd.names <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
              "FRic.fruit.length", "FDis.all.traits", "FDis.stem.height",
              "FDis.blade.length", "FDis.fruit.length"
              )
fric <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
          "FRic.fruit.length")
fdis <- c("FDis.all.traits", "FDis.stem.height", "FDis.blade.length",
          "FDis.fruit.length")
for (index in fd.names) {
  fd.indices [[index]] <-
    read.csv(file = paste0("output/FD_summary_", index, ".csv"), row.names = 1)
}


cat("Preparing predictor variables...\n")
# ---------------------------------------
env <- read.csv(file = "output/tdwg3_predictors_complete.csv",
                 row.names = 1
                 )

# The following are the predictors of primary interest:
pred.names <- c("alt_range", "soilcount", "CH_SD", "bio1_sd", "bio12_sd",
                "bio4_mean", "bio15_mean", "lgm_Tano", "lgm_Pano"
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
# alt_range is the only variable that was excluded.
# The correlogram shows it is strongly correlated with bio1_sd, which makes sense.
# Restricting the VIF threshold to <2.5 additionaly removes bio1_sd.
# With bio1_sd also removed, the highest remaining VIF is <1.4, which is awesome.
# This should be considered given Cade (2015), depending on how strong th effect
# of bio1_sd is in the single-predictor models (i.e. is it likely significant?)
env.complete <- env[, pred.names]



################################
# Automated best model selection
################################
cat("Automated best SAR model selection:\n")

DredgeBestMod <- function(response, predictors, mod.formula, tdwg.map,
                          dist.weight = FALSE,
                          mc.cores = getOption("mc.cores", 2L)) {
# Construct a global SAR error model, run model selection with dredge(), and return
# the best model as fitted sarlm object.
#
# mod.formula should be a character string, NOT a formula object

  # Combine input data and prepare spatial weights matrix
  mod.data <- 
    data.frame(predictors, response = response) %>%
    ParseData(., response = "response", standardize = TRUE, numeric.only = TRUE)
  ind.complete <- complete.cases(predictors) & complete.cases(response)
  tdwg.map.subset <- tdwg.map[which(ind.complete), ]
  nb <- SoiNB(tdwg.map.subset)
  swmat <- SWMat(nb, tdwg.map.subset, dist.weight = dist.weight, style = "W")

  # Fit global model
  mod.formula <<- as.formula(mod.formula)
  # Yes that's high-level assignment. dredge() needs it.
  sar.mod.global <- errorsarlm(formula = mod.formula,
                               data = mod.data,
                               listw = swmat,
                               na.action = na.fail
                               )

  # Initiate cluster:
  cluster <- makeCluster(mc.cores)
  # Export environment to cluster
  clusterExport(cluster,
                varlist = c("mod.formula", "mod.data", "swmat"),
                envir = environment()
                )
  clusterEvalQ(cluster, library(spdep))
  # Fit all models
  cat("\tFitting all models...\n")
  fits <- pdredge(global.model = sar.mod.global,
                  cluster = cluster,
                  beta = "none",
                  rank = "AICc",
                  extra = "PseudoRsq"
                  )
  # Extract model coefficients
  cat("\tGet model coefficients...\n")
  fits.table <- model.sel(fits, beta = "none")
  fits.table %<>% as.data.frame()
  # model averaging
  cat("\tmodel averaging...\n")
  fits.avg <- model.avg(fits, beta = "none", cluster = cluster)
  cat("\tCalculating statistics of results...\n")
  fits.coefs <- coefTable(fits.avg, full = TRUE, adjust.se = FALSE)
  CI95 <- confint(fits.avg, level = 0.95, full = TRUE)
  stats <- cbind(fits.coefs, CI95)[, -3]

  # Gracefully end cluster
  stopCluster(cluster)

  # Return output as list
  list(global.mod = sar.mod.global,
       fits = fits,
       coefs = fits.table,
       mod.avg = fits.avg,
       coef.avg = stats
       )
}




# cases:
#   - global + 3 realm subsets (4)
#   - FRic 4 null models + FDis observed (5)
#   - all traits + 3 single traits (4)
#   - with/without distance-weighted swmat (2)
# 4 * 5 * 4 * 2 = 160 runs!

cat("Test run: FRic adf all traits in New World\n")
fdis.data <-
  GetFD(fd.indices[fric], "adf.SES") %>%
  RealmSubset()
env.complete.realms <- RealmSubset(env.complete)

test <- DredgeBestMod(response = fdis.data[["NewWorld"]] [, "FRic.all.traits"],
                      predictors = env.complete.realms[["NewWorld"]],
                      mod.formula = "response ~ soilcount + CH_SD + bio1_sd + bio12_sd + bio4_mean + bio15_mean + lgm_Tano + lgm_Pano",
                      tdwg.map = tdwg.map,
                      dist.weight = TRUE,
                      mc.cores = num.cores
                      )
# Using Histogram(), all coefficient distributions appear to be unimodal.
# Some better than others, but meh.



cat("Done.\n")


stop("enough")

# Multimodel averaging
# --------------------
# Adapted/expanded from https://drewtyre.rbind.io/post/rebutting_cade/ and also
# https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
#
#   (1) fit all possible models:
fits <- dredge(global.model, beta = "partial.sd")
# Standardization by partial sd is reccommended by Cade (2015).
# TODO: does this work in conjunction with (regularly) standardized predictors?
# TODO: 'partial.sd' does not work for sarlm models :(
#  (2) Obtain model selection table:
fits.table <- model.sel(fits, beta = "partial.sd")
fits.table %<>% as.data.frame()
#  (3) From this table, extract estimated coefficients and investigate their
#      normality. They should have UNIMODAL normal distributions (Cade, 2015).
# < nonstandard >
#  (4) Do model averaging:
fits.avg <- model.avg(fits, beta = "partial.sd")
fits.coefs <- coefTable(fits.avg, full = TRUE, adjust.se = TRUE)
#  (5) Get a 95% confidence interval for coefficients
fits.conf <- conf.int(fits)
# or directly from the 'fits.coefs' object as "estimate' +/- 1.96 * "Std.Error"
# >>> Compare these two: they should be equal (in which case, use conf.int)
#
# Naive variable importance is obtained with importance(fits), based on the
# Akaike weights.
# THIS IS VERY POOR INFORMATION. DON'T USE IT.
#
# You might consider extracting the best model and running moran.mc() on it
#
# We are NOT getting into variance partitioning. Not enough time, and I am not
# certain it is helpful: I think the standardized model-averaged coefficients
# are sufficiently informative, especially when paired with single-predictor
# model Rsq.
#
#  (6) Make a plot (only for hand-picked cases)
#
#
#
#  (7) For perspective, you could report the pseudo-Rsq of the full model,
#      and possibly even of the single 'best' model, although that can be
#      questionable if you don't put it in the right context.
