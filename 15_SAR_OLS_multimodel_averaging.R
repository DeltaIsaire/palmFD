###########################################################
# Palm FD project: Multimodel Averaging on OLS + SAR models
###########################################################
#
# In which automated model selection and  multimodel averaging is applied for OLS
# and SAR-error models correlating palm functional diversity indices to
# environmental variation.
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
plot.dir <- "graphs/multimodels/"

# Single-predictor model output directory (with trailing slash!)
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
# alt_range is the only variable that was excluded at threshold 2.5.
# The correlogram shows it is strongly correlated with bio1_sd, which makes sense.
# Restricting the VIF threshold to <2.5 additionaly removes bio1_sd.
# With bio1_sd also removed, the highest remaining VIF is <1.4, which is awesome.
# Removing bio1_sd as well should be considered given Cade (2015), depending on how
# strong th effect of bio1_sd is in the single-predictor models (i.e. is it likely
# significant?)
env.complete <- env[, pred.names]



################################
# Automated best model selection
################################
cat("Automated best model selection:\n")


DredgeBestMod <- function(global.model, beta = "none", cluster, ...) {
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
#   ...: additional arguments to pass to dredge()

  # Fit all models
  cat("\tFitting all models...\n")
  fits <- pdredge(global.model = global.model,
                  cluster = cluster,
                  beta = beta,
                  ...
                  )
  # Extract model coefficients
  cat("\tGet model coefficients...\n")
  fits.table <-
     model.sel(fits, beta = beta) %>%
     as.data.frame()
  # model averaging
  cat("\tmodel averaging...\n")
  fits.avg <- model.avg(fits, beta = beta, cluster = cluster)
  cat("\tCalculating statistics of results...\n")
  fits.coefs <- coefTable(fits.avg, full = TRUE, adjust.se = FALSE)
  CI95 <- confint(fits.avg, level = 0.95, full = TRUE)
  stats <-
    cbind(fits.coefs, CI95)[, -3] %>%
    as.data.frame()
  best.coefs <- unlist(fits.table[1, ])
  stats[, "best.model"] <- best.coefs[match(rownames(stats), names(best.coefs))]

  # Return output as list
  list(global.model = global.model,
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

cat("Test run: FRic global all traits in New World\n")
fdis.data <-
  GetFD(fd.indices[fric], "global.SES") %>%
  RealmSubset()
env.complete.realms <- RealmSubset(env.complete)

response <- fdis.data[["NewWorld"]] [, "FRic.all.traits"]
predictors <- env.complete.realms[["NewWorld"]]


stop("test stop")

# Fit global models
# -----------------
global.ols <- FitGlobalOLS(response = response, predictors = predictors)
global.sar <- FitGlobalSAR(response = response,
                           predictors = predictors,
                           tdwg.map = tdwg.map,
                           dist.weight = TRUE
                           )

# Create cluster for parallel processing
# --------------------------------------
# Initiate cluster:
cluster <- makeCluster(num.cores)
# Export environment to cluster
clusterExport(cluster,
              varlist = ls(name = .GlobalEnv),
              envir = .GlobalEnv
              )
clusterEvalQ(cluster, library(spdep))


# Perform multimodel averaging
# ----------------------------
mod.avg.ols <- DredgeBestMod(global.ols, beta = "none", cluster = cluster)
mod.avg.ols[["coef.avg"]]
mod.avg.ols.std <- DredgeBestMod(global.ols, beta = "partial.sd", cluster = cluster)
mod.avg.ols.std[["coef.avg"]]
mod.avg.sar <- DredgeBestMod(global.sar, beta = "none", cluster = cluster)
mod.avg.sar[["coef.avg"]]


# Gracefully end cluster
# ----------------------
stopCluster(cluster)




cat("Done.\n")

