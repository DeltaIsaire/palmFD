##################################################################
# Palm FD project: Relation between trait distributions and realms
##################################################################
#
# In which the statistical relationship between community average trait values
# and realm is explored.
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)
library(spdep)
library(spatialreg)
library(sf)

source(file = "functions/base_functions.R")
source(file = "functions/OLS_regression_functions.R")
source(file = "functions/SAR_regression_functions.R")



#######################
# Load and prepare data
#######################
cat("preparing data...\n")

# tdwg3 info:
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")

tdwg.map <- read_sf(dsn = "/home/delta/R_Projects/palm_FD/data/tdwg",
                    layer = "TDWG_level3_Coordinates")

# traits, based on stochastic gapfilled trait matrix:
cwm.traits <- read.csv(file = "output/observed_FD/community_trait_means_stochastic_mean_of.csv",
                                    row.names = 1
                                    )

# all traits FD indices:
fdis <- read.csv(file = "output/FD_summary_FDis.all.traits.csv", row.names = 1)
fric <- read.csv(file = "output/FD_summary_FRic.all.traits.csv", row.names = 1)

# Palm occurrence data
palm.dist <- read.csv(file = "data/palms_in_tdwg3.csv")
pres.abs <-
  read.csv(file = "output/palm_tdwg3_pres_abs_gapfilled.csv",
           row.names = 1
           ) %>%
   as.matrix()


# List of TDWG3 units with complete env data, for subsetting.
env <- read.csv(file = "output/tdwg3_predictors_complete_noch.csv",
                row.names = 1
                )
tdwg3.complete <- rownames(env) [complete.cases(env)]


# Merge trait data with realm data
# --------------------------------
cwm.traits[, "realm"] <-
  tdwg3.info[match(rownames(cwm.traits), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()
# Subset
cwm.traits %<>% .[rownames(.) %in% tdwg3.complete, ]


# Merge FD data with realm data and subset
# ----------------------------------------
fdis[, "realm"] <-
  tdwg3.info[match(rownames(fdis), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()
fdis %<>% .[rownames(.) %in% tdwg3.complete, ]

fric[, "realm"] <-
  tdwg3.info[match(rownames(fric), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()
fric %<>% .[rownames(.) %in% tdwg3.complete, ]



########################################
# assemblage mean traits and FD vs realm
########################################
# Using SAR error models,
# because ANOVA does not account for spatial autocorrelation.
cat("Modelling trait and FD differences between realms...\n")

# For a full analysis we will construct a dataframe with the following columns:
# - SAR model intercept
# - SAR model slope
# - SAR model Pseudo-Rsquared
# - SAR model pseudo-R-squared of realm
# - SAR model significance of realm
# - SAR model lambda
# - SAR model significance of lambda
# - SAR model SAC in residuals
# - SAR model pbtest of homoskedasticity

# Gather realm data and response variables
# ----------------------------------------
div.data <-
  data.frame(realm            = cwm.traits$realm,
             cwm.traits[, c("stem.height", "blade.length", "fruit.length")],
             fric.observed    = fric$observed,
             fric.global      = fric$global.SES,
             fric.realm       = fric$realm.SES,
             fric.realm.noMDG = fric$realm.SES.noMDG,
             fric.adf         = fric$adf.SES,
             fdis.observed    = fdis$observed,
             fdis.global      = fdis$global.SES,
             fdis.realm       = fdis$realm.SES,
             fdis.realm.noMDG = fdis$realm.SES.noMDG,
             fdis.adf      = fdis$adf.SES,
             row.names = rownames(cwm.traits)
             )


# Init output table
# -----------------
div.results <-
  matrix(data = NA,
         nrow = (ncol(div.data) - 1) * 5,
         ncol = 8,
         dimnames = list(NULL,
                         c("Response", "Coefficient", "Estimate", "P",
                           "PsRsq.KC", "PsRsq.realm", "resid.SAC.P", "bptest.P")
                         )
         ) %>%
  as.data.frame()
div.results[, "Response"] <- rep(colnames(div.data)[-1], each = 5)


# Fit SAR-error model and extract results
# ---------------------------------------
for (predictor in seq_along(colnames(div.data)[-1])) {
  # Fit SAR model with SOI neighbourhood spatial weights matrix
  # We standardize the response variable to mean 0 and unit variance
  cat(colnames(div.data)[predictor + 1], "...\n")
  mod <- FitGlobalSAR(response = div.data[, (predictor + 1)],
                      predictors = div.data[, "realm", drop = FALSE],
                      tdwg.map = tdwg.map,
                      dist.weight = TRUE,
                      double.std = TRUE,
                      numeric.only = FALSE
                      )
  # Extract coefficient estimate and P-value. We have five coefficients:
  # Lambda, intercept, the null realm (new world), OWEast andOWWest
  rows <- (1:5 + 5 * (predictor - 1))
  div.results[rows, "Coefficient"] <-
    c("Lambda", "Intercept", "New World", "Old World West", "Old World East")
  div.results[rows, "Estimate"] <-
    c(mod$lambda, mod$coefficients[1], NA, mod$coefficients[3], mod$coefficients[2])
  div.results[rows, "P"] <-
    c(summary(mod)$LR1$p.value[[1]],
      summary(mod)$Coef[1, 4],
      NA,
      summary(mod)$Coef[3, 4],
      summary(mod)$Coef[2, 4]
      )
  # Calculate pseudo-R-squared
  row <- 1 + (predictor - 1) * 5
  div.results[row, "PsRsq.KC"] <-
    cor(x = predict(mod), y = mod[["y"]], method = "pearson") ^ 2
  div.results[row,"PsRsq.realm"] <- PseudoRsq(mod)
  # Evaluate model validity: spatial autocorrelation in residuals and test
  # for homoskedasticity
  # Note that the 'sar.swmat' object is created by FitGlobalSar()
  div.results[row, "resid.SAC.P"] <-
    moran.test(resid(mod), listw = sar.swmat)[["p.value"]]
  div.results[row, "bptest.P"] <- bptest.sarlm(mod)[["p.value"]]
}


# Save result
# -----------
write.csv(div.results,
          file = "output/SAR_models_realm_effect.csv",
          eol = "\r\n",
          row.names = FALSE
          )





################################
# Species overlap between realms
################################
cat("Comparing species overlap between realms...\n")

# (1) Based on all species in database
# ------------------------------------
# init lists
realms <- levels(tdwg3.info[, "realm"])
areas <- vector("list", length = 3)
names(areas) <- realms
species <- areas
# populate lists
for (i in 1:3) {
  areas[[i]] <- tdwg3.info[, "tdwg3.code"] [tdwg3.info[, "realm"] == realms[i]]
  species[[i]] <- palm.dist[palm.dist[, "Area_code_L3"] %in% areas[[i]], "SpecName"]
}
species <- llply(species, unique)
# Find number of shared species, i.e. species not unique to a realm
shared <- numeric(length = 3)
for (i in 1:3) {
  shared[i] <- sum(species[[i]] %in% unlist(species[-i]))
}


realm.species <-
  data.frame(realm          = c("NewWorld", "OldWorldEast", "OldWorldWest"),
             total.species  = laply(species, length),
             shared.species = shared
             )

write.csv(realm.species,
          file = "output/realm_species_count_shared.csv",
          eol = "\r\n",
          row.names = FALSE
          )


# (2) Based on species included in the final dataset
# --------------------------------------------------
# subset presabs matrix to included TDWG3 units
pres.abs %<>% .[rownames(.) %in% tdwg3.complete, ]
# Find species in each realm
realm.spec <- vector("list", length = 3)
names(realm.spec) <- levels(tdwg3.info$realm)
for (i in seq_along(realm.spec)) {
  tdwg.ind <- tdwg3.info[, "tdwg3.code"] %in% rownames(pres.abs)
  realm.ind <- tdwg3.info[tdwg.ind, "realm"] == names(realm.spec)[i]
  subset <- pres.abs[realm.ind, ]
  subset %<>% .[, colSums(.) > 0]
  realm.spec[[i]] <- colnames(subset)
}
# Find number of shared species
shared <- numeric(length = 3)
for (i in seq_along(realm.spec)) {
  shared[i] <- sum(realm.spec[[i]] %in% unlist(realm.spec[-i]))
}

realm.species.2 <- realm.species
realm.species.2[, "total.species"] <- laply(realm.spec, length)
realm.species.2[, "shared.species"] <- shared

write.csv(realm.species.2,
          file = "output/realm_species_count_shared_final.csv",
          eol = "\r\n",
          row.names = FALSE
          )



cat("Done.\n")

