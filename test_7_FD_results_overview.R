##########################################################
# Palm FD project: Summary and visualization of FD results
##########################################################
#
# In which the calculated FD indices are summarized and visualized.
#
# Input files:
#   output/test/test_fd_indices_filled.csv
#   output/test/test_fd_indices_unfilled.csv
#   output/test/tdwg3_info.csv
#   output/test/test_palm_tdwg3_pres_abs_filled.csv
#   output/test/test_palm_tdwg3_pres_abs_unfilled.csv
# Generated output files:
#


cat("Loading required packages and functions...\n")
library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


# ---------------
# Load FD results
# ---------------
cat("Preparing data...\n")

# FD indices based on all traits:
fd.filled <- read.csv(file = "output/test/test_fd_indices_filled.csv")
fd.unfilled <- read.csv(file = "output/test/test_fd_indices_unfilled.csv")
# FD indices based on single traits:
fd.single.filled <- 
  read.csv(file = "output/test/test_fd_indices_single_traits_filled.csv")
fd.single.unfilled <-
  read.csv(file = "output/test/test_fd_indices_single_traits_unfilled.csv")

# Info data on botanical countries:
tdwg3.info <- read.csv(file = "output/test/tdwg3_info.csv")
# The palm.richness column in this dataframe is the complete known palm richness.
# However, the datasets used to calculate FD indices use only a subset of all
# palm species. Consequently, the effective palm richness for the purpose of
# FD calculation differs depending on the trait matrix used.
#
# Hence, a separate dataframe of palm.richness.
# Effective richness can be extracted from the presence/absence matrices.
pres.abs.filled <-
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_filled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
pres.abs.unfilled <-
  read.csv(file = "output/test/test_palm_tdwg3_pres_abs_unfilled.csv",
           row.names = 1,
           check.names = FALSE
           ) %>%
  as.matrix()
# Number of tdwg3 units differs between sources, so some magic is required.
palm.richness <- merge(x = tdwg3.info[, c("tdwg3.code", "palm.richness")],
                       y = data.frame(tdwg3.code = rownames(pres.abs.filled),
                                      filled     = rowSums(pres.abs.filled)
                                      ),
                       by = "tdwg3.code",
                       all.x = TRUE
                       )
palm.richness <- merge(x = palm.richness,
                       y = data.frame(tdwg3.code = rownames(pres.abs.unfilled),
                                      unfilled   = rowSums(pres.abs.unfilled)
                                      ),
                       by = "tdwg3.code",
                       all.x = TRUE
                       )
# Known absences, hence NAs can be set to 0.
palm.richness[is.na(palm.richness)] <- 0
#
# Richness can/should also be added to the fd dataframes:
fd.filled[, "richness"] <- rowSums(pres.abs.filled)
fd.single.filled[, "richness"] <- rowSums(pres.abs.filled)
fd.unfilled[, "richness"] <- rowSums(pres.abs.unfilled)
fd.single.unfilled[, "richness"] <- rowSums(pres.abs.unfilled)


#############################################
# Functional Richness of palms in tdwg3 units
#############################################

# Correlation of single trait FRic (range) with palm richness
# -----------------------------------------------------------
fric.list <- 
  list(all.filled      = fd.filled[, c("FRic", "richness")],
       all.unfilled    = fd.unfilled[, c("FRic", "richness")],
       height.filled   = fd.single.filled[, c("height.FRic", "richness")],
       blade.filled    = fd.single.filled[, c("blade.FRic", "richness")],
       fruit.filled    = fd.single.filled[, c("fruit.FRic", "richness")],
       height.unfilled = fd.single.unfilled[, c("height.FRic", "richness")],
       blade.unfilled  = fd.single.unfilled[, c("blade.FRic", "richness")],
       fruit.unfilled  = fd.single.unfilled[, c("fruit.FRic", "richness")]
       )
fric.mods <- llply(fric.list,
                   function(x) {
                     lm(paste(names(x)[1], "~", names(x)[2]),
                        data = x
                        )
                   }
                   )
fric.mods.summary <-
  data.frame(FRic.trait = names(fric.list),
             p.value    = laply(fric.mods,
                                function(x) {
                                  summary(x)$coefficients[8]
                                }
                                ),
             R.squared  = laply(fric.mods,
                                function(x) {
                                  summary(x)$r.squared
                                }
                                )
             )
# All significant, but much lower R-squared for the single-trait fric compared to
# all-trait fric.
# Not much of a difference between filled and unfilled.

# graphs of correlations
# ----------------------
# Gap-filled data:
Scatterplot(x = fd.filled$richness,
            y = fd.filled$FRic,
            xlab = "Species richness",
            ylab = "Functional richness"
            )
lines(x = fd.filled$richness, y = predict(mod.richness.filled))
# Unfilled data:
Scatterplot(x = fd.unfilled$richness,
            y = fd.unfilled$FRic,
            xlab = "Species richness",
            ylab = "Functional richness"
            )
lines(x = fd.unfilled$richness, y = predict(mod.richness.unfilled))
# WORD OF CAUTION: From these graphs the correlation appears to be curvilinear,
# not linear. There appears to be a saturation effect, which makes sense.

# Single traits:
singlesPlot <- function() {
  

# A single trait sample plot:
Scatterplot(x = fd.single.filled$richness,
            y = fd.single.filled$height.FRic,
            xlab = "Species richness",
            ylab = "Functional richness (stem height)"
            )
lines(x = fd.single.filled$richness, y = predict(fric.mods[[3]]))




cat("Done.\n")

