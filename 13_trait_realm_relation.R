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
library(multcomp)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")


###########
# Load data
###########
cat("preparing data...\n")

# tdwg3 info:
tdwg3.info <- read.csv(file = "output/tdwg3_info_v2.csv")

# traits, based on stochastic gapfilled trait matrix:
cwm.traits <- read.csv(file = "output/observed_FD/community_trait_means_stochastic_mean_of.csv",
                                    row.names = 1
                                    )

#########################################
# Analysis of variance on traits vs realm
#########################################

# Merge trait data with realm data
cwm.traits[, "realm"] <-
  tdwg3.info[match(rownames(cwm.traits), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()

# Run models:
mod.height <- aov(stem.height ~ realm, data = cwm.traits)
mod.blade <- aov(blade.length ~ realm, data = cwm.traits)
mod.fruit <- aov(fruit.length ~ realm, data = cwm.traits)

# Check assumptions:
if (FALSE) {
  # Normality of residuals
  Histogram(resid(mod.height))  # decent
  Histogram(resid(mod.blade))  # perfect
  Histogram(resid(mod.fruit))  # decent
  # Homogeneity of variance
  Scatterplot(x = fitted(mod.height), y = resid(mod.height))  # decent
  Scatterplot(x = fitted(mod.blade), y = resid(mod.blade))  # good
  Scatterplot(x = fitted(mod.fruit), y = resid(mod.fruit))  # decent
}

# Evaluation: p < 0.001 for all three traits. We have ourselves a pattern.

# Post-hoc Tukey tests:
post.height <- glht(model = mod.height, linfct = mcp("realm" = "Tukey"))
post.blade <- glht(model = mod.blade, linfct = mcp("realm" = "Tukey"))
post.fruit <- glht(model = mod.fruit, linfct = mcp("realm" = "Tukey"))

# Evaluation:
#   Stem height: All three realms are significantly different from each other
#   Blade length: OWWest is different from OWEast and NewWorld, but OWEast and
#                 NewWorld do not differ.
#   Fruit length: OWWest is different from OWEast and NewWorld, while OwEast and
#                 NewWorld are only marginally different from each other (p = 0.023)

# For perspective, we might want to look at the realm means:
means <- ddply(cwm.traits,
               "realm",
               summarize,
               mean.stem.height = mean(stem.height),
               mean.blade.length = mean(blade.length),
               mean.fruit.length = mean(fruit.length)
               )
# Keep in mind these are log10-transformed values.


# Visualization of trait differences between realms
# -------------------------------------------------
if (FALSE)
  boxplot(stem.height ~ realm, data = cwm.traits)
  boxplot(blade.length ~ realm, data = cwm.traits)
  boxplot(fruit.length ~ realm, data = cwm.traits)
}
# Looks about like you'd expect.


cat("Done.\n")

