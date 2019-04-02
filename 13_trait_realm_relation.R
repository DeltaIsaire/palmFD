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
library(ggplot2)
library(gridExtra)

theme_set(theme_bw())

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

# all traits FD indices:
fdis <- read.csv(file = "output/FD_summary_FDis.all.traits.csv", row.names = 1)
fric <- read.csv(file = "output/FD_summary_FRic.all.traits.csv", row.names = 1)


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
CWM_Plot <- function(x, x.var = "realm", y.var, title = NULL, subtitle = NULL,
                     xlab = x.var, ylab = paste0("Log10 (", y.var, ")")) {
  # Plot the community mean of each TDWG3 unit as a dot, then superimpose a boxplot.
  # In both cases, the points are split by realm.
  y.max <- max(x[, y.var])
  ggplot(x) +
    geom_boxplot(aes_string(x = x.var, y = y.var),
                 outlier.size = 0,
                 outlier.color = "white"
                 ) +
    geom_point(aes_string(x = x.var, y = y.var), pch = 21, size = 2) +
    expand_limits(y = 0) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab)
}

for (trait in c("stem.height", "blade.length", "fruit.length")) {
  ggsave(plot = CWM_Plot(cwm.traits,
                         x.var = "realm",
                         y.var = trait,
                         xlab = "Realm",
                         ylab = paste0("Log10 (", trait, ")")
                         ),
         filename = paste0("graphs/anova_", trait, "_vs_realm.png"),
         width = 4,
         height = 4
         )
}


stem.height <-  CWM_Plot(cwm.traits, y.var = "stem.height", title = "Stem height")
blade.length <- CWM_Plot(cwm.traits, y.var = "blade.length", title = "Blade length")
fruit.length <- CWM_Plot(cwm.traits, y.var = "fruit.length", title = "Fruit length")

ggsave(plot = arrangeGrob(stem.height, blade.length, fruit.length, ncol = 3),
       filename = "graphs/anova_traits_vs_realm.png",
       width = 9,
       height = 4
       )

# Looks about like you'd expect.



#############################################
# Analysis of variance on FD indices vs realm
#############################################

# Merge FD data with realm data
# -----------------------------
fdis[, "realm"] <-
  tdwg3.info[match(rownames(fdis), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()
# Merge trait data with realm data
fric[, "realm"] <-
  tdwg3.info[match(rownames(fric), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()

# Run models:
#  ----------
mod.fdis.global <- aov(global.SES ~ realm, data = fdis)
mod.fric.global <- aov(global.SES ~ realm, data = fric)
mod.fdis.realm <- aov(realm.SES ~ realm, data = fdis)
mod.fric.realm <- aov(realm.SES ~ realm, data = fric)
mod.fdis.adf <- aov(adf.SES ~ realm, data = fdis)
mod.fric.adf <- aov(adf.SES ~ realm, data = fric)

# Check assumptions:
if (FALSE) {
  # Normality of residuals
  Histogram(resid(mod.fdis.global))  # perfect
  Histogram(resid(mod.fric.global))  # decent
  Histogram(resid(mod.fdis.realm))  # decent
  Histogram(resid(mod.fric.realm))  # decent
  Histogram(resid(mod.fdis.adf))  # decent
  Histogram(resid(mod.fric.adf))  # decent
  # Homogeneity of variance
  Scatterplot(x = fitted(mod.fdis.global), y = resid(mod.fdis.global))  # good
  Scatterplot(x = fitted(mod.fric.global), y = resid(mod.fric.global))  # good
  Scatterplot(x = fitted(mod.fdis.realm), y = resid(mod.fdis.realm))  # good
  Scatterplot(x = fitted(mod.fric.realm), y = resid(mod.fric.realm))  # decent
  Scatterplot(x = fitted(mod.fdis.adf), y = resid(mod.fdis.adf))  # decent
  Scatterplot(x = fitted(mod.fric.adf), y = resid(mod.fric.adf))  # good
}

# Evaluation: p < 0.001 for both FDis and FRic. We have ourselves a pattern.
# Post-hoc Tukey tests:
post.fdis.global <- glht(model = mod.fdis.global, linfct = mcp("realm" = "Tukey"))
post.fric.global <- glht(model = mod.fric.global, linfct = mcp("realm" = "Tukey"))
post.fdis.realm <- glht(model = mod.fdis.realm, linfct = mcp("realm" = "Tukey"))
post.fric.realm <- glht(model = mod.fric.realm, linfct = mcp("realm" = "Tukey"))
post.fdis.adf <- glht(model = mod.fdis.adf, linfct = mcp("realm" = "Tukey"))
post.fric.adf <- glht(model = mod.fric.adf, linfct = mcp("realm" = "Tukey"))

# Evaluation:
#   FDis: OWWest is significantly different from OWEast and NewWorld, but OWEast
#         and NewWorld do not differ.
#   FRic: NewWorld is significantly different from OWWest and OWEast, but OWWest
#         and OWEast are only marginally different (p = 0.042)

# Visualize the results
# ---------------------
for (null.mod in c("global", "realm", "adf")) {
  plot.fric <- CWM_Plot(fric[complete.cases(fric), ],
                        y.var = paste0(null.mod, ".SES"),
                        title = "Functional Richness (FRic)",
                        subtitle = paste("Null model:", null.mod),
                        ylab = "FRic"
                        )

  plot.fdis <- CWM_Plot(fdis[complete.cases(fdis), ],
                        y.var = paste0(null.mod, ".SES"),
                        title = "Functional Dispersion (FDis)",
                        subtitle = paste("Null model:", null.mod),
                        ylab = "FDis"
                        )

  ggsave(plot = arrangeGrob(plot.fric, plot.fdis, ncol = 2),
         filename = paste0("graphs/anova_FD_vs_realm_", null.mod, ".png"),
         width = 6,
         height = 4
         )
}

cat("Done.\n")

