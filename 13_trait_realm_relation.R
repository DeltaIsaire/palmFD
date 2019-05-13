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

# Palm occurrence data
palm.dist <- read.csv(file = "data/palms_in_tdwg3.csv")

# List of TDWG3 units with complete env data, for subsetting.
env <- read.csv(file = "output/tdwg3_predictors_complete_noch.csv",
                row.names = 1
                )
tdwg3.complete <- rownames(env) [complete.cases(env)]


#########################################
# Analysis of variance on traits vs realm
#########################################
cat("Performing anova on traits vs realm...\n")

# Merge trait data with realm data
cwm.traits[, "realm"] <-
  tdwg3.info[match(rownames(cwm.traits), tdwg3.info[, "tdwg3.code"]), "realm"] %>%
  as.factor()
# Subset
cwm.traits %<>% .[rownames(.) %in% tdwg3.complete, ]

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
cat("performing anova on FD indices vs realm...\n")

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


# Run models:
#  ----------
mod.fric.global <- aov(global.SES ~ realm, data = fric)
mod.fric.realm <- aov(realm.SES ~ realm, data = fric)
mod.fric.realm.noMDG <- aov(realm.SES.noMDG ~ realm, data = fdis)

mod.fdis.realm.noMDG <- aov(realm.SES.noMDG ~ realm, data = fdis)


mod.fric.adf <- aov(adf.SES ~ realm, data = fric)
mod.fdis.adf <- aov(adf.SES ~ realm, data = fdis)

mod.fdis.observed <- aov(observed ~ realm, data = fdis)
mod.fric.observed <- aov(observed ~ realm, data = fric)

mod.fdis.global <- aov(global.SES ~ realm, data = fdis)


# Check assumptions:
# ------------------
if (FALSE) {
  # Normality of residuals
  Histogram(resid(mod.fric.global))  # decent
  Histogram(resid(mod.fric.realm))  # decent
  Histogram(resid(mod.fric.realm.noMDG))  # perfect
  Histogram(resid(mod.fric.adf))  # decent

  Histogram(resid(mod.fdis.observed))  # Good

  # Homogeneity of variance
  Scatterplot(x = fitted(mod.fric.global), y = resid(mod.fric.global))  # good
  Scatterplot(x = fitted(mod.fric.realm), y = resid(mod.fric.realm))  # decent
  Scatterplot(x = fitted(mod.fric.realm.noMDG),
              y = resid(mod.fric.realm.noMDG)
              )  # decent
  Scatterplot(x = fitted(mod.fric.adf), y = resid(mod.fric.adf))  # good

  Scatterplot(x = fitted(mod.fdis.observed), y = resid(mod.fdis.observed))  # good
}


# Evaluation:
# -----------
# summary(mod)
# p < 0.001 for all cases. We have ourselves a pattern.


# Post-hoc Tukey tests:
# ---------------------
post.fric.global <- glht(model = mod.fric.global, linfct = mcp("realm" = "Tukey"))
post.fdis.global <- glht(model = mod.fdis.global, linfct = mcp("realm" = "Tukey"))

post.fric.realm <- glht(model = mod.fric.realm, linfct = mcp("realm" = "Tukey"))
post.fric.realm.noMDG <-
  glht(model = mod.fric.realm.noMDG, linfct = mcp("realm" = "Tukey"))

post.fdis.realm.noMDG <-
  glht(model = mod.fdis.realm.noMDG, linfct = mcp("realm" = "Tukey"))

post.fric.adf <- glht(model = mod.fric.adf, linfct = mcp("realm" = "Tukey"))
post.fdis.adf <- glht(model = mod.fdis.adf, linfct = mcp("realm" = "Tukey"))

post.fdis.observed <-
  glht(model = mod.fdis.observed, linfct = mcp("realm" = "Tukey"))
post.fric.observed <-
  glht(model = mod.fric.observed, linfct = mcp("realm" = "Tukey"))

# Post-hoc Evaluation:
# --------------------
# summary(post.mod)
#   fric.global: NewWorld significantly different from other 2 realms (p < 0.001),
#                marginal difference between OWWest and OWEast (p = 0.0417
#   fric.realm: NewWorld significantly different from other two realm (p < 0.001),
#               no difference between OWWest and OWEast
#   fric.realm.noMDG: same story as above
#
#   fdis.observed: OWWest significantly different from NewWorld and OWEast
#                  (p < 0.001), marginal difference between OWEast and NewWorld
#                  (p = 0.045)


# Visualize the results
# ---------------------
# First prepare FDis data
fdis.subset <- fdis[, c("observed", "realm")]
fdis.subset[, "realm"] %<>% as.character()
fdis.noMDG <- fdis.subset[fdis.subset$realm == "OWWest", ]
fdis.noMDG %<>% .[!rownames(.) == "MDG", ]
fdis.noMDG[, "realm"] <- rep("OWWest.noMDG", nrow(fdis.noMDG))
fdis.subset <- rbind(fdis.subset, fdis.noMDG)

# then plot with CWM_Plot()
ggsave(CWM_Plot(x = fdis.subset[complete.cases(fdis.subset), ],
                x.var = "realm",
                y.var = "observed",
                title = NULL,
                subtitle = "Functional Dispersion (FDis) comparison between realms",
                xlab = "Realm",
                ylab = "FDis"
                ),
       filename = "graphs/anova_realm_vs_FDis_observed.png",
       width = 5,
       height = 4
       )

# Next, plots for FRic
for (null.mod in c("global", "realm", "adf")) {
  plot.fric <- CWM_Plot(fric[complete.cases(fric), ],
                        y.var = paste0(null.mod, ".SES"),
                        title = "Functional Richness (FRic)",
                        subtitle = paste("Null model:", null.mod),
                        ylab = "FRic"
                        )
  ggsave(plot = plot.fric,
         filename = paste0("graphs/anova_realm_vs_FRic_", null.mod, ".png"),
         width = 4,
         height = 4
         )
}

# And a special plot for realm
fric.subset <- fric[, ]
fric.subset[, "realm"] %<>% as.character()
fric.subset[, "realm.SES.dual"] <- fric.subset[, "realm.SES"]
fric.noMDG <- fric.subset[fric.subset$realm == "OWWest", ]
fric.noMDG[, "realm.SES.dual"] <- fric.noMDG[, "realm.SES.noMDG"]
fric.noMDG %<>% .[!rownames(.) == "MDG", ]
fric.noMDG[, "realm"] <- rep("OWWest.noMDG", nrow(fric.noMDG))
fric.subset <- rbind(fric.subset, fric.noMDG)

ggsave(CWM_Plot(x = fric.subset[complete.cases(fric.subset), ],
                x.var = "realm",
                y.var = "realm.SES.dual",
                title = NULL,
                subtitle = "Functional Richness (FRic) comparison between realms",
                xlab = "Realm",
                ylab = "FRic"
                ),
       filename = "graphs/anova_realm_vs_FRic_realm_noMDG.png",
       width = 5,
       height = 4
       )



################################
# Species overlap between realms
################################
cat("Comparing species overlap between realms...\n")

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


cat("Done.\n")

