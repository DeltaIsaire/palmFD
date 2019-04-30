########################################
# Palm FD project: Graphs of all results
########################################
#
# One single script to generate all relevant graphs of the final results.
#
#
# Input files:
#   d
# Generated output files:
#   d


cat("Loading required packages and functions...\n")
library(plyr)
library(magrittr)
library(car)
library(sf)
library(ggplot2)
library(scales)
library(viridis)
library(corrplot)
library(leaps)
library(gridExtra)
library(spdep)
library(ncf)

theme_set(theme_bw())

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")
source(file = "functions/plotting_functions_ggplot2.R")
source(file = "functions/OLS_regression_functions.R")
source(file = "functions/SAR_regression_functions.R")


# Graphics output directory (with trailing slash!)
plot.dir <- "graphs/main_results/"


if (!dir.exists(plot.dir)) {
  cat("creating directory:", plot.dir, "\n")
  dir.create(plot.dir)
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


# Community trait means
# ---------------------
Temp <- function(x) {
  x[, "tdwg3.code"] <- rownames(x)
  merge(x = tdwg3.info[, "tdwg3.code", drop = FALSE],
        y = x,
        by = "tdwg3.code",
        all.x = TRUE,
        sort = TRUE
        )
}
cwm.observed <- 
  read.csv(file = "output/observed_FD/community_trait_means_stochastic_mean_of.csv",
           row.names = 1
           ) %>%
  Temp()
# Keep in mind these are log10-transformed (after adding one):
#   trans = log10(data + 1)
# The reverse operation would be
#   original = (10 ^ trans) - 1


# Environmental data
# ------------------
env.complete <- read.csv(file = "output/tdwg3_predictors_complete.csv",
                         row.names = 1
                         )

# Palm species presence/absence matrix
# ------------------------------------
pres.abs <- read.csv(file = "output/palm_tdwg3_pres_abs_gapfilled.csv",
                     row.names = 1
                     )


# Species trait values
# --------------------
palm.traits <- read.csv(file = "output/trait_matrices/palm_trait_matrix_filled_stochastic_mean_of.csv",
                        row.names = 1)



########################################################
cat("Plotting correlograms of predictor variables...\n")
########################################################

# 1. Full set of predictors, including species, endemism and realm
GraphSVG(Correlogram(env.complete),
         file = paste0(plot.dir, "Correlogram_predictors_ALL.svg"),
         width = 7,
         height = 7
         )

# 2. All environmental predictors
exclude <- c("palm.richness", "endemism", "realm")
GraphSVG(Correlogram(env.complete[, !colnames(env.complete) %in% exclude]),
         file = paste0(plot.dir, "Correlogram_predictors_ENV.svg"),
         width = 7,
         height = 7
         )

# 3. Environmental predictors excluding collinear variables,
#    i.e. without alt_range
exclude <- c("palm.richness", "endemism", "realm", "alt_range")
GraphSVG(Correlogram(env.complete[, !colnames(env.complete) %in% exclude]),
         file = paste0(plot.dir, "Correlogram_predictors_ENV_noncoll.svg"),
         width = 7,
         height = 7
         )



#######################################################
cat("Plotting community mean trait distributions...\n")
#######################################################

MultiTraitPlot <- function(tdwg.map, cwm.observed) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  MakePlot <- function(trait, name, title) {
    SpatialPlot(tdwg.map = tdwg.map[no.ANT, ],
                vector = cwm.observed[no.ANT, trait],
                vector.name = name,
                title = title,
                legend.position = "bottom",
                labels = TraitTrans,
                colors = "inferno"
                ) +
                theme(legend.key.width = unit(1.5, "cm"))
  }

  height <- MakePlot("stem.height",
                     "Stem height (m)",
                     "Mean stem height in palm communities"
                     )
  blade <- MakePlot("blade.length",
                    "Blade length (m)",
                    "Mean blade length in palm communities"
                    )
  fruit <- MakePlot("fruit.length",
                    "Fruit length (cm)",
                    "Mean fruit length in botanical countries"
                    )

  arrangeGrob(height, blade, fruit, ncol = 1)
}

# Full resolution
ggsave(plot = MultiTraitPlot(tdwg.map, cwm.observed),
       filename = paste0(plot.dir, "TDWG3_palm_cwm_traits.png"),
       width = 7,
       height = 12
       )
# Low resolution
ggsave(plot = MultiTraitPlot(tdwg.map, cwm.observed),
       filename = paste0(plot.dir, "TDWG3_palm_cwm_traits_lowres.png"),
       width = 7,
       height = 12,
       dpi = 100
       )



#######################################################
cat("Plotting functional diversity distributions...\n")
#######################################################

# First the all-traits data. 
# A single huge grid. Two columns: FDis and FRic.
# Five rows: observed, global, realm, realmnoMDG, ADF.
MultiFDPlot <- function(tdwg.map, fd.indices) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  MakePlot <- function(index, case, name, title) {
    SpatialPlot(tdwg.map = tdwg.map[no.ANT, ],
                vector = GetFD(fd.indices[index], case) [no.ANT, 1],
                vector.name = name,
                title = title,
                legend.position = "bottom",
                colors = "inferno"
                ) +
                theme(legend.key.width = unit(1.5, "cm"))
  }

  fdis.obs <- MakePlot("FDis.all.traits",
                       "observed",
                       "FDis (observed)",
                       "Functional Dispersion (observed)"
                       )
  fdis.global <- MakePlot("FDis.all.traits",
                          "global.SES",
                          "FDis (SES)",
                          "Functional Dispersion (global)"
                          )
  fdis.realm <- MakePlot("FDis.all.traits",
                         "realm.SES",
                         "FDis (SES)",
                         "Functional Dispersion (realm)"
                         )
  fdis.realm.noMDG <- MakePlot("FDis.all.traits",
                               "realm.SES.noMDG",
                               "FDis (SES)",
                               "Functional Dispersion (realm, removed Madagascar)"
                               )
  fdis.adf <- MakePlot("FDis.all.traits",
                       "adf.SES",
                       "FDis (SES)",
                       "Functional Dispersion (ADF)"
                       )

  fric.obs <- MakePlot("FRic.all.traits",
                       "observed",
                       "FRic (observed)",
                       "Functional Richness (observed)"
                       )
  fric.global <- MakePlot("FRic.all.traits",
                          "global.SES",
                          "FRic (SES)",
                          "Functional Richness (global)"
                          )
  fric.realm <- MakePlot("FRic.all.traits",
                         "realm.SES",
                         "FRic (SES)",
                         "Functional Richness (realm)"
                         )
  fric.realm.noMDG <- MakePlot("FRic.all.traits",
                               "realm.SES.noMDG",
                               "FRic (SES)",
                               "Functional Richness (realm, removed Madagascar)"
                               )
  fric.adf <- MakePlot("FRic.all.traits",
                       "adf.SES",
                       "FRic (SES)",
                       "Functional Richness (ADF)"
                       )

  arrangeGrob(fdis.obs, fric.obs,
              fdis.global, fric.global,
              fdis.realm, fric.realm,
              fdis.realm.noMDG, fric.realm.noMDG,
              fdis.adf, fric.adf,
              ncol = 2
              )
}

# Full resolution
ggsave(plot = MultiFDPlot(tdwg.map, fd.indices),
       filename = paste0(plot.dir, "TDWG3_FD_distributions.png"),
       width = 12,
       height = 18
       )
# Low resolution
ggsave(plot = MultiFDPlot(tdwg.map, fd.indices),
       filename = paste0(plot.dir, "TDWG3_FD_distributions_lowres.png"),
       width = 12,
       height = 18,
       dpi = 100
       )



######################################################
cat("Plotting Sphere of Influence neighbourhood...\n")
######################################################

# The global SOI neighbourhood, and partial SOI neighbourhood for each realm.
# This is for posterity: the neighbourhoods of realms do not overlap.
# The global neighbourhood is saved seperately as well.

# List with the tdwg3 units in each realm
realm.tdwg3 <- vector("list", length = 3)
names(realm.tdwg3) <- levels(tdwg3.info[, "realm"])
for (i in 1:3) {
  realm.tdwg3[[i]] <-
    tdwg3.info[tdwg3.info[, "realm"] == names(realm.tdwg3)[i], "tdwg3.code"]
}

MultiSOIPlot <- function(tdwg.map, env.complete, realm.tdwg3) {

  tdwg.map.subset <- tdwg.map[complete.cases(env.complete), ]
  nb.soi <- SoiNB(tdwg.map.subset)
  global <- NbPlot(nb.soi,
                   tdwg.map,
                   tdwg.map.subset,
                   title = "Global dataset",
                   subtitle = "Sphere of Influence neighbourhood",
                   filename = paste0(plot.dir, "NB_SOI_global.png")
                   )

  tdwg.newworld <-
    tdwg.map.subset[tdwg.map.subset$LEVEL_3_CO %in% realm.tdwg3[["NewWorld"]], ]
  tdwg.owwest <-
    tdwg.map.subset[tdwg.map.subset$LEVEL_3_CO %in% realm.tdwg3[["OWWest"]], ]
  tdwg.oweast <-
    tdwg.map.subset[tdwg.map.subset$LEVEL_3_CO %in% realm.tdwg3[["OWEast"]], ]

  newworld <- NbPlot(SoiNB(tdwg.newworld),
                     tdwg.map,
                     tdwg.newworld,
                     title = "New World",
                     subtitle = "Sphere of Influence neighbourhood",
                     filename = NULL
                     )

  owwest <- NbPlot(SoiNB(tdwg.owwest),
                   tdwg.map,
                   tdwg.owwest,
                   title = "Old World West",
                   subtitle = "Sphere of Influence neighbourhood",
                   filename = NULL
                   )

  oweast <- NbPlot(SoiNB(tdwg.oweast),
                     tdwg.map,
                     tdwg.oweast,
                     title = "Old World East",
                     subtitle = "Sphere of Influence neighbourhood",
                     filename = NULL
                     )

  arrangeGrob(global, newworld, owwest, oweast, ncol = 2)
}

# Full resolution
ggsave(plot = MultiSOIPlot(tdwg.map, env.complete, realm.tdwg3),
       filename = paste0(plot.dir, "NB_SOI_realms.png"),
       width = 12,
       height = 6
       )
# Low resolution
ggsave(plot = MultiSOIPlot(tdwg.map, env.complete, realm.tdwg3),
       filename = paste0(plot.dir, "NB_SOI_realms_lowres.png"),
       width = 12,
       height = 6,
       dpi = 100
       )



################################################################
cat("Plotting trait values as a function of species range...\n")
################################################################

# Plots where x is the # of TDWG3 units in which a species occurs,
# and y is the trait value of that species.
# This will require some magic.

# Prepare species 'ranges' and trait data
# ---------------------------------------
# We use the presence/absence matrix loaded above, subsetted to the TDWG3 units for
# which we have complete data
tdwg3.complete <- rownames(env.complete)[complete.cases(env.complete)]
pres.abs.complete <- pres.abs[rownames(pres.abs) %in% tdwg3.complete, ]
pres.abs.complete %<>% .[, !colSums(.) == 0]
# The trait matrix is loaded above.
# However, we need to fix some name mismatches:
rownames(palm.traits) <-
  gsub(pattern = "-", replacement = ".", x = rownames(palm.traits))
# We can directly construct a dataframe:
spec.traits <-
  data.frame(species = colnames(pres.abs.complete),
             range   = colSums(pres.abs.complete),
             palm.traits[rownames(palm.traits) %in% colnames(pres.abs.complete), ]
             )
# NOTE: this trait data is log10-transformed!

DoPlot <- function(spec.traits) {
  par(mfrow = c(1, 3))
  Scatterplot(x = spec.traits[, "range"],
              y = spec.traits[, "stem.height"],
              xlab = "no. of TWDG3 units occupied",
              ylab = "Stem height ( log10(m) )",
              title = "Stem Height vs species range"
              )
  Scatterplot(x = spec.traits[, "range"],
              y = spec.traits[, "blade.length"],
              xlab = "no. of TWDG3 units occupied",
              ylab = "Blade length ( log10(m) )",
              title = "Blade length vs species range"
              )
  Scatterplot(x = spec.traits[, "range"],
              y = spec.traits[, "fruit.length"],
              xlab = "no. of TWDG3 units occupied",
              ylab = "Fruit length ( log10(cm) )",
              title = "Fruit length vs species range"
              )
}

GraphPNG(DoPlot(spec.traits),
         file = paste0(plot.dir, "traits_vs_species_range.png"),
         width = 1800,
         height = 600,
         pointsize = 28
         )



#####################################################
cat("Plotting comparisons of FD between realms...\n")
#####################################################

# Function for plotting
# ---------------------
RealmPlot <- function(data, x.var, y.var, title = NULL, subtitle = NULL,
                      xlab = NULL, ylab = NULL, breaks = NBreaks,
                      labels = TraitTrans) {
  # Plot the community mean of each TDWG3 unit as a dot, then superimpose a boxplot.
  # In both cases, the points are split by realm.
  y.max <- max(data[, y.var])
  ggplot(data) +
    geom_jitter(aes_string(x = x.var, y = y.var),
                pch = 21,
                bg = "#D0D0D0",
                size = 2,
                height = 0,
                width = 0.15,
                alpha = 0.8
                ) +
    geom_boxplot(aes_string(x = x.var, y = y.var),
                 outlier.size = 0,
                 outlier.color = alpha("white", 0.01),
                 alpha = 0.7,
                 lwd = 0.75
                 ) +
    expand_limits(y = 0) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    scale_y_continuous(breaks = breaks, labels = labels)
}


# 1. Realm vs community trait means 
# ---------------------------------
# Merge trait data with realm data
cwm.observed[, "realm"] <-
  tdwg3.info[match(cwm.observed[, "tdwg3.code"],
                   tdwg3.info[, "tdwg3.code"]
                   ),
             "realm"
             ] %>%
  as.factor()
cwm.data <- cwm.observed[complete.cases(cwm.observed), ]

DoPlot <- function(cwm.data) {
  height <- RealmPlot(data = cwm.data,
                      x.var = "realm",
                      y.var = "stem.height",
                      title = "Stem height",
                      xlab = "Realm",
                      ylab = "Stem height (m)"
                      )
  blade <- RealmPlot(data = cwm.data,
                     x.var = "realm",
                     y.var = "blade.length",
                     title = "Blade length",
                     xlab = "Realm",
                     ylab = "Blade length (m)"
                     )
  fruit <- RealmPlot(data = cwm.data,
                     x.var = "realm",
                     y.var = "fruit.length",
                     title = "Fruit length",
                     xlab = "Realm",
                     ylab = "Fruit length (cm)"
                     )

  arrangeGrob(height, blade, fruit, ncol = 3)
}

ggsave(DoPlot(cwm.data),
       filename = paste0(plot.dir, "realm_comparison_traits.png"),
       width = 12,
       height = 4
       )


# 2. FD for different null models
# -------------------------------
GetData <- function(index, case) {
  data <- data.frame(realm = cwm.observed[, "realm"],
                     index = GetFD(fd.indices[index], case)
                     )
  data[complete.cases(data), ]
}

FDRealmPlot <- function(index, case, title, ylab) {
  RealmPlot(data = GetData(index, case),
          x.var = "realm",
          y.var = index,
          title = title,
          xlab = "Realm",
          ylab = ylab,
          labels = waiver()
          )
}

DoPlot <- function() {
  fdis.obs <-
    FDRealmPlot("FDis.all.traits", "observed", "FDis (observed)", "FDis")
  fdis.global <-
    FDRealmPlot("FDis.all.traits", "global.SES", "FDis (global)", "FDis (SES)")
  fdis.realm <-
    FDRealmPlot("FDis.all.traits", "realm.SES", "FDis (realm)", "FDis (SES)")
  fdis.realm.noMDG <-
    FDRealmPlot("FDis.all.traits", "realm.SES.noMDG",
                "FDis (realm, removed Madagascar)", "FDis (SES)"
                )
  fdis.adf <-
    FDRealmPlot("FDis.all.traits", "adf.SES", "FDis (ADF)", "FDis (SES)")

  fric.obs <-
    FDRealmPlot("FRic.all.traits", "observed", "FRic (observed)", "FRic")
  fric.global <-
    FDRealmPlot("FRic.all.traits", "global.SES", "FRic (global)", "FRic (SES)")
  fric.realm <-
    FDRealmPlot("FRic.all.traits", "realm.SES", "FRic (realm)", "FRic (SES)")
  fric.realm.noMDG <-
    FDRealmPlot("FRic.all.traits", "realm.SES.noMDG",
                "FRic (realm, removed Madagascar)", "FRic (SES)"
                )
  fric.adf <-
    FDRealmPlot("FRic.all.traits", "adf.SES", "FRic (ADF)", "FRic (SES)")

arrangeGrob(fdis.obs, fric.obs,
            fdis.global, fric.global,
            fdis.realm, fric.realm,
            fdis.realm.noMDG, fric.realm.noMDG,
            fdis.adf, fric.adf,
            ncol = 2
            )
}

ggsave(DoPlot(),
       filename = paste0(plot.dir, "realm_comparison_FD.png"),
       width = 8,
       height = 16
       )







cat("Done.\n")

