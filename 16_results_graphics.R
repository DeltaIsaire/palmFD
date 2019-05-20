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
# Using the noch case!
env.complete <- read.csv(file = "output/tdwg3_predictors_complete_noch.csv",
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


# Propagate NAs in env.complete to trait data
# -------------------------------------------
cwm.observed[!complete.cases(env.complete), ] <- NA



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
exclude <- c("palm.richness", "endemism", "realm", "alt_range", "bio12_mean",
             "bio1_mean", "mio_Tano", "mio_Pano")
GraphSVG(Correlogram(env.complete[, !colnames(env.complete) %in% exclude]),
         file = paste0(plot.dir, "Correlogram_predictors_ENV_noncoll.svg"),
         width = 7,
         height = 7
         )



#######################################################
cat("Plotting community mean trait distributions...\n")
#######################################################

# The Legend labels are transformed to meters.
# NOTE: this makes the community mean trait values GEOMETRIC means, not
# arithmetic means!
MultiTraitPlot <- function(tdwg.map, cwm.observed) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  MakePlot <- function(trait, name, title) {
    SpatialPlotBins(tdwg.map = tdwg.map[no.ANT, ],
                    vector = cwm.observed[no.ANT, trait],
                    vector.name = name,
                    vector.size = cwm.observed[no.ANT, trait],
                    title = title,
                    legend.position = "bottom",
                    labels = TraitTrans,
                    colors = "viridis"
                    )
  }

  height <- MakePlot("stem.height",
                     "Geometric mean (m)",
                     "Mean stem height in palm communities"
                     )
  blade <- MakePlot("blade.length",
                    "Geometric mean (m)",
                    "Mean blade length in palm communities"
                    )
  fruit <- MakePlot("fruit.length",
                    "Geometric mean (cm)",
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
    fd.index <- GetFD(fd.indices[index], case)
    fd.index[!complete.cases(env.complete), ] <- NA
    fd.index %<>% .[no.ANT, 1]

    SpatialPlotBins(tdwg.map = tdwg.map[no.ANT, ],
                    vector = fd.index,
                    vector.name = name,
                    vector.size = fd.index,
                    title = title,
                    legend.position = "bottom",
                    colors = "viridis"
                    )
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
  y.min <- min(data[, y.var])
  ggplot(data) +
    geom_hline(yintercept = 0, color = "#808080", linetype = "dashed") +
    geom_jitter(aes_string(x = x.var, y = y.var, fill = "realm"),
                pch = 21,
#                bg = "#D0D0D0",
                size = 2,
                height = 0,
                width = 0.15,
                alpha = 0.8,
                show.legend = FALSE
                ) +
    geom_boxplot(aes_string(x = x.var, y = y.var),
                 outlier.size = 0,
                 outlier.color = alpha("white", 0.01),
                 alpha = 0.7,
                 lwd = 0.75
                 ) +
    expand_limits(y = (0 - 0.05 * (y.max - y.min))) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    scale_y_continuous(breaks = breaks, minor_breaks = NULL, labels = labels) +
    scale_fill_viridis(discrete = TRUE, option = "plasma")
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
                      ylab = "A. Stem height (m)"
                      ) +
              annotate(geom = "text",
                       x = c(1, 2, 3),
                       y = 0.05,
                       label = c("A", "B", "C")
                       )
  blade <- RealmPlot(data = cwm.data,
                     x.var = "realm",
                     y.var = "blade.length",
                     title = "Blade length",
                     xlab = "Realm",
                     ylab = "B. Blade length (m)"
                     ) +
              annotate(geom = "text",
                       x = c(1, 2, 3),
                       y = 0.05,
                       label = c("A", "A", "B")
                       )
  fruit <- RealmPlot(data = cwm.data,
                     x.var = "realm",
                     y.var = "fruit.length",
                     title = "Fruit length",
                     xlab = "Realm",
                     ylab = "C. Fruit length (cm)"
                     ) +
              annotate(geom = "text",
                       x = c(1, 2, 3),
                       y = 0.05,
                       label = c("A", "B", "C")
                       )

  arrangeGrob(height, blade, fruit, ncol = 3)
}

ggsave(DoPlot(cwm.data),
       filename = paste0(plot.dir, "realm_comparison_traits.png"),
       width = 12,
       height = 4,
       dpi = 100
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



##############################################################
cat("Plotting distributions of environmental predictors...\n")
##############################################################

MultiEnvPlot <- function(tdwg.map, env.complete) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  MakePlot <- function(predictor, name, title) {
    SpatialPlotBins(tdwg.map = tdwg.map[no.ANT, ],
                    vector = env.complete[, predictor] [no.ANT],
                    vector.name = name,
                    vector.size = env.complete[, predictor] [no.ANT],
                    title = title,
                    legend.position = "bottom",
                    colors = "viridis"
                    )
  }

  soilcount <- MakePlot("soilcount", "no. soiltypes", "Soilcount")
  temp.sd <- MakePlot("bio1_sd", "°C", "Mean annual temperature (sd)")
  precip.sd <- MakePlot("bio12_sd", "mm/year", "Annual Precipitation (sd)")
  temp.seas <-
    MakePlot("bio4_mean", "°C", "Temperature Seasonality")
  precip.seas <-
    MakePlot("bio15_mean", "CV", "Precipitation Seasonality")
  lgmt <- MakePlot("lgm_Tano", "°C", "LGM Temperature Anomaly")
  lgmp <- MakePlot("lgm_Pano", "mm/year", "LGM Precipitation Anomaly")
  pliot <- MakePlot("plio_Tano", "°C", "Pliocene Temperature Anomaly")
  pliop <- MakePlot("plio_Pano", "mm/year", "Pliocene Precipitation Anomaly")
  miot <- MakePlot("mio_Tano", "°C", "Miocene Temperature Anomaly")
  miop <- MakePlot("mio_Pano", "mm/year", "Miocene Temperature Anomaly")

  arrangeGrob(temp.sd, precip.sd,
              temp.seas, precip.seas,
              lgmt, lgmp,
              pliot, pliop,
              miot, miop,
              soilcount,
              ncol = 2
              )
}

# Full resolution
ggsave(plot = MultiEnvPlot(tdwg.map, env.complete),
       filename = paste0(plot.dir, "TDWG3_predictor_distributions.png"),
       width = 12,
       height = 20
       )
# Low resolution
ggsave(plot = MultiEnvPlot(tdwg.map, env.complete),
       filename = paste0(plot.dir, "TDWG3_predictor_distributions_lowres.png"),
       width = 12,
       height = 20,
       dpi = 100
       )


####################################################
cat("Realm division and community trait means...\n")
####################################################

# The Legend labels for trait values are transformed to meters.
# NOTE: this makes the community mean trait values GEOMETRIC means, not
# arithmetic means!
TraitRealmPlot <- function(tdwg.map, cwm.observed) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  threerealm <- SpatialPlotFill(tdwg.map =tdwg.map[no.ANT, ],
                                vector = tdwg3.info[no.ANT, "realm"],
                                vector.name = "Realm",
                                title = "A. Division of TDWG3 units into realms",
                                legend.position = "bottom",
                                colors = "plasma",
                                begin = 0.25
                                )

  MakePlot <- function(trait, name, title) {
    SpatialPlotBins(tdwg.map = tdwg.map[no.ANT, ],
                    vector = cwm.observed[no.ANT, trait],
                    vector.name = name,
                    vector.size = cwm.observed[no.ANT, trait],
                    title = title,
                    legend.position = "bottom",
                    labels = TraitTrans,
                    colors = "viridis"
                    )
  }

  # Add text labels indicating ANOVA post-hoc results
  height <- MakePlot("stem.height",
                     "Geometric mean (m)",
                     "B. Stem height"
                     )
  blade <- MakePlot("blade.length",
                    "Geometric mean (m)",
                    "C. Blade length"
                    )
  fruit <- MakePlot("fruit.length",
                    "Geometric mean (cm)",
                    "D. Fruit length"
                    )

  arrangeGrob(threerealm, height, blade, fruit, ncol = 2)
}

# Full resolution
ggsave(plot = TraitRealmPlot(tdwg.map, cwm.observed),
       filename = paste0(plot.dir, "TDWG3_Trait_Realm.png"),
       width = 12,
       height = 7
       )
# Low resolution
ggsave(plot = TraitRealmPlot(tdwg.map, cwm.observed),
       filename = paste0(plot.dir, "TDWG3_Trait_Realm_lowres.png"),
       width = 12,
       height = 7,
       dpi = 100
       )



###########################################
cat("FD distributions observed + ADF...\n")
###########################################

FDPlot <- function(tdwg.map, fd.indices) {

  no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"

  MakePlot <- function(index, case, name, title) {
    fd.index <- GetFD(fd.indices[index], case)
    fd.index[!complete.cases(env.complete), ] <- NA
    fd.index %<>% .[no.ANT, 1]

    SpatialPlotBins(tdwg.map = tdwg.map[no.ANT, ],
                    vector = fd.index,
                    vector.name = name,
                    vector.size = fd.index,
                    title = title,
                    legend.position = "bottom",
                    colors = "viridis"
                    )
  }

  fdis.obs <- MakePlot("FDis.all.traits",
                       "observed",
                       "FDis (unitless)",
                       "A. Functional Dispersion (observed)"
                       )

  fdis.adf <- MakePlot("FDis.all.traits",
                       "adf.SES",
                       "FDis (SES)",
                       "C. Functional Dispersion (ADF null model)"
                       )

  fric.obs <- MakePlot("FRic.all.traits",
                       "observed",
                       "FRic (unitless)",
                       "B. Functional Richness (observed)"
                       )

  fric.adf <- MakePlot("FRic.all.traits",
                       "adf.SES",
                       "FRic (SES)",
                       "D. Functional Richness (ADF null model)"
                       )

  arrangeGrob(fdis.obs, fric.obs,
              fdis.adf, fric.adf,
              ncol = 2
              )
}

# Full resolution
ggsave(plot = FDPlot(tdwg.map, fd.indices),
       filename = paste0(plot.dir, "TDWG3_FD_observed_ADF.png"),
       width = 12,
       height = 7
       )
# Low resolution
ggsave(plot = FDPlot(tdwg.map, fd.indices),
       filename = paste0(plot.dir, "TDWG3_FD_observed_ADF_lowres.png"),
       width = 12,
       height = 7,
       dpi = 100
       )



##############################################
cat("FD Realm comparison observed + ADF...\n")
##############################################

# For observed + ADF null models
# ------------------------------
DoPlot <- function() {
  fdis.obs <-
    FDRealmPlot("FDis.all.traits", "observed", "A. FDis (observed)", "FDis") +
    annotate(geom = "text", x = c(1, 2, 3), y = 0.06, label= c("A", "B", "C"))

  fdis.adf <-
    FDRealmPlot("FDis.all.traits", "adf.SES", "C. FDis (ADF null model)",
                "FDis (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -7, label= c("A", "B", "AB"))

  fric.obs <-
    FDRealmPlot("FRic.all.traits", "observed", "B. FRic (observed)", "FRic") +
    annotate(geom = "text", x = c(1, 2, 3), y = -5, label= c("AB", "A", "B"))

  fric.adf <-
    FDRealmPlot("FRic.all.traits", "adf.SES", "D. FRic (ADF null model)",
                "FRic (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -5, label= c("A", "B", "A"))

arrangeGrob(fdis.obs, fric.obs,
            fdis.adf, fric.adf,
            ncol = 2
            )
}

ggsave(DoPlot(),
       filename = paste0(plot.dir, "realm_comparison_FD_obs_ADF.png"),
       width = 9,
       height = 8,
       dpi = 100
       )


# For all three null models (using noMDG for realm)
# -------------------------------------------------

# First FRic
DoPlot <- function() {
  fric.global <-
    FDRealmPlot("FRic.all.traits", "global.SES", "A. FRic (global)", "FRic (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -3.4, label = c("A", "B", "B"))

  fric.realm <-
    FDRealmPlot("FRic.all.traits", "realm.SES", "B. FRic (realm)", "FRic (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -4, label = c("A", "B", "C"))

  fric.adf <-
    FDRealmPlot("FRic.all.traits", "adf.SES", "D. FRic (ADF null model)",
                "FRic (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -3.5, label= c("A", "B", "A"))
  
  arrangeGrob(fric.global, fric.realm, fric.adf, ncol = 1)
}

ggsave(DoPlot(),
       filename = paste0(plot.dir, "realm_comparison_FRic_nullmodels.png"),
       width = 4,
       height = 7,
       dpi = 100
       )


# Second FDis
DoPlot <- function() {
  fdis.global <-
    FDRealmPlot("FDis.all.traits", "global.SES", "A. FDis (global)", "FDis (SES)") +
    annotate(geom = "text", x = c(1.2, 2, 3), y = -5, label = c("A", "A", "B"))

  fdis.realm <-
    FDRealmPlot("FDis.all.traits", "realm.SES", "B. FDis (realm)", "FDis (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -6, label = c("A", "B", "C"))

  fdis.adf <-
    FDRealmPlot("FDis.all.traits", "adf.SES", "C. FDis (ADF null model)",
                "FDis (SES)") +
    annotate(geom = "text", x = c(1, 2, 3), y = -6, label= c("A", "B", "AB"))

  arrangeGrob(fdis.global, fdis.realm, fdis.adf, ncol = 1)
}

ggsave(DoPlot(),
       filename = paste0(plot.dir, "realm_comparison_FDis_nullmodels.png"),
       width = 4,
       height = 7,
       dpi = 100
       )



#######################################################
cat("Comparison of all-traits with single-trait FD:\n")
#######################################################

# Function for plotting
# ---------------------
MakePlot <- function(index, case, name, title) {
  fd.index <- GetFD(fd.indices[index], case)
  fd.index[!complete.cases(env.complete), ] <- NA
  fd.index %<>% .[no.ANT, 1]
  SpatialPlotBins(tdwg.map = tdwg.map[no.ANT, ],
                  vector = fd.index,
                  vector.name = name,
                  vector.size = fd.index,
                  title = title,
                  legend.position = "bottom",
                  colors = "viridis"
                  )
}

no.ANT <- !tdwg.map$LEVEL_3_CO == "ANT"


FDPlot <- function(tdwg.map, fd.indices, fd.index, null.model, unit, titles) {
  all.traits <- MakePlot(index = fd.index[1],
                         case = null.model,
                         name = unit,
                         title = titles[1]
                         )
  height <- MakePlot(index = fd.index[2],
                     case = null.model,
                     name = unit,
                     title = titles[2]
                     )
  blade <- MakePlot(index = fd.index[3],
                    case = null.model,
                    name = unit,
                    title = titles[3]
                    )
  fruit <- MakePlot(index = fd.index[4],
                    case = null.model,
                    name = unit,
                    title = titles[4]
                    )

  arrangeGrob(all.traits, height,
              blade, fruit,
              ncol = 2
              )
}

SavePlots <- function(plot, filename, width, height, dpi.lowres) {
  ggsave(plot = plot,
         filename = paste0(filename, ".png"),
         width = width,
         height = height
         )
  # Low resolution
  ggsave(plot = plot,
         filename = paste0(filename, "_lowres.png"),
         width = width,
         height = height,
         dpi = dpi.lowres
         )
}

cat("(1) FRic observed...\n")
# ---------------------------
SavePlots(plot = FDPlot(tdwg.map = tdwg.map,
                        fd.indices = fd.indices,
                        fd.index = fric,
                        null.model = "observed",
                        unit = "FRic (unitless)",
                        titles = c("A. FRic of all traits (observed)",
                                   "B. FRic of stem height (observed)",
                                   "C. FRic of blade length (observed)",
                                   "D. FRic of fruit length (observed)"
                                   )
                        ),
          filename = paste0(plot.dir, "TDWG3_FRic_observed_traits"),
          width = 12,
          height = 7,
          dpi.lowres = 100
          )

cat("(2) FDis observed...\n")
# ---------------------------
SavePlots(plot = FDPlot(tdwg.map = tdwg.map,
                        fd.indices = fd.indices,
                        fd.index = fdis,
                        null.model = "observed",
                        unit = "FDis (unitless)",
                        titles = c("A. FDis of all traits (observed)",
                                   "B. FDis of stem height (observed)",
                                   "C. FDis of blade length (observed)",
                                   "D. FDis of fruit length (observed)"
                                   )
                        ),
          filename = paste0(plot.dir, "TDWG3_FDis_observed_traits"),
          width = 12,
          height = 7,
          dpi.lowres = 100
          )

cat("(3) FRic ADF...\n")
# ---------------------------
SavePlots(plot = FDPlot(tdwg.map = tdwg.map,
                        fd.indices = fd.indices,
                        fd.index = fric,
                        null.model = "adf.SES",
                        unit = "FRic (SES)",
                        titles = c("A. FRic of all traits (ADF null model)",
                                   "B. FRic of stem height (ADF null model)",
                                   "C. FRic of blade length (ADF null model)",
                                   "D. FRic of fruit length (ADF null model)"
                                   )
                        ),
          filename = paste0(plot.dir, "TDWG3_FRic_ADF_traits"),
          width = 12,
          height = 7,
          dpi.lowres = 100
          )

cat("(4) FDis ADF...\n")
# ---------------------------
SavePlots(plot = FDPlot(tdwg.map = tdwg.map,
                        fd.indices = fd.indices,
                        fd.index = fdis,
                        null.model = "adf.SES",
                        unit = "FDis (SES)",
                        titles = c("A. FDis of all traits (ADF null model)",
                                   "B. FDis of stem height (ADF null model)",
                                   "C. FDis of blade length (ADF null model)",
                                   "D. FDis of fruit length (ADF null model)"
                                   )
                        ),
          filename = paste0(plot.dir, "TDWG3_FDis_ADF_traits"),
          width = 12,
          height = 7,
          dpi.lowres = 100
          )



########################################
cat("Multimodel averaging results...\n")
########################################

PlotModAvg <- function(filename, title, subtitle = NULL) {
  predictors <- c("lambda", "(Intercept)", "soilcount", "bio1_sd",
                  "bio12_sd", "bio4_mean", "bio15_mean", "lgm_Tano", "lgm_Pano",
                  "plio_Tano", "plio_Pano", "mio_Tano", "mio_Pano")
  df <- read.csv(file = filename)
  df %<>% { .[match(predictors, .[, "X"]), ] }
  df %<>% { .[!.[, "X"] %in% c("lambda", "(Intercept)"), ] }
  df[, "X"] <- c("Soil", "Temp.sd", "Prec.sd", "T.seas", "P.seas",
                 "lgm.Tano", "lgm.Pano", "pli.Tano", "pli.Pano", "mio.Tano",
                 "mio.Pano")
  df[, 1] %<>% { factor(., levels = .) }

  ggplot(df) +
    geom_hline(yintercept = 0, color = "#808080", linetype = "dashed") +
    geom_linerange(aes(x = X, ymin = X2.5.., ymax = X97.5..),
                   color = "black",
                   size = 1
                   ) +
    geom_point(aes(x = X, y = Estimate), color = "red", size = 3) +
    labs(x = NULL,
         y = "Estimate (+/- 95% CI)",
         title = title,
         subtitle = subtitle
         ) #+
#    theme(axis.text.x = element_text(vjust = unit(c(0, -0.5), "points")))
}

fdis.obs <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FDis.all.traits_observed_full_SAR.csv",
                       title = "A. Functional Dispersion (observed)",
                       subtitle = "Model averaging coefficient estimates"
                       )
fric.obs <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FRic.all.traits_observed_full_SAR.csv",
                       title = "B. Functional Richness (observed)",
                       subtitle = "Model averaging coefficient estimates"
                       )

fdis.global <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FDis.all.traits_global.SES_full_SAR.csv",
                       title = "Functional Dispersion (global)",
                       subtitle = "Model averaging coefficient estimates"
                       )
fric.global <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FRic.all.traits_global.SES_full_SAR.csv",
                       title = "Functional Richness (global)",
                       subtitle = "Model averaging coefficient estimates"
                       )

fdis.realm <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FDis.all.traits_realm.SES.noMDG_full_SAR.csv",
                       title = "Functional Dispersion (realm)",
                       subtitle = "Model averaging coefficient estimates"
                       )
fric.realm <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FRic.all.traits_realm.SES.noMDG_full_SAR.csv",
                       title = "Functional Richness (realm)",
                       subtitle = "Model averaging coefficient estimates"
                       )

fdis.adf <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FDis.all.traits_adf.SES_full_SAR.csv",
                       title = "C. Functional Dispersion (ADF null model)",
                       subtitle = "Model averaging coefficient estimates"
                       )
fric.adf <- PlotModAvg(filename = "output/multimodel_averaging/multimod_avg_FRic.all.traits_adf.SES_full_SAR.csv",
                       title = "D. Functional Richness (ADF null model)",
                       subtitle = "Model averaging coefficient estimates"
                       )

ggsave(plot = arrangeGrob(fdis.obs, fric.obs, fdis.adf, fric.adf, ncol = 2),
       filename = paste0(plot.dir, "mod_avg_full_SAR.png"),
       width = 12,
       height = 7,
       dpi = 100
       )

ggsave(plot = arrangeGrob(fdis.obs, fric.obs, ncol = 2),
       filename = paste0(plot.dir, "mod_avg_full_SAR_observed.png"),
       width = 12,
       height = 4,
       dpi = 100
       )

ggsave(plot = arrangeGrob(fdis.global, fric.global,
                          fdis.realm, fric.realm,
                          fdis.adf, fric.adf,
                          ncol = 2
                          ),
       filename = paste0(plot.dir, "mod_avg_full_SAR_nullmodels.png"),
       width = 12,
       height = 10,
       dpi = 100
       )




cat("Done.\n")

