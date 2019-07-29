library(magrittr)
library(plyr)

source(file = "functions/plotting_functions.R")


# Get names of all palm species, divided by genus
palm.traits <- read.csv(file = "data/PalmTraits_10.csv")
palm.spec <- palm.traits[, c("accGenus", "accSpecies")]
palm.spec$specname <- paste(palm.spec$accGenus, palm.spec$accSpecies, sep = "_")
palm.spec$specname %<>% sub("-", ".", .)

# Get list of species included in final dataset.
# 	(1) get list of included tdwg3 units
env <- read.csv(file = "output/tdwg3_predictors_complete_noch.csv", row.names = 1)
tdwg3.included <- rownames(env)[complete.cases(env)]
#	(2) crosscheck with species occurence data
pres.abs <-
  read.csv(file = "output/palm_tdwg3_pres_abs_gapfilled.csv", row.names = 1) %>%
  .[rownames(.) %in% tdwg3.included, ] %>%
  .[, colSums(.) > 0]

# Subset palm.spec to included species
palm.spec %<>% .[.[, "specname"] %in% colnames(pres.abs), ]

# For each tdwg3 unit, extract present species
speclist <- alply(pres.abs, .margins = 1, function(x) { colnames(x)[x > 0] } )
# Then convert species to genus
genlist <- llply(speclist,
                 function(x) { palm.spec$accGenus[palm.spec$specname %in% x] }
                 )
# Get frequency of each genus in each tdwg3 unit
genfreq <- llply(genlist, count)
names(genfreq) <- rownames(pres.abs)
# convert frequency to proportion
Prop <- function(x) {
  colsum <- sum(x$freq)
  x$fraction <- round(x$freq / colsum, digits = 3)
  x
  }
genfreq <- llply(genfreq, Prop)

# Let's make a Histogram of all the frequencies, in order to select a cutoff
# point for what we consider local radiations
allfreqs <- unlist(llply(genfreq, function(x) { x$freq } ))
# Histogram(allfreqs)
# Histogram(log10(allfreqs))
# yeah this is not working. There are no clear gaps.
# But we could say 15 is a decent cutoff.
summary(allfreqs)
sd(allfreqs)
# mean is 3.5, sd is 7.5.
# we could use a cutoff of two standard deviations above the mean,
# or 3.5 + 2 * 7.5 = 18.5, or 18.
radiations <- llply(genfreq, function(x) { x[x$freq >= 18, ] } )
radiations %<>% .[laply(radiations, nrow) > 0]
# So that is 31 (out of 126) tdwg3 units that have unusually abundant genera

# Let's summarize that
Unwrap <- function(list, col) {
  unlist(llply(list, function(x) { x[, col] } ))
}

genus.rad <- list(genera = unique(Unwrap(radiations, "x")),
                  freqs = Unwrap(radiations, "freq"),
                  fractions = Unwrap(radiations, "fraction")
                  )
# 21 genera. Mean radiation is 32 species per tdwg3 unit, making up 25.4% of all
# palm species in the tdwg3 unit.

# Let's get some genus-level statistics
genus.count <-
  count(palm.spec$accGenus) %>%
  .[order(.$freq, decreasing = TRUE), ]
# 158 included genera, with 2421 species in total
# But the 9 most abundant genera already account for >50% of species (52.3%)
# species per genus is average 15.32 +/- sd of 39.27
# Mean + 2 SD is ~94 species, a threshold reached by six genera:
# Calamus, Dypsis, Licuala, Pinanga, Daemonorops and Chamaedorea
# (together 1051 species or 43.4% of all included species)

# Dypsis on Madagascar is the principal outlier, with 158 species of Dypsis
# accounting for 81% of all palm species in Madagascar.
# In terms of % total species, MDG is the only outlier.
# In terms of number of species per radiation, Borneo is also special with 83
# species of Calamus, and Malaya (Malaysia) is mildly special with 63 species
# of Calamus.

# What these three tdwg3 units have in common is high palm species richness
# (195, 306 and 232 species, respectively) and high endemism (98%, 74% and 39%,
# respectively).
# But only MDG has very low functional diversity.
