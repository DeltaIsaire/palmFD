#####################################################
# Palm FD project: trait_filling_comparison test code
#####################################################
#
# This is a companion file listing the R code used to verify the code used for
# trait filling comparisons, and to explore/verify inputs and outputs.
# The code in this file is non-essential and provided only for reference and
# bug-fixing purposes.
# Input files:
#   none.
# Generated output files:
#   none.


library(plyr)

source(file="2_trait_filling_comparison.R")


# ----------------
# Data preparation
# ----------------
# Find out how many species are not shared between filled trait matrices:
MultiCC(list(mean  = traits.mean$species,
             BHPMF = traits.BHPMF$species),
        presence=FALSE, value=FALSE)


# ---------------------------------
# Comparing 5 gap-filling scenarios
# ---------------------------------
# We begin with the base unimputed trait matrix, which has a single 'oberved'
# trait value for each species, with gaps (NAs).
# Subsequently we have generated:
# 1. Matrix gap-filled using genus means.
filled.one <- read.csv(file="output/palm.traits.genus.mean.csv")
# 2. Matrix gap-filled using BHPMF with our three traits of interest, extracting
#    from the BHPMF output only the estimates for the gaps in the original matrix.
filled.two <- read.csv(file="output/palm.traits.BHPMF.csv")
# 3. Matrix gap-filled using BHPMF with our three traits of interest,
#    using the full output of BHPMF (including estimates of 'observed' trait data)
filled.three <- read.csv(file="output/TEST.BHPMF_filled_3.csv")
# 4. Same as #2, but including a fourth dummy trait, with all values = 1
filled.four <- read.csv(file="output/TEST.BHPMF_filled_4.csv")
# 5. Same as #2, but including the following additional traits:
#    MaxStemDia_cm, MaxLeafNumber, Max_Rachis_Length_m, Max_Petiole_length_m,
#    AverageFruitWidth_cm
filled.five <- read.csv(file="output/TEST.BHPMF_filled_5.csv")

# Completeness
# ------------
# As a first step, compare the number of species in the filled matrices,
# as a measure of completeness.
# For comparison we load the original unfilled matrix as well.
original <- read.csv(file="output/palm.traits.csv")
lengths <- ldply(list(original, filled.one, filled.two, filled.three,
                      filled.four, filled.five),
                 function(x) { length(x[, 1]) } )
rownames(lengths) <- c("original", "one", "two", "three", "four", "five")
names(lengths) <- "species.count"
# Save a plot for visual inspection:
GraphSVG(barplot(lengths[, 1], names.arg=rownames(lengths), ylab="species count"),
         file="graphs/gapfill_completeness.svg", width=6, height=4)

# Boxplots
# --------
# Next, for each trait we can make a boxplot to compare the spread of the
# original and gap-filled trait values.

# First, for each trait, we need the list of shared.species whose trait value
# was estimated. Thankfully we only need to look at the most conservative of
# the BHPMF methods, which is method two, and crossCheck it with genus means
# (method one)
shared.species <- CrossCheck(filled.two$species,
                             filled.one$ species,
                             presence=TRUE, value=TRUE)
# Then, for each trait, we extract the predicted trait values for these species.
# We use the power of functions for the sake of code legibility:
objects <- list(original, filled.one, filled.two, filled.three, filled.four,
                filled.five)
Extract <- function(objects, trait) {
  values <- lapply(objects, function(x) {
                    spec <- CrossCheck(x[, "species"], shared.species,
                                       presence=TRUE, value=FALSE)
                    return (x[, trait][spec])
                    }
                  )
 frame <- data.frame(species = shared.species, values)
 names(frame) <- c("species", "original", "one", "two", "three", "four", "five")
 return (frame)
}
height.estimates <- Extract(objects=objects, trait="stem.height")
blade.estimates <- Extract(objects=objects, trait="blade.length")
fruit.estimates <- Extract(objects=objects, trait="fruit.length")

# Now we can make the boxplots
Boxplot <- function(x, ylab) {
  # set margins to be nicely narrow
  par(mar=c(4.1, 4.1, .5, .5))
  boxplot(log10(x[, -1]), border="black", col="lightgrey", ylab=ylab,
          xlab="log10(values)")
  # Add horizontal line for comparison
  abline(h=median(log10(x[, "original"]), na.rm=TRUE))
}
GraphSVG(Boxplot(x=height.estimates, ylab="log10(stem height)"),
         file="graphs/boxplot.test.height.svg", width=6, height=9)
GraphSVG(Boxplot(x=blade.estimates, ylab="log10(blade length)"),
         file="graphs/boxplot.test.blade.svg", width=6, height=9)
GraphSVG(Boxplot(x=fruit.estimates, ylab="log10(fruit length)"),
         file="graphs/boxplot.test.fruit.svg", width=6, height=9)

# Next, it would be cool to boxplot only the gap-filled values.
# For this we have to subset the data.
# For each trait, subset the trait.estimates dataframes
# to only species with original NA value.
height.estimates.subset <- 
  height.estimates[which(is.na(height.estimates$original)), ]
blade.estimates.subset <- 
  blade.estimates[which(is.na(blade.estimates$original)), ]
fruit.estimates.subset <- 
  fruit.estimates[which(is.na(fruit.estimates$original)), ]
# And produce the boxplots:
GraphSVG(Boxplot(x=height.estimates.subset, ylab="log10(stem height)"),
         file="graphs/boxplot.test.height.subset.svg", width=6, height=9)
GraphSVG(Boxplot(x=blade.estimates.subset, ylab="log10(blade length)"),
         file="graphs/boxplot.test.blade.subset.svg", width=6, height=9)
GraphSVG(Boxplot(x=fruit.estimates.subset, ylab="log10(fruit length)"),
         file="graphs/boxplot.test.fruit.subset.svg", width=6, height=9)

# Scatterplots
# ------------
# 1 on 1 comparison of methods. Most relevant is comparing method 2 (standard
# BHPMF) to all the other methods. That's a lot of graphs. 
# First for stem height:
GraphSVG(MultiScatter(height.estimates$one, height.estimates$two, x.name="mean",
         y.name="BHPMF"),
         file="graphs/scatter.test.height.2.vs.1.svg",
         width=12, height=4)
GraphSVG(MultiScatter(height.estimates$three, height.estimates$two, x.name="three",
         y.name="two"),
         file="graphs/scatter.test.height.2.vs.3.svg",
         width=12, height=4)
GraphSVG(MultiScatter(height.estimates$four, height.estimates$two, x.name="four",
         y.name="two"),
         file="graphs/scatter.test.height.2.vs.4.svg",
         width=12, height=4)
GraphSVG(MultiScatter(height.estimates$five, height.estimates$two, x.name="five",
         y.name="two"),
         file="graphs/scatter.test.height.2.vs.5.svg",
         width=12, height=4)
# Next for Blade length:
GraphSVG(MultiScatter(blade.estimates$one, blade.estimates$two, x.name="mean",
         y.name="BHPMF"),
         file="graphs/scatter.test.blade.2.vs.1.svg",
         width=12, height=4)
GraphSVG(MultiScatter(blade.estimates$three, blade.estimates$two, x.name="three",
         y.name="two"),
         file="graphs/scatter.test.blade.2.vs.3.svg",
         width=12, height=4)
GraphSVG(MultiScatter(blade.estimates$four, blade.estimates$two, x.name="four",
         y.name="two"),
         file="graphs/scatter.test.blade.2.vs.4.svg",
         width=12, height=4)
GraphSVG(MultiScatter(blade.estimates$five, blade.estimates$two, x.name="five",
         y.name="two"),
         file="graphs/scatter.test.blade.2.vs.5.svg",
         width=12, height=4)
# And finally for fruit length:
GraphSVG(MultiScatter(fruit.estimates$one, fruit.estimates$two, x.name="mean",
         y.name="BHPMF"),
         file="graphs/scatter.test.fruit.2.vs.1.svg",
         width=12, height=4)
GraphSVG(MultiScatter(fruit.estimates$three, fruit.estimates$two, x.name="three",
         y.name="two"),
         file="graphs/scatter.test.fruit.2.vs.3.svg",
         width=12, height=4)
GraphSVG(MultiScatter(fruit.estimates$four, fruit.estimates$two, x.name="four",
         y.name="two"),
         file="graphs/scatter.test.fruit.2.vs.4.svg",
         width=12, height=4)
GraphSVG(MultiScatter(fruit.estimates$five, fruit.estimates$two, x.name="five",
         y.name="two"),
         file="graphs/scatter.test.fruit.2.vs.5.svg",
         width=12, height=4)


# --------------------------------
# In-depth gap-filling comparisons
# --------------------------------
# above, we compared 5 gap-filling scenarios.
# Here, we do a few more in-depth comparisons.

# BHPMF: original vs estimated values
# -----------------------------------
# Here we compare the original values (dataframe 'original') to estimates from
# BHPMF (dataframe 'filled.three')
# First prepare the data. For each trait, data for the subset of species whose
# trait value was known/available.
temp <- function(trait) {
  temp <- original[!is.na(original[, trait]), c("species", trait)]
  temp <- data.frame(temp,
                     estimate = filled.three[CrossCheck(filled.three$species,
                                                        temp$species,
                                                        presence=TRUE,
                                                        value=FALSE),
                                             trait]
                     )
  return (temp)
}
height <- temp("stem.height")
blade <- temp("blade.length")
fruit <- temp("fruit.length")

# Then we can compare originals with estimates using MultiScatter:
GraphSVG(MultiScatter(height[, 2], height[, 3],
                      x.name="observed stem height", 
                      y.name="estimated stem height"),
         file="graphs/scatter.test.height.original.vs.estimated.svg",
         width=12, height=4)
GraphSVG(MultiScatter(blade[, 2], blade[, 3],
                      x.name="observed blade length", 
                      y.name="estimated blade length"),
         file="graphs/scatter.test.blade.original.vs.estimated.svg",
         width=12, height=4)
GraphSVG(MultiScatter(fruit[, 2], fruit[, 3],
                      x.name="observed fruit length", 
                      y.name="estimated fruit length"),
         file="graphs/scatter.test.fruit.original.vs.estimated.svg",
         width=12, height=4)


