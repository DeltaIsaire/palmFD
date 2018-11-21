###################################################
# Palm FD project: Trait Filling Methods Comparison
###################################################
#
# In which the output of different trait gap-filling methods are compared.
#
# Input files:
#   output/palm.traits.genus.mean.csv
#   output/palm.traits.BHPMF.csv
# Generated output files:
#   graphs/scatter.stem.height.estimates.pdf
#   graphs/scatter.blade.length.estimates.pdf
#   graphs/scatter.fruit.length.estimates.pdf


source(file="functions/base.functions.R")
source(file="functions/plotting.functions.R")


# ---------------
# Read input data
# ---------------
cat("Preparing data...", "\n")
traits.mean <- read.csv(file="output/palm.traits.genus.mean.csv")
traits.BHPMF <- read.csv(file="output/palm.traits.BHPMF.csv")

# The different methods had to exclude different species, so first we must
# subset the filled trait matrices to the list of shared species. Otherwise
# a fair and statistically accurate comparison is not possible.
# Extract list of shared species:
shared.species <- CrossCheck(traits.mean$species,
                             traits.BHPMF$ species,
                             presence=TRUE, value=TRUE)
# Subset filled trait matrices:
traits.mean.subset <- traits.mean[CrossCheck(traits.mean$species,
                                             shared.species, presence=TRUE,
                                             value=FALSE), ]
traits.BHPMF.subset <- traits.BHPMF[CrossCheck(traits.BHPMF$species,
                                             shared.species, presence=TRUE,
                                             value=FALSE), ]


# ---------------------------------------
# Scatterplots of gap-filled trait values
# ---------------------------------------
cat("Creating scatterplots...", "\n")
# When plotting, always store default par values:
opar <- par()

# Plot the result of genus-mean gap-filling against the results of BHPMF
# gap-filling, by invoking the power of functions:
scatter.stem.height <- function() {
  MultiScatter(traits.mean.subset$stem.height, traits.BHPMF.subset$stem.height,
               x.name="Genus Mean", y.name="BHPMF", title="Stem height")
}
scatter.blade.length <- function() {
  MultiScatter(traits.mean.subset$blade.length, traits.BHPMF.subset$blade.length,
               x.name="Genus Mean", y.name="BHPMF", title="Blade length")
}
scatter.fruit.length <- function() {
  MultiScatter(traits.mean.subset$fruit.length, traits.BHPMF.subset$fruit.length,
               x.name="Genus Mean", y.name="BHPMF", title="Fruit length")
}

# Save graphs to pdf:
GraphPDF(scatter.stem.height(), file="graphs/scatter.stem.height.estimates.pdf",
         width=12, height=4)
GraphPDF(scatter.blade.length(), file="graphs/scatter.blade.length.estimates.pdf",
         width=12, height=4)
GraphPDF(scatter.fruit.length(), file="graphs/scatter.fruit.length.estimates.pdf",
         width=12, height=4)

# Just in case:
dev.off()


# Remember the standard plotting workflow:
# 1. compose a decent plot
# 2. Open a plot device, for example with pdf()
# 3. Run the plot, as a function
# 4. close the plot device

