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
# none


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
# Plot the result of genus-mean gap-filling against the results of BHPMF
# gap-filling.

scatter.stem.height <- function() {
  Scatterplot(x    = traits.mean.subset$stem.height, 
              y    = traits.BHPMF.subset$stem.height,
              xlab = "Genus Mean estimated stem height",
              ylab = "BHPMF estimated stem height")
  # If mean and BHPMF estimated the same value for a species, then that
  # data point should be on a straight line (x=y). Add this line to the plot:
  lines(x=c(-10, 1000), y=c(-10, 1000))
  # Reset graphical parameters to default
  on.exit(par(mar=c(5, 4, 4, 2) + 0.1))
}
GraphPDF(scatter.stem.height(), file="graphs/scatter.stem.height.estimates.pdf",
         width=6, height=4)



# Remember the standard plotting workflow:
# 1. compose a decent plot
# 2. Open a plot device, for example with pdf()
# 3. Run the plot, as a function
# 4. close the plot device

