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


# ---------------
# Read input data
# ---------------
palm.trait.genus.mean <- read.csv(file="output/palm.traits.genus.mean.csv")
palm.trait.BHPMF <- read.csv(file="output/palm.traits.BHPMF.csv")


# ---------------------------------
# Scatterplots of gap-filled traits
# ---------------------------------
# TODO: implement this.
# Step 1 will be merging data by species, which is not trivial.
MultiCC(list(mean  = palm.trait.genus.mean$species,
             BHPMF = palm.trait.BHPMF$species),
        presence=FALSE, value=FALSE)


# Remember the standard plotting workflow:
# 1. compose a decent plot
# 2. Open a plot device, for example with pdf()
# 3. Run the plot, as a function
# 4. close the plot device

