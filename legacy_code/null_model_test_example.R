################
# Jun's test:
################

x <- rnorm(100, 3, 1)  # A single trait for 100 species
x_rand_m <- rep(NA, 1000)  # Initialize a statistic for 1000 samples

# For 1000 samples, select 100 random trait values, and calculate the mean
for (i in seq_len(1000)) {
  rand <- sample(x, size = length(x), replace = TRUE)
  x_rand_m[i] <- mean(rand)
}

# Compare original trait mean with the sampled means
(mean(x) - mean(x_rand_m)) / sd(x_rand_m)




















#######################
# Sebastiaan's examples
#######################

# Old situation: trait standardization on null community values
# -------------------------------------------------------------
# 1. 1000 species with one trait
traits.observed <- rnorm(n = 1000, mean = 10, sd = 1)

# 2. Sample 10 null communities of 100 species each, for 2 iterations
#    plus an 11th community that always gets the same species
null.species <- vector("list", length = 11)
null.species.one <- 
  lapply(null.species, sample, x = sample.int(1000), size = 100, replace = FALSE)
null.species.one[[11]] <- sample(x = sample.int(20), size = 20, replace = FALSE)
null.species.two <- 
  lapply(null.species, sample, x = sample.int(1000), size = 100, replace = FALSE)
null.species.two[[11]] <- sample(x = sample.int(20), size = 20, replace = FALSE)


# 3. Subset traits_observed to only species in null.species
traits.subset.one <- traits.observed[sort(unique(unlist(null.species.one)))]
traits.subset.two <- traits.observed[sort(unique(unlist(null.species.two)))]
length(traits.subset.one)  # less than 1000: not all species are in here

# 4. Standardize trait values, using the SUBSETTED trait matrix
traits.std.one <- scale(traits.subset.one, center = TRUE, scale = TRUE)
traits.std.two <- scale(traits.subset.two, center = TRUE, scale = TRUE)

# Evaluation: species 1:20 are always here, because they constitute community #11.
isTRUE(all.equal(traits.std.one[1:20], traits.std.two[1:20]))
# FALSE
# Because the scaling is applied to different subsets of the trait data.


# New situation: trait standardization over the full trait matrix
# ---------------------------------------------------------------
# 1. standardize full trait matrix
traits.observed.std <- scale(traits.observed, center = TRUE, scale = TRUE)

# 2. Using the same null communities, subset the trait matrix
traits.subset.one <- traits.observed.std[sort(unique(unlist(null.species.one)))]
traits.subset.two <- traits.observed.std[sort(unique(unlist(null.species.two)))]

# Evaluation
isTRUE(all.equal(traits.subset.one[1:20], traits.subset.two[1:20]))
# TRUE
# As it should be: standardized trait values do not differ depending on the null
# community










################################
# Why are the z-scores not zero?
################################
# Let's see if the problem can be reproduced with artificial data.# 
# We still focus on community #11.

source("functions/faster_FD_function.R")

# Trial 1
# -------
identical(mean(traits.observed.std[1:20]), mean(traits.subset.one[1:20]))
# TRUE, unsurprisingly.

# But FD calculation is not that simple. 

trait.matrix <- matrix(data = traits.observed.std[1:20],
                       ncol = 1,
                       dimnames = list(1:20, "trait")
                       )
pres.abs.matrix <- matrix(data = rep(1, 20), nrow = 1, dimnames = list(11, 1:20))

run.one <- FastFD(trait.matrix, pres.abs.matrix)
run.two <- FastFD(trait.matrix, pres.abs.matrix)

identical(run.one, run.two)
# TRUE

# Hmmm. 
# ADF null model z-scores for single-trait FRic/FDis are non-zero,
# so multiple traits is not the issue.
# In other words: The problem is not limited to convex hull calculation.
# 'joggling' undoubtedly aggravates the problem, but is not the sole cause.
#
# This example used the FastFD function without finding the error, suggesting
# the source of the error is not in the FastFD function.
# Makes sense: different function output is caused by different inputs.

source("functions/weighted_ADF_null_model_functions.R")

# Test the null community sampler:
species.pool <- data.frame(species = 1:20,
                           weights = rep(1/20, 20)
                           )

test <- SampleADF(species.pool, n = 20, trait.matrix)
# This is fine.



# Test the null model function
test <- NullADF(trait.matrix,
                pres.abs.matrix,
               species.pools = list(species.pool),
                process.dir = "output/test/test/",
                iterations = 100,
                mc.cores = 1,
                subset = FALSE,
                verbose = TRUE,
                single.traits = TRUE,
                fast = TRUE
                )




#########################
# Testing on the real data reveals that non-zero stdevs occur do not occur 
# when the data is subsetted to just HAW(aii), but do occur when a second community
# is included.


# Parallel debug testing revealed the PcoA in the FastFD function is where
# the changes originate, because the resulting eigenvalues differ depending
# on the size of the euclidean distance matrix, which depends on how many
# species are included across all communities.
# However, this only explains FRic, not FDis.
# Then again, FDis is calculated also from the euclidean distance matrix, via
# fdisp.
# Now, the function fdisp also calculates eigenvalues and stuff, which again differ
# depending on the size of the euclidean distance matrix.


# Conclusion: it is possible to get z-scores of exactly zero, if and only if you
# completely seperate FD calculation for each community. I.e., pass to the FastFD
# function the data of only one community at a time.
#
# Using the FastFD function that way is not how it (and the FD::dbFD function) was
# designed, but it is possible. 
#
# Doing that might be a good way to ensure that endemism does not adversely affect
# our z-scores.
