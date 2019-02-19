# ------------------------------------------------
# Compare the list of palm species of each dataset
# ------------------------------------------------
# 2. cross-check the three lists of species
missing.table <- MultiCC(species, presence=FALSE, value=FALSE)
missing.list <- MultiCC(species, presence=FALSE, value=TRUE)
# palm.dist and trait.data are the same, but palm.tree differs from the others.

# Test to check that this is true:
new.species <- list(palm.dist = present.list$palm.dist$palm.tree,
                    palm.tree = present.list$palm.tree$palm.dist)
MultiCC(new.species, presence=FALSE, value=FALSE)
# Yup, nothing missing.

# ---------------------------------
# Gap-fill trait matrix using BHPMF
# Initial tests
# ---------------------------------
# To speed up code development or bug-testing you can run BHPMF with a test
# dataset:
test.trait.matrix <- trait.matrix[1:100, ]
test.hierarchy.matrix <- hierarchy.matrix[1:100, ]

# It is good to test how BHPMF works with different input trait matrices,
# Especially because the function doesn't like it when a species has NA
# for all provided traits.

# standard three trait matrix, including all-3 NAs, with substituted zero values;
# and associated hierarchy matrix:
trait.matrix <- as.matrix(palm.traits[, c("stem.height", "blade.length",
                                          "fruit.length")
                                      ]
                          )
rownames(trait.matrix) <- sub(pattern=" ", replacement="_", x=trait.data$SpecName)
trait.matrix[, "stem.height"][which(trait.matrix[, "stem.height"] == 0)] <- 0.0001

hierarchy.matrix <- as.matrix(data.frame(species   = sub(pattern=" ",
                                                         replacement="_",
                                                         x=trait.data$SpecName),
                                         genus     = trait.data$accGenus,
                                         tribe     = trait.data$PalmTribe,
                                         subfamily = trait.data$PalmSubfamily)
                              )

# Optional subsetting to cases without all NAs:
# First extract these cases to a seperate dataframe:
to.remove <- which(rowSums(is.na(trait.matrix)) == ncol(trait.matrix))
BHPMF.missing <- palm.traits[to.remove, ]
# Then remove them from the matrices:
trait.matrix <- trait.matrix[-to.remove, ]
hierarchy.matrix <- hierarchy.matrix[-to.remove, ]
rm(to.remove)
# Result: BHPMF works.

# Optional: add a 4th dummy trait with the same value (1) for all species.
# This could be another way to deal with all NA species.
trait.matrix <- as.matrix(data.frame(trait.matrix,
                                     dummy=rep(1, times=length(trait.matrix[, 1]))
                          )          )
# Result: BHPMF works, even including species with NA for the three real traits.
# TODO: However, we have to check the results for accuracy.

# Optional: remove all but the first column from the trait matrix:
trait.matrix <- as.matrix(trait.matrix[, "stem.height"])
# Result: gives the same error "there are observations with all features missing"


# BHPMF GapFilling function call, with subset for quick processing
# and verbose=TRUE.
# Remember to always clear the preprocessing directory before use:
unlink("output/BHPMF_preprocessing_test", recursive=TRUE)
dir.create("output/BHPMF_preprocessing_test")
# Function disabled to make the script work
#GapFilling(X=trait.matrix, hierarchy.info=hierarchy.matrix,
#           prediction.level=4, used.num.hierarchy.levels=3,
#           mean.gap.filled.output.path="output/TEST.BHPMF.mean.txt",
#           std.gap.filled.output.path="output/TEST.BHPMF.std.txt",
#           tmp.dir="output/BHPMF_preprocessing_test", 
#           rmse.plot.test.data=FALSE, verbose=TRUE)


# Find out how many species are not shared between filled trait matrices:
MultiCC(list(mean  = traits.mean$species,
             BHPMF = traits.BHPMF$species),
        presence=FALSE, value=FALSE)



CrossCheck <- function(x, y, presence = TRUE, value = TRUE) {
# For each value in x, check if it is present or absent in y, then return
# the subset of x values that were present/absent, or their indices.
# TODO: try instead Reduce(intersect, (list()))
# 
# Args:
#   x: Vector whose values should be cross-referenced with the values in y
#   y: Vector against which values in x are compared
#   presence: Logical indicating whether to check for presence (TRUE) or
#             absence (FALSE) of the values of x in y. Default is TRUE.
#   value: Logical indicating whether to return the values of x (TRUE) which
#          are present/absent in y, or their index in x (FALSE). Default is TRUE.
#
# Returns:
#  Vector of the same type as x containing the subset of x present/absent in y,
#  or an integer vector with the indices of x whose values are present/absent
#  in y.
  if (!IsOneDimensional(x)) {
    stop("argument x is not one-dimensional")
  }
  if (!IsOneDimensional(y)) {
    stop("argument y is not one-dimensional")
  }
  if (is.factor(x)) {
    x %<>% as.character()
  }
  if (is.factor(y)) {
    y %<>% as.character()
  }
  checklist <- lapply(x, function(x) { which(x == y) } )
  if (isTRUE(presence)) {
    present <- lapply(checklist, function(x) { length(which(!is.na(x))) } )
    if (isTRUE(value)) {
      return (x[which(present >= 1)])
    } else {
      return (which(present >= 1))
    }
  } else {
    absent <- lapply(checklist, function(x) { which(length(x) == 0) } )
    if (isTRUE(value)) {
      return (x[which(absent >= 1)])
    } else {
      return (which(absent >= 1))
    }
  }
}

MultiCC <- function(x, presence = TRUE, value = TRUE) {
# Wrapper function for CrossCheck. Takes a list and cross-references all
# elements of the list with each other.
#
# Args:
#   x: List with elements to cross-reference.
#   presence: Logical indicating whether to check for presence (TRUE) or
#             absence (FALSE) of the values being cross-referenced.
#             Default is TRUE.
#   value: Logical (see CrossCheck)
#
# Returns:
#   When value = TRUE, returns a list of the values present/absent in each
#   pairwise comparison. When value = FALSE, returns a matrix listing the
#   number of values missing in each pairwise comparison.
  if (IsOneDimensional(x)) {
    stop("Argument x is one-dimensional")
  }
  if (isTRUE(value)) {
    list <- plyr::llply(x,
                        function(a) {
                         plyr::llply(x, 
                                     CrossCheck, 
                                     y = a, 
                                     presence = presence,
                                     value=TRUE
                                     )
                        }
                        )
    return (list)
  } else {
    indices <- plyr::llply(x,
                           function(a) { 
                             plyr::llply(x,
                                         CrossCheck, 
                                         y = a,
                                         presence = presence,
                                         value = FALSE
                                         )
                           }
                           )
    lengths <- plyr::laply(simplify2array(indices), length)
    matrix <- matrix(data = lengths,
                     nrow = length(x),
                     ncol = length(x),
                     byrow = FALSE,
                     dimnames=list(x = names(x),
                                   y = names(x)
                                   )
                     )
    return (matrix)
  }
}


# First, remove genera with only 1 member. A gap in the trait data of such
# genera means estimating a genus mean is not possible.
genera <-
  ddply(complete.traits, "genus", nrow) %>%
  .[!.[, 2] <= 1, ] %>%
  .[, 1] %>%
  as.character()
complete.traits <-
  CrossCheck(complete.traits$genus,
             genera,
             presence = TRUE,
             value = FALSE,
             unique = FALSE
             ) %>%
  complete.traits[., ]




# In other words, a combined plot with 6 frames in 3
# columns, where each column is a trait and each row a method.
SixScatter <- function() {
  par(mfrow = c(2, 3))
  for (df in seq_along(all.estimates)) {
    Scatterplot(x = all.estimates[[df]][, "original"],
                y = all.estimates[[df]][, "std.BHPMF"],
                xlab = paste("Original",
                             trait.names[df]
                             ),
                ylab = paste("Est.",
                             trait.names[df],
                             "(standard BHPMF)"
                             )
                )
    # Add 1:1 line for reference:
    lines(x = c(par("usr")[1], par("usr")[2]), 
          y = c(par("usr")[1], par("usr")[2])
          )
  }
  for (df in seq_along(all.estimates)) {
    Scatterplot(x = all.estimates[[df]][, "original"],
                y = all.estimates[[df]][, "growthform.BHPMF"],
                xlab = paste("Original",
                             trait.names[df]
                             ),
                ylab = paste("Est.",
                             trait.names[df],
                             "(growthform BHPMF)"
                             )
                )
    # Add 1:1 line for reference:
    lines(x = c(par("usr")[1], par("usr")[2]), 
          y = c(par("usr")[1], par("usr")[2])
          )
  }
}









# TDWG3 units HAW, NFK, NWC and SEY need to be excluded. These areas have 100%
# endemic palm assemblages, for which a weighted ADF cannot be defined.
#
# That's Hawaii, Norfolk Island (Oceania), New Caledonia and Seychelles.
# All of them islands, and HAW + NWC have fairly high palm species richness.
# All HAW species are in genus Pritchardia, which is unique to HAW.
# NFK has 5 species in 4 genera, all genera unique to NFK.
# NWC has 39 species in 10 genera, not all genera unique to NWC.
# SEY has 6 species in 6 genera, most genera are unique to SEY.
#
# Subset the datasets:
pres.abs.matrix %<>% .[!row.names(.) %in% c("HAW", "NFK", "NWC", "SEY", "MAU", "MDG", "REU"), ]
pres.abs.matrix %<>% .[, colSums(.) > 0]
traits.mean %<>% .[rownames(.) %in% colnames(pres.abs.matrix), ]
traits.gapfilled <-
  llply(traits.gapfilled,
        function(x) {
          x[rownames(x) %in% colnames(pres.abs.matrix), ]
        }
        )




  # Community weighted mean trait values
  # Not weighted in our case, so just the mean
  # KEEP IN MIND THESE ARE LOG10-TRANSFORMED VALUES
  community.means <- output.combined[[1]]$CWM
  write.csv(community.means,
            file = paste0(fd.dir, "test_community_trait_means_", id, ".csv"),
            eol = "\r\n",
            row.names = TRUE
            )

