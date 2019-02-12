

############################################
# Determine optimal iterations for final run
############################################
cat("Assessing z-score convergence...\n")
# How many iterations is enough to get convergence on a z-score?
#
# We can test this quickly by running the null model n times for a few hand-picked
# communities.
# Which communities? The ones with highest/lowest raw FRic and FDis.
# That's MLY, WAU, SCZ and TCI, respectively.

# Subset data to these tdwg3 units
# --------------------------------
pres.abs.subset <- 
  pres.abs.gapfilled[rownames(pres.abs.gapfilled) %in% c("MLY",
                                                         "WAU",
                                                         "SCZ",
                                                         "TCI"
                                                          ),
                     ]
fd.test <- fd.gapfilled[rownames(fd.gapfilled) %in% rownames(pres.abs.subset),
                        c(1, 5)
                        ]
orphaned.species <-
  which(colSums(pres.abs.subset) == 0) %>%
  colnames(pres.abs.subset)[.]
pres.abs.subset %<>% .[, -which(colSums(.) == 0)]
indices <- CrossCheck(x = rownames(traits.gapfilled),
                      y = orphaned.species,
                      presence = TRUE,
                      value = FALSE
                      )
traits.subset <- traits.gapfilled[-indices, ]
realm.test <- realm.tdwg3[c(1, 3)]

# Run regional null model for n iterations and compute z-values
# -------------------------------------------------------------
# We should test for iterations in (100 * 1:10)
# To save processing time we can cheat: instead of running ten times for ten
# different iterations, we run it once for 1000 iterations, and then
# calculate z-scores based on partial output.
# We do still need say 10 z-scores for each number of iterations, which means
# running 1000 iterations ten times.
# To further save on processing time, exclude single-trait FD and focus on
# all-traits FD.
#
# Only run code if the output doesn't exist
if (!file.exists("output/test/nullmodel_regional_z_convergence_FRic_sd.csv")) {
  test.results <- vector("list", length = 10)

  for (run in 1:10) {
    # init partial results list
    partial.results <- unlist(list(rep(list(fd.test), 10)), recursive = FALSE)
    names(partial.results) <- paste0("z.score.", 100 * 1:10)
    # clean processing dir
    if (dir.exists("output/test/nullmodel_test/")) {
      unlink("output/test/nullmodel_test/", recursive = TRUE)
    }
    # run nullmodel
    regional.test <-
      NullModel(trait.matrix = traits.subset,
                pres.abs.matrix = pres.abs.subset,
                groups = realm.test,
                process.dir = "output/test/nullmodel_test/",
                iterations = 1000,
                mc.cores = num.cores,
                subset = subset,
                verbose = verbose,
                random.groups = TRUE,
                single.traits = FALSE
                )
    # compute z-scores
    for (i in 1:10) {
      # subset nullmodel output
      regional.test.subset <-
        llply(regional.test, function(x) { x[[1]] [1:(100 * i), ] } )
      names(regional.test.subset) <- c("FRic.all.traits", "FDis.all.traits")
      partial.results[[i]] <- NullTransform(fd.test, regional.test.subset)
    }
    test.results[[run]] <- partial.results
  }
  
  # Compute sd of z-values for each number of iterations
  # ----------------------------------------------------
  # Init result list
  mat <- matrix(ncol = 10,
                nrow = nrow(fd.test),
                dimnames = list(row = rownames(fd.test),
                                col = paste0("iterations.", 100 * 1:10)
                                )
                )
  sd.results <- list(FRic = mat, FDis = mat)
  mean.results <- list(FRic = mat, FDis = mat)
  # apply magic
  values <- numeric(length(test.results))
  for (run in seq_along(test.results)) {
    for (index in 1:2) {
      for (iterations in 1:10) {
        for (area in seq_len(nrow(fd.test))) {
          values[i] <- test.results[[run]] [[iterations]] [area, index]
          sd.results[[index]] [area, iterations] <- sd(values)
          mean.results[[index]] [area, iterations] <- mean(values)
        }
      }
    }
  }
  sd.results
  mean.results

  # Save results
  # ------------
  write.csv(sd.results[[1]],
            file = "output/test/nullmodel_regional_z_convergence_FRic_sd.csv",
            eol = "\r\n",
            row.names = TRUE
            )
  write.csv(sd.results[[2]],
            file = "output/test/nullmodel_regional_z_convergence_FDis_sd.csv",
            eol = "\r\n",
            row.names = TRUE
            )
  write.csv(mean.results[[1]],
            file = "output/test/nullmodel_regional_z_convergence_FRic_mean.csv",
            eol = "\r\n",
            row.names = TRUE
            )
  write.csv(mean.results[[2]],
            file = "output/test/nullmodel_regional_z_convergence_FDis_mean.csv",
            eol = "\r\n",
            row.names = TRUE
            )
} else {
  # Read the saved output:
  sd.results <- 
    list(FRic = read.csv("output/test/nullmodel_regional_z_convergence_FRic_sd.csv",
                         header = TRUE,
                         row.names = 1
                         ),
         FDis = read.csv("output/test/nullmodel_regional_z_convergence_FDis_sd.csv",
                         header = TRUE,
                         row.names = 1
                         )
         )
  mean.results <- 
    list(FRic = read.csv("output/test/nullmodel_regional_z_convergence_FRic_mean.csv",
                         header = TRUE,
                         row.names = 1
                         ),
         FDis = read.csv("output/test/nullmodel_regional_z_convergence_FDis_mean.csv",
                         header = TRUE,
                         row.names = 1
                         )
         )
}


# Interpretation of sd.results:
# FRic - WAU is stable after 200 iterations.
#        TCI is stable after 300 iterations.
#        SCZ is stable after 400 iterations.
#        MLY is stable after 100 iterations.
# FDis - WAU is stable after 100 iterations.
#        TCI is stable after 100 iterations.
#        SCZ is stable after 100-200 iterations.
#        MLY is stable after 400 iterations.
# Overall, something like 300 should be quite sufficient.
# In fact the difference between 100 iterations sd or 500 iterations sd is very
# small (less than 10% on average). A few hundred iterations will be absolutely
# fine. 500 or more would be overkill for all practical purposes.
#
# Interpretation of mean.results:
# FRic - WAU is accurate at 200 iterations
#        TCI is accurate at 100 iterations
#        SCZ is accurate at 200 iterations
#        MLY is accurate at 100 iterations
# FDis - WAU is accurate at 100 iterations
#        TCI is accurate at 100-200 iterations
#        SCZ is accurate at 100-200 iterations
#        MLY is accurate at 400 iterations
# Keep in mind these are averages over 10 runs, so being accurate after x iterations
# means the value is based on x * 10 iterations.
# The slowest to converge was FDis for MLY, which has the highest FDis Z-score.
#
# Final conclusion: You need at most 400 iterations. 
# 300 would also be fine in the majority of cases.


#######################################
# Code to visualize z-score convergence
#######################################
cat("Visualizing z-score convergence...\n")
# For a few carefully chosen communities, we can visualize the convergence
# of their z-scores as the number of iterations increases.
# How? By computing the z-score incrementally.
# Which communities? The ones with highest/lowest raw FRic and FDis.
# That's MLY, WAU, SCZ and TCI, respectively.
areas <- sort(c("MLY", "WAU", "SCZ", "TCI"))
#
# For practical purposes, we can focus only on the all-traits FD for the gapfilled
# dataset. The unfilled and single-trait data is assumed to converge in a similar
# fashion.

# Subset data to a form usable by ZScore()
# FRic:
temp <- fd.gapfilled[rownames(fd.gapfilled) %in% areas, 1]
raw.FRic.subset <- structure(temp, names = areas)
null.FRic.subset <-
  regional.gapfilled[[1]] [[1]] [, colnames(regional.gapfilled[[1]] [[1]]) %in%
                                   areas]
# FDis:
temp <- fd.gapfilled[rownames(fd.gapfilled) %in% areas, 5]
raw.FDis.subset <- structure(temp, names = areas)
null.FDis.subset <-
  regional.gapfilled[[2]] [[1]] [, colnames(regional.gapfilled[[2]] [[1]]) %in%
                                   areas]

# Initialize output
# z.scores require at least 2 iterations, so the first one doesn't count
z.scores.FRic <- matrix(ncol = length(areas),
                        nrow = nrow(null.FRic.subset) - 1,
                        dimnames = list(row = seq_len(nrow(null.FRic.subset) - 1),
                                        col = areas
                                        )
                        )
z.scores.FDis <- matrix(ncol = length(areas),
                        nrow = nrow(null.FDis.subset) - 1,
                        dimnames = list(row = seq_len(nrow(null.FDis.subset) - 1),
                                        col = areas
                                        )
                        )

# calculate incremental means
# Instead of NullTransform() we use ZScore() directly.
# z.scores require at least 2 iterations, so the first one doesn't count
for (i in seq_len(nrow(null.FRic.subset))[-1]) {
  z.scores.FRic[(i - 1), ] <-
    ZScore(raw.FRic.subset, null.FRic.subset[(1:i), ,drop = FALSE])
}
for (i in seq_len(nrow(null.FDis.subset))[-1]) {
  z.scores.FDis[(i - 1), ] <-
    ZScore(raw.FDis.subset, null.FDis.subset[(1:i), ,drop = FALSE])
}

# Visualize results
MyPlot <- function(z.scores, index) {
  par(mfrow = c(2, 2))
  for (i in seq_along(areas)) {
    Scatterplot(x = seq_len(nrow(z.scores)),
                y = z.scores[, i],
                xlab = "iterations",
                ylab = "z.score",
                title = paste(index,
                              "z.score convergence in",
                              colnames(z.scores)[i]
                              )
                )
  }
}
GraphSVG(MyPlot(z.scores.FRic, "FRic"),
         file = "graphs/test/test_nullmodel_regional_z_convergence_FRic.svg",
         width = 8,
         height = 8
         )
GraphSVG(MyPlot(z.scores.FDis, "FDis"),
         file = "graphs/test/test_nullmodel_regional_z_convergence_FDis.svg",
         width = 8,
         height = 8
         )




cat("Done.\n")

