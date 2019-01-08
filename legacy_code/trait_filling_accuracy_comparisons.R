
# Scatterplots of standard BHPMF vs dummy BHPMF estimates
# -------------------------------------------------------
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "std.BHPMF"],
                        all.estimates[[df]][, "dummy.BHPMF"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "std.BHPMF"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "dummy.BHPMF"
                                       )
                        ),
             file=paste0("graphs/test_sparse_scatter_",
                         trait.names[df],
                         "_std.BHPMF_vs_dummy.BHPMF.svg"
                         ),
             width = 12,
             height = 4
             )
}

# Scatterplots of genus mean vs standard BHPMF estimates
# ------------------------------------------------------
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "mean"],
                        all.estimates[[df]][, "std.BHPMF"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "mean"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "std.BHPMF"
                                       )
                        ),
             file=paste0("graphs/test_sparse_scatter_",
                         trait.names[df],
                         "_genus_mean_vs_std.BHPMF.svg"
                         ),
             width = 12,
             height = 4
             )
}

# Scatterplots of standard BHPMF vs growthform BHPMF estimates
# ------------------------------------------------------------
for (df in seq_along(all.estimates)) {
  GraphSVG(MultiScatter(all.estimates[[df]][, "std.BHPMF"],
                        all.estimates[[df]][, "growthform.BHPMF"],
                        x.name = paste("Est.",
                                       trait.names[df],
                                       "std.BHPMF"
                                       ),
                        y.name = paste("Est. ",
                                       trait.names[df],
                                       "growthform.BHPMF"
                                       )
                        ),
             file=paste0("graphs/test_sparse_scatter_",
                         trait.names[df],
                         "_std.BHPMF_vs_growthform.BHPMF.svg"
                         ),
             width = 12,
             height = 4
             )
}

# Scatterplots of original vs estimates for BHPMF: standard + growthform
# ----------------------------------------------------------------------
# For each trait, plot the original values vs the estimates from standard BHPMF
# and growthform BHPMF. 
OneScatterGrowthform <- function(index) {
  Scatterplot(x = all.estimates[[index]][, "original"],
              y = all.estimates[[index]][, "std.BHPMF"],
              xlab = paste("Original",
                           trait.names[index]
                           ),
              ylab = paste("Estimated",
                           trait.names[index]
                           ),
              pch = 17
              )
    points(x = all.estimates[[index]][, "original"],
           y = all.estimates[[index]][, "growthform.BHPMF"],
           col = "red",
           pch = 19
           )
  # Add 1:1 line for reference:
  lines(x = c(par("usr")[1], par("usr")[2]), 
        y = c(par("usr")[1], par("usr")[2])
        )
  # legend
  legend("bottomright",
         legend = c("Standard BHPMF", "Growthform BHPMF"),
         col = c("black", "red"),
         pch = c(17, 19),
         bg = "white"
         )
}

for (i in seq_along(all.estimates)) {
  GraphSVG(OneScatterGrowthform(i),
           file = paste0("graphs/test_sparse_scatter_",
                         trait.names[i],
                         "_original_vs_std_and_growthform_BHPMF.svg"
                         ),
           width = 6,
           height = 4
           )
}

# Scatterplots of original vs estimates for BHPMF: standard + binary
# ------------------------------------------------------------------
# For each trait, plot the original values vs the estimates from standard BHPMF
# and growthform BHPMF. 
OneScatterBinary <- function(index) {
  Scatterplot(x = all.estimates[[index]][, "original"],
              y = all.estimates[[index]][, "std.BHPMF"],
              xlab = paste("Original",
                           trait.names[index]
                           ),
              ylab = paste("Estimated",
                           trait.names[index]
                           ),
              pch = 17
              )
    points(x = all.estimates[[index]][, "original"],
           y = all.estimates[[index]][, "binary.BHPMF"],
           col = "red",
           pch = 19
           )
  # Add 1:1 line for reference:
  lines(x = c(par("usr")[1], par("usr")[2]), 
        y = c(par("usr")[1], par("usr")[2])
        )
  # legend
  legend("bottomright",
         legend = c("Standard BHPMF", "Binary BHPMF"),
         col = c("black", "red"),
         pch = c(17, 19),
         bg = "white"
         )
}

for (i in seq_along(all.estimates)) {
  GraphSVG(OneScatterBinary(i),
           file = paste0("graphs/test_sparse_scatter_",
                         trait.names[i],
                         "_original_vs_std_and_binary_BHPMF.svg"
                         ),
           width = 6,
           height = 4
           )
}

cat("Done.\n")

# Combined 2x3 plot of standard BHPMF vs growthform and binary, for all traits
SixScatter <- function() {
  par(mfrow = c(3, 2))
  for (i in seq_along(all.estimates)) {
    OneScatterGrowthform(i)
    OneScatterBinary(i)
  }
}
GraphSVG(SixScatter(),
         file = "graphs/test_sparse_scatter_std_vs_growthform_binary_BHPMF.svg",
         width = 8,
         height = 8
         )

# Extra comparisons of original vs mean vs std.BHPMF
# --------------------------------------------------
MultiHist(data = all.estimates[[1]][, 2:4],
          id = names(all.estimates[[1]])[2:4],
          xlab = "stem height"
          )
# Original is highlyl skewed, estimates less so
MultiHist(data = all.estimates[[2]][, 2:4],
          id = names(all.estimates[[2]])[2:4],
          xlab = "blade length"
          )
# Distributions much more similar
MultiHist(data = all.estimates[[2]][, 2:4],
          id = names(all.estimates[[2]])[2:4],
          xlab = "fruit length"
          )
# Like Blade length, distributions much more similar.

# Same thing for residuals:
all.residuals <- llply(all.estimates,
                       function(df) {
                         data.frame(original    = df$original,
                                    resid.mean  = df$mean - df$original,
                                    resid.BHPMF = df$std.BHPMF - df$original
                                    )
                       }
                       )
MultiHist(data = all.residuals[[1]][, 2:3],
          id = names(all.residuals[[1]])[2:3],
          xlab = trait.names[1]
          )
MultiHist(data = all.residuals[[2]][, 2:3],
          id = names(all.residuals[[2]])[2:3],
          xlab = trait.names[2]
          )
MultiHist(data = all.residuals[[3]][, 2:3],
          id = names(all.residuals[[3]])[2:3],
          xlab = trait.names[3]
          )

# Find those stem height outliers:
palm.traits[which(palm.traits$stem.height > 100), ]
palm.traits[which(palm.traits$genus == "Calamus"), ]
palm.traits[which(palm.traits$genus == "Eremospatha"), ]
# In both cases, the 150m climbing palm resides in a genus filled with much
# shorter companions.
# A tiny speck of hope: most species for which we do not know the height
# are not likely to be >100m.

