
#############################################
# Functional Richness of palms in tdwg3 units
#############################################

# -------------------------------------------
# Correlating FRic (range) with palm richness
# -------------------------------------------
cat("Correlating FRic with palm richness...\n")

fric.list <- 
  list(all.filled      = fd.filled[, c("FRic", "richness")],
       all.unfilled    = fd.unfilled[, c("FRic", "richness")],
       height.filled   = fd.single.filled[, c("height.FRic", "richness")],
       blade.filled    = fd.single.filled[, c("blade.FRic", "richness")],
       fruit.filled    = fd.single.filled[, c("fruit.FRic", "richness")],
       height.unfilled = fd.single.unfilled[, c("height.FRic", "richness")],
       blade.unfilled  = fd.single.unfilled[, c("blade.FRic", "richness")],
       fruit.unfilled  = fd.single.unfilled[, c("fruit.FRic", "richness")]
       )
fric.mods <- llply(fric.list,
                   function(x) {
                     lm(paste(names(x)[1], "~", names(x)[2]),
                        data = x
                        )
                   }
                   )
fric.mods.summary <-
  data.frame(FRic.trait = names(fric.list),
             p.value    = laply(fric.mods,
                                function(x) {
                                  summary(x)$coefficients[8]
                                }
                                ),
             R.squared  = laply(fric.mods,
                                function(x) {
                                  summary(x)$r.squared
                                }
                                )
             )
# All significant, but much lower R-squared for the single-trait fric compared to
# all-trait fric.
# Not much of a difference between filled and unfilled.

# graphs of correlations
# ----------------------
RichRichPlot <- function() {
  par(mfrow = c(1, 2))
  # Gap-filled data:
  Scatterplot(x = fd.filled$richness,
              y = fd.filled$FRic,
              xlab = "Species richness",
              ylab = "Functional richness",
              title = "Gap-filled dataset"
              )
  lines(x = fd.filled$richness, y = predict(fric.mods[[1]]))
  # Unfilled data:
  Scatterplot(x = fd.unfilled$richness,
              y = fd.unfilled$FRic,
              xlab = "Species richness",
              ylab = "Functional richness",
              title = "Unfilled dataset"
              )
  lines(x = fd.unfilled$richness, y = predict(fric.mods[[2]]))
}
GraphSVG(RichRichPlot(),
         file = "graphs/test/test_scatter_fric_vs_richness.svg",
         width = 12,
         height = 4
         )
# WORD OF CAUTION: From these graphs the correlation appears to be curvilinear,
# not linear. There appears to be a slight saturation effect, which makes sense.

# Single traits:
SinglesPlot <- function() {
  par(mfrow = c(2, 3))
  for (i in c(3, 4, 5, 6, 7, 8)) { 
    Scatterplot(x = fric.list[[i]][ ,2],
                y = fric.list[[i]][ ,1],
                xlab = names(fric.list[[i]])[2],
                ylab = names(fric.list[[i]])[1],
                title = names(fric.list)[i]
                )
    lines(x = fric.list[[i]][ ,2], y = predict(fric.mods[[i]]))
  }
}
GraphSVG(SinglesPlot(),
         file = "graphs/test/test_scatter_fric_vs_richness_singles.svg",
         width = 8,
         height = 4
         )
# Those are curvilinear relationships with much stronger saturation!

# -------------------------------------------
# Plotting Functional Richness in tdwg3 units
# -------------------------------------------
cat("Plotting FRic in tdwg3 units...\n")

plot.fric.filled <-
  ggplot(data = spatial.filled) +
         geom_sf(color = "black", aes(fill = FRic)) +
         scale_fill_viridis_c(option = "plasma") +
         ggtitle("Distribution of FRic (gap-filled dataset)")
ggsave(plot = plot.fric.filled,
       filename = "graphs/test/test_distribution_FRic_filled.png",
       width = 8,
       height = 4
       )


###############################################
# Functional Dispersion of palms in tdwg3 units
###############################################

# -----------------------------------
# Correlating FDis with palm richness
# -----------------------------------
cat("Correlating FDis with palm richness...\n")

fdis.list <- 
  list(all.filled      = fd.filled[, c("FDis", "richness")],
       all.unfilled    = fd.unfilled[, c("FDis", "richness")],
       height.filled   = fd.single.filled[, c("height.FDis", "richness")],
       blade.filled    = fd.single.filled[, c("blade.FDis", "richness")],
       fruit.filled    = fd.single.filled[, c("fruit.FDis", "richness")],
       height.unfilled = fd.single.unfilled[, c("height.FDis", "richness")],
       blade.unfilled  = fd.single.unfilled[, c("blade.FDis", "richness")],
       fruit.unfilled  = fd.single.unfilled[, c("fruit.FDis", "richness")]
       )
fdis.mods <- llply(fdis.list,
                   function(x) {
                     lm(paste(names(x)[1], "~", names(x)[2]),
                        data = x
                        )
                   }
                   )
fdis.mods.summary <-
  data.frame(FDis.trait = names(fdis.list),
             p.value    = laply(fdis.mods,
                                function(x) {
                                  summary(x)$coefficients[8]
                                }
                                ),
             R.squared  = laply(fdis.mods,
                                function(x) {
                                  summary(x)$r.squared
                                }
                                )
             )
# No significant correlations, except for blade.unfilled which
# has p-value < 0.01, but R-squared fo only 0.07.
# So functional dispersion is not correlated with species richness.
# Of course FDis isn't correlated with richness by design, so these tests
# are merely a checksum anyway.

# -------------------------------------------
# Plotting Functional Dispersion in tdwg3 units
# -------------------------------------------
cat("Plotting FDis in tdwg3 units...\n")

plot.fdis.filled <-
  ggplot(data = spatial.filled) +
         geom_sf(color = "black", aes(fill = FDis)) +
         scale_fill_viridis_c(option = "plasma") +
         ggtitle("Distribution of FDis (gap-filled dataset)")
ggsave(plot = plot.fdis.filled,
       filename = "graphs/test/test_distribution_FDis_filled.png",
       width = 8,
       height = 4
       )

