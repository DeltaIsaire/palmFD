library(magrittr)
library(plyr)

source(file = "functions/base_functions.R")
source(file = "functions/plotting_functions.R")



# Functional diversity indices
# ----------------------------
# A list with all the FD index data.
# Note these dataframes include richness and endemism variables.
fd.indices <- vector("list", length = 0)
fd.names <- c("FRic.all.traits", "FRic.stem.height", "FRic.blade.length",
              "FRic.fruit.length", "FDis.all.traits", "FDis.stem.height",
              "FDis.blade.length", "FDis.fruit.length"
              )
for (index in fd.names) {
  fd.indices [[index]] <-
    read.csv(file = paste0("output/FD_summary_", index, ".csv"), row.names = 1)
}



#####################################################################
cat("Checking linearity of relation between FD and environment...\n")
#####################################################################

MakePlot <- function(fd.names, null.model, predictors) {
# fd.indices: character vector of fd index names
# null.model: character string of null model name
# predictors: dataframe of predictor variables

  par(mfrow = c(length(fd.indices), ncol(predictors)))
  for (index in fd.names) {
    for (predictor in colnames(predictors)) {
    Scatterplot(x = predictors[, predictor, drop = FALSE],
                y = GetFD(fd.indices[index], null.model),
                ylab = paste0(index, " (", null.model, ")"),
                pch = 16,
                col = tdwg3.info[, "realm"][complete.cases(predictors[, predictor])]
                )
    }
  }

  return (0)
}

for (null.model in c("observed", "global.SES", "realm.SES.noMDG", "adf.SES")) {
  GraphPNG(MakePlot(fd.names, null.model, predictors),
           file = paste0("graphs/scatter_predictors_vs_FD_", null.model, ".png"),
           width = ncol(predictors) * 600,
           height = length(fd.names) * 600,
           pointsize = 24
           )
}
# No obvious nonlinearities in any of these.
# No obvious linearities either. Most of these are pointclouds.



