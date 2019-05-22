sarmod <-
  FitGlobalSAR(response = GetFD(fd.indices["FDis.all.traits"], "observed")[[1]],
               predictors = env[, c("soilcount"), drop = FALSE],
               tdwg.map = tdwg.map,
               dist.weight = TRUE,
               double.std = TRUE,
               numeric.only = TRUE
               )

summary(sarmod, Nagelkerke = TRUE)
PseudoRsq(sarmod)

Scatterplot(x = env[, "soilcount"],
            y = GetFD(fd.indices["FDis.all.traits"], "observed")[[1]],
            )

