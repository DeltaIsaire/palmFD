# And do it all again with BIO1_SD removed to minimize collinearity
# -----------------------------------------------------------------
env.complete %<>% .[, !colnames(.) %in% "bio1_sd"]
env.complete.realms <- RealmSubset(env.complete)
n <- 1

for (index in fd.names) {
  for (null.model in null.models) {
    # Prepare data for cases full + realm subsets
    index.data <- GetFD(fd.indices[index], null.model)
    index.data.realms <- RealmSubset(index.data)
    response.list <- c(list(full = index.data[, index]), index.data.realms)
    predictors.list <- c(list(full = env.complete), env.complete.realms)

    for (case in names(response.list)) {
      # verbose output
      nmax <- length(fd.names) * length(null.models) * length(names(response.list))
      cat("\t",
          paste0("(", n, " of ", nmax, ")"),
          paste(index, null.model, "for", case, "dataset"),
          "\n"
          )
      n <- n + 1

      # Handle process flow: for FDis we only need the 'observed' null model
#      if (isTRUE(index %in% fdis) & isTRUE(!null.model == "observed")) {
#        cat("\tSKIP - for FDis we use only the 'observed' data\n")
#        next
#      }
      # Handle process flow: skip if output already exists
      header <-
        paste0(output.dir, "multimod_avg_noBIO1_", index, "_", null.model, "_", case, "_")
      if (file.exists(paste0(header, "SAR.csv"))) {
        cat("\tSKIP - output files already exist\n")
        next
      }

      # Fit global models
      response <- response.list[[case]]
      predictors <- predictors.list[[case]]
      global.ols <- FitGlobalOLS(response = response, predictors = predictors)
      global.sar <- suppressMessages(FitGlobalSAR(response = response,
                                                  predictors = predictors,
                                                  tdwg.map = tdwg.map,
                                                  dist.weight = TRUE
                                                  )
                                     )

      # Export GLOBAL environment to cluster
      clusterExport(cluster,
                    varlist = ls(name = .GlobalEnv),
                    envir = .GlobalEnv
                    )

      # Perform multimodel averaging
      cat("\tOLS:\n")
      mod.avg.ols <-
        DredgeBestMod(global.ols, beta = "none", cluster = cluster)
      cat("\tOLS with standardization:\n")
      mod.avg.ols.std <-
        DredgeBestMod(global.ols, beta = "partial.sd", cluster = cluster)
      cat("\tSAR error:\n")
      mod.avg.sar <-
        DredgeBestMod(global.sar, beta = "none", cluster = cluster)

      # Save results summary to disk
      write.csv(mod.avg.ols[["coef.avg"]],
                file = paste0(header, "OLS.csv"),
                eol = "\r\n"
                )
      write.csv(mod.avg.ols.std[["coef.avg"]],
                file = paste0(header, "OLS_std.csv"),
                eol = "\r\n"
                )
      write.csv(mod.avg.sar[["coef.avg"]],
                file = paste0(header, "SAR.csv"),
                eol = "\r\n"
                )
    }
  }
}

avg.FDis.OLS.noT <- ParseDredge(header = paste0(output.dir, "multimod_avg_noBIO1_"),
                                fd.names = fdis,
                                null.models = null.models,
                                cases = cases,
                                model = "OLS"
                                )

avg.FDis.OLS.std.noT <- ParseDredge(header = paste0(output.dir,
                                                    "multimod_avg_noBIO1_"
                                                    ),
                                    fd.names = fdis,
                                    null.models = null.models,
                                    cases = cases,
                                    model = "OLS_std"
                                    )

avg.FDis.SAR.noT <- ParseDredge(header = paste0(output.dir, "multimod_avg_noBIO1_"),
                                fd.names = fdis,
                                null.models = null.models,
                                cases = cases,
                                model = "SAR"
                                )

avg.FRic.OLS.noT <- ParseDredge(header = paste0(output.dir, "multimod_avg_noBIO1_"),
                                fd.names = fric,
                                null.models = null.models,
                                cases = cases,
                                model = "OLS"
                                )

avg.FRic.OLS.std.noT <- ParseDredge(header = paste0(output.dir,
                                                    "multimod_avg_noBIO1_"
                                                    ),
                                    fd.names = fric,
                                    null.models = null.models,
                                    cases = cases,
                                    model = "OLS_std"
                                    )

avg.FRic.SAR.noT <- ParseDredge(header = paste0(output.dir, "multimod_avg_noBIO1_"),
                                fd.names = fric,
                                null.models = null.models,
                                cases = cases,
                                model = "SAR"
                                )

write.csv(adply(avg.FDis.OLS.noT, .margins = c(3, 4, 5), CheckCI),
          file = paste0(output.dir, "00_95CI_FDis_OLS_noT.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )
write.csv(adply(avg.FDis.OLS.std.noT, .margins = c(3, 4, 5), CheckCI),
          file = paste0(output.dir, "00_95CI_FDis_OLS_std_noT.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )
write.csv(adply(avg.FDis.SAR.noT, .margins = c(3, 4, 5), CheckCI),
          file = paste0(output.dir, "00_95CI_FDis_SAR_noT.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )
write.csv(adply(avg.FRic.OLS.noT, .margins = c(3, 4, 5), CheckCI),
          file = paste0(output.dir, "00_95CI_FRic_OLS_noT.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )
write.csv(adply(avg.FRic.OLS.std.noT, .margins = c(3, 4, 5), CheckCI),
          file = paste0(output.dir, "00_95CI_FRic_OLS_std_noT.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )
write.csv(adply(avg.FRic.SAR.noT, .margins = c(3, 4, 5), CheckCI),
          file = paste0(output.dir, "00_95CI_FRic_SAR_noT.csv"),
          eol = "\r\n",
          quote = FALSE,
          row.names = FALSE
          )
