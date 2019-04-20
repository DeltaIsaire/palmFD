# - spatial weights matrix with and without distance-weighting
# - all single predictors
# - FRic + FDis
# - null models global, realm, realmnoMDG, adf, and observed
# - global dataset and subset for each three realms
# That is a lot of combinations!
# Fortunately, we can borrow code from the OLS regressions, so developing functions
# for single-predictor SAR models will is relatively easy.


cat("test 1...\n")
model.data <- env.complete
model.data[, "null"] <- runif(n = nrow(model.data))
model.data[, "response"] <- fd.indices[[1]] [, "global.SES"]

result <- SingleSAR(model.data,
                    listw = nb.soi.swmat,
                    response = "response"
                    )


cat("test 2...\n")
responses <- fd.indices[[1]] [, c("global.SES", "realm.SES", "adf.SES")]
# TODO: if you use 'realm.SES.noMDG' you need a special neighbourhood, where MDG
# is removed!

result2 <- MultiSingleSAR(responses = responses,
                          predictors = env.complete,
                          listw = nb.soi.swmat
                          )

cat("test 3...\n")
result3 <- RunMSSAR(name = paste0(output.dir, "SAR_single_test"),
                    responses = responses,
                    predictors = env.complete,
                    listw = nb.soi.swmat
                    )

cat("test 4...\n")
AllSARSingles(fd.indices[fric],
              colname = "global.SES",
              predictors = env.complete,
              tdwg.map = tdwg.map,
              dist.weight = FALSE,
              name.all = paste0(output.dir, "SAR_single_test_FRic_global"),
              standardize = TRUE,
              numeric.only = FALSE
              )

cat("test 5...\n")
result5 <-
  ParseSARSingle(name.all = paste0(output.dir, "SAR_single_test_FRic_global"),
                 cases = c("full", "NewWorld", "OWWest", "OWEast"),
                 )
