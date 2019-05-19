require(ggplot2)
require(spdep)
require(relaimpo)
require(MuMIn)
require(care) # to calculate CAR scores

summarizeRelImportance <- function(x){
  # Plots relative
  # x = "averaging" class from MuMIn package
  summaryStats <- data.frame(coefficient = colnames(x$coefficients),
                             fullAvgCoef = summary(x)$coefmat.full[,1],
                             condAvgCoef = summary(x)$coefmat.subset[,1],
                             fullAvgSE = summary(x)$coefmat.full[,2],
                             condAvgSE = summary(x)$coefmat.subset[,2] )
  
  relimportStats <- data.frame(importance = x$importance, coefficient = names(x$importance))
  fullconfintStats <- data.frame(confint(x, full = TRUE))
  names(fullconfintStats) <- c("fulllower2.5CI", "fullupper97.5")
  subsetconfintStats <- data.frame(confint(x, full = FALSE))
  names(subsetconfintStats) <- c("subsetlower2.5CI", "subsetupper97.5")
  confintStats <- cbind(fullconfintStats,subsetconfintStats)
  summaryStats <- cbind(summaryStats, confintStats)
  summaryStats <- merge(summaryStats, relimportStats, by = "coefficient", all = T)
  rownames(summaryStats) <- summaryStats$coefficient
  return(summaryStats)
}

plotRelImportance <- function(x){
  # Plots basic diagnostic output of the summarizeRelImportance function
  ggplot(data = x) + geom_point(aes(y = fullAvgCoef, x= coefficient, size = importance)) + geom_segment(aes(y = lower2.5, yend = upper97.5, x = coefficient, xend = coefficient)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
}

dredge2 <- function(full_mod, data, sw){
  # Fit all possible combinations of predictor variables
  #
  # Arguments:
  #     full_mod, model object with all predictor variables fit
  #     data, data with all predictor variables. Should be the same data.frame used in the model object
  #     sw, spatial weights object
  # 
  # Returns:
  #     list, containing standardized (Cade 2015) and non-standardized model coefficients, coefficient standard errors, AICc, delta AIC and model weights
  
  # Testing
  # full_mod <- knear1_nsw_mod; sw = knear1_sw; data = tdwg_final_glob
  # full_mod <- glob_curr_medBS_mod; data = tdwg_final_glob
  # Extract terms from model object
  data <- full_mod$model
  
  if(class(full_mod) == "lm"){
    full_pred_terms <- attr(full_mod$terms, "term.labels")
  }
  
  lhs_term <- as.character(full_mod$call$formula)[2]
  nterms <- length(full_pred_terms)
  
  # Create all possible combinations of predictor variables
  # Remove first row as no predictor variables are included
  mod_list <- expand.grid(rep(list(0:1), nterms))[-1,] 
  names(mod_list) <- full_pred_terms
  
  # Calculate standard deviations of predictors
  pred_sd <- sapply(full_pred_terms, FUN = function(x){ sd(full_mod$model[[x]]) } )
  
  # Create empty lists to hold results
  mod_obj <- list()
  vif_list <- list()
  mod_AICc <- vector()
  vif_list <- list()
  sd_list <- list()
  pSD_list <- list()
  stdCoef_list <- list()
  coef_list <- list()
  SE_list <- list()
  stdSE_list <- list()

  for(i in 1:nrow(mod_list)){
    rhs_term <- names(mod_list)[mod_list[i,] == 1]
    mod_formula <- as.formula( paste(lhs_term, "~",
                                     paste(rhs_term, collapse = " + ")) )
    if(class(full_mod) == "lm"){
      mod_obj[[i]] <- lm(mod_formula, data = data)  
    }
    nobs <- nobs(mod_obj[[i]])
    npred <- length(rhs_term)
    npara <- length(coef(mod_obj[[i]]))+1# includes intercept and standard dev. in AIC calculations
    
    # Calculate model AICc
    mod_AICc[i] <- AIC(mod_obj[[i]]) + ((2*npara ^2 + 2*npara) / (nobs-npara-1))
    
    sd_list[[i]] <- pred_sd[rhs_term] # sd of the variable
    coef_list[[i]] <- coef( mod_obj[[i]] )[rhs_term]
    
    # Extract standard errors
    if(class(full_mod) == "lm"){
      SE_list[[i]] <- summary(mod_obj[[i]])$coefficients[,2][rhs_term]  
    }
    if(class(full_mod) == "sarlm"){
      SE_list[[i]] <- summary(mod_obj[[i]])$Coef[,2][rhs_term]
    }
    
    if(length(rhs_term) > 1){
      vif_list[[i]] <- vif(mod_obj[[i]])
      pSD_list[[i]] <- sd_list[[i]] * vif_list[[i]]^-0.5 * ((nobs - 1) / (nobs - npred))^0.5
      stdCoef_list[[i]] <- coef( mod_obj[[i]] )[rhs_term] * pSD_list[[i]]
      stdSE_list[[i]] <- (SE_list[[i]]^2 * pSD_list[[i]]^2)^0.5
    } else {
      vif_list[[i]] <- setNames(1, rhs_term)
      pSD_list[[i]] <- NA # Partial SD undefined for single variable models
      stdCoef_list[[i]] <- coef( mod_obj[[i]] )[rhs_term] * sd_list[[i]]
      stdSE_list[[i]] <- (SE_list[[i]]^2 * sd_list[[i]]^2)^0.5
    }
  }
  mod_dAICc <- mod_AICc - min(mod_AICc)
  mod_wt <- exp(-0.5 * mod_dAICc) / sum(exp(-0.5 * mod_dAICc))
  list("Models" = mod_list,
              "Model obj" = mod_obj,
              "AICc" = mod_AICc,
              "dAICc" = mod_dAICc,
              "Weights" = mod_wt,
              "SD" = sd_list,
              "VIF" = vif_list,
              "pSD" = pSD_list,
              "Coef" = coef_list,
              "stdCoef" = stdCoef_list, 
              "SE" = SE_list,
              "stdSE" = stdSE_list)
}


model.avg2 <- function(model_list){
  # Calculates standardized model averaged coefficients
  # Reference: Cade, B.S. (2015) Model averaging and muddled multimodel inferences. Ecology, 96 (9), 2370 - 2382.
  # model_list <- full_moddr
  # Note htatthat the new model.avg does have this functionality under the beta argument
  
  predVars <- names(model_list$Models)
  
  avgStdCoef <- vector()
  avgUnstdCoef <- vector()
  avgStdSE <- vector()
  avgUnstdSE <- vector()
  var_impt <- vector()
  for(i in 1:length(predVars)){
    
    modInd <- which(model_list$Models[predVars[i]] == 1) # models with target pred. variable
    weights <- model_list$Weights[modInd]
    weight_sum <- sum(weights)
    
    var_impt[i] <- weight_sum
    
    StdCoefs <- unlist(lapply(model_list$stdCoef[modInd], function(x) { x[predVars[i]]} ))
    avgStdCoef[i] <- sum(StdCoefs * (weights/weight_sum)) 
    
    UnstdCoefs <- unlist(lapply(model_list$Coef[modInd], function(x) { x[predVars[i]]} ))
    avgUnstdCoef[i] <- sum(UnstdCoefs * (weights/weight_sum)) 
    
    # See Burnham & Anderson 2002, p180 for SE calculation
    StdSEs <- unlist(lapply(model_list$stdSE[modInd], function(x) { x[predVars[i]]} ))
    avgStdSE[i] <- sum(( (StdCoefs - avgStdCoef[i])^2 + (StdSEs)^2)^0.5 * (weights / weight_sum))
    
    UnstdSEs <- unlist(lapply(model_list$SE[modInd], function(x) { x[predVars[i]]} ))
    avgUnstdSE[i] <- sum(( (UnstdCoefs - avgUnstdCoef[i])^2 + (UnstdSEs)^2)^0.5 * (weights / weight_sum))
    
  }
  
  data.frame("Variable" = predVars,
             "avgStdCoef" = avgStdCoef,
             "avgStdSE" = avgStdSE,
             "avgStdCIlower" = avgStdCoef + qnorm(p = 0.025) * avgStdSE,
             "avgStdCIupper" = avgStdCoef + qnorm(p = 0.975) * avgStdSE,
             "avgUnstdCoef" = avgUnstdCoef,
             "avgUnstdSE" = avgUnstdSE,
             "avgUnstdCIlower" = avgUnstdCoef + qnorm(p = 0.025) * avgUnstdSE,
             "avgUnstdCIupper" = avgUnstdCoef + qnorm(p = 0.925) * avgUnstdSE,
             "var_impt" = var_impt)
}

# gpa_data <- data.frame(gpa = c(1.97, 2.74, 2.19, 2.60, 2.98, 1.65, 1.89, 2.38, 2.66, 1.96, 3.14, 1.96, 2.20, 3.90, 2.02, 3.61, 3.07, 2.63, 3.11, 3.20), 
#                        sat_math = c(321, 718, 358, 403, 640, 237, 270, 418, 443, 359, 669, 409, 582, 750, 451, 645, 791, 521, 594, 653),
#                        sat_verb = c(247, 436, 578, 447, 563, 342, 472, 356, 327, 385, 664, 518, 364, 632, 435, 704, 341, 483, 665, 606),
#                        hs_math = c(2.3, 3.8, 2.98, 3.58, 3.38, 1.48, 1.67, 3.73, 3.09, 1.54, 3.21, 2.77, 1.47, 3.14, 1.54, 3.50, 3.20, 3.59, 3.42, 3.69),
#                        hs_eng = c(2.63, 3.57, 2.57, 2.21, 3.48, 2.14, 2.64, 2.52, 3.2, 3.46, 3.37, 2.60, 2.90, 3.49, 3.20, 3.74, 2.93, 3.32, 2.70, 3.52))
# full_mod <- lm(gpa ~ scale(sat_math) + scale(sat_verb) + scale(hs_math) + scale(hs_eng), data = gpa_data, na.action = "na.fail")
# vif(full_mod)
# full_moddr <- dredge2(full_mod)
# full_modavg <- model.avg2(full_moddr)
# 
# full_moddr2 <- dredge(full_mod)
# full_modavg2 <- model.avg(full_moddr2)

# summary(full_modavg2)
# confint(full_modavg2)
# 0.4241378 + 0.3607670 * 1.96
# 0.4241378 + 0.3607670 * -1.96
# qnorm(p = c(0.975, 0.025))


computeModelAvg <- function(x, ...){
  # Takes a model object and performs model averaging
  # Calculates relative variable importance in two ways:
  # 1) Sum of akaike weights of models that include the variable (Burnham & Anderson 2002)
  # 2) Decomposes variance explained by each predictor variable by calculating marginal
  #   correlations adjusted for correlation among explanatory variable(Zuber & Strimmer 2011)
  # 
  # Args:
  #         x: "lm" or "errorsarlm" object
  #       ...: arguments passed to dredge
  #
  # Returns:
  #         data.frame containing coefficients, variance explained and total R2
  
  # Perform model averaging
  moddr <- dredge(x, ...)
  df <- summarizeRelImportance(model.avg(moddr))
  
  # Calculate variable importance sensu Zuber & Strimmer 2011
  if(class(x) == "lm"){
    modvarimpo <- calc.relimp(x, type = "car")
    modvarimpodf <- data.frame(varexp = as.vector(modvarimpo@car),
                               totalR2 = modvarimpo@R2,
                               coefficient = names(modvarimpo@car))
  }
  if(class(x) == "sarlm"){
    mod_spatialonly <- update(x, ~1) # fit a spatial model with only the intercept 
    # mod_spatialonly$lambda <- x$lambda # residuals don't appear to be updated automatically like in lm model class
    # cor(fitted(mod_spatialonly), tdwg_final_glob$logMax95FS_scl, method = "pearson")^2 # 0.453
    
    # Calculate CAR scores, remove first column which is intercept
    mod_car <- carscore(Ytrain = residuals(mod_spatialonly), Xtrain = x$X[,-1], lambda = 0)
    modvarimpo <- mod_car^2
    #cor(print.sarlm.pred(predict.sarlm(x))$trend,
    #    tdwg_final_glob$logMax95FS_scl, method = "pearson")^2 # very different here, but this is for the full model. Acutall much larger in the full model for some reason.
    totalR2 <- cor(fitted(x), x$y, method = "pearson")^2
    # Very similar results when using the predict function
    # cor(print.sarlm.pred(predict.sarlm(x))$fit,
    #     x$y, method = "pearson")^2
    modvarimpodf <- data.frame(varexp = as.vector(modvarimpo),
                               totalR2 = totalR2,
                               coefficient = names(modvarimpo))
  }
  return(merge(df, modvarimpodf, by = "coefficient", all = T))
}

computeModelAvg2 <- function(x){
  # Takes a full model and performs model averaging but calculates adjusted coefficients (following Cade 2015)
  moddr <- dredge2(x)
  modvarimpo <- calc.relimp(x, type = "car")
  modvarimpodf <- data.frame(varexp = as.vector(modvarimpo@car),
                             totalR2 = modvarimpo@R2,
                             Variable = names(modvarimpo@car))
  df <- model.avg2(moddr)
  return(merge(df, modvarimpodf, by = "Variable", all = T))
}

roundNumbers <- function(df, digits = 3){
  # For a given dataframe, round all numeric columns to 3 significant digits
  # Arguments:
  #     df = data.frame
  #     digits = number of digits for rounding
  # Returns:
  #     data.frame
  temp <- sapply(df, FUN = is.numeric)
  cols <- names(temp)[temp]
  for(i in 1:length(cols)){
    df[[cols[i]]] <- round(df[[cols[i]]], digits = digits)
  }
  return(df)
}