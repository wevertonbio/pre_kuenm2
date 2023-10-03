# ####Import calibration results####
# cr <- read.csv("Models/Piper_fuligineum/calibration_results_v2.csv")
# dir <- "Models/Piper_fuligineum/"
#
# #Test function
# cand_models <- cr

####Function to select best models####
sel_best_models <- function(cand_models,
                     test_concave = TRUE,
                     omrat = 5,
                     omrat_threshold = 5,
                     allow_tolerance = T,
                     tolerance = 0.01,
                     AIC = "japones",
                     significance = 0.05,
                     verbose = TRUE,
                     delta_aic = 2,
                     save_file = T,
                     output_dir = NULL){

  if(AIC == "japones") {
    AIC <-   "AIC"} else {
      AIC <-  "AIC_Warren"}

  #Which omission rate?
  om_thr <- paste0("Omission_rate_at_", omrat, ".mean")

  #How many models are being filtered?
  if(verbose){
    message("Filtering ", nrow(cand_models), " models")
  }

  #If test concave = TRUE, remove concave curves
  if(test_concave) {
    if(verbose){
      message("Removing ", nrow(subset(cand_models, is_concave == TRUE)),
              " models with concave curves")
    }
    cand_models <- subset(cand_models, is_concave == FALSE)
  }

  #Remove NAs from results
  if(verbose){
    message("Removing ", nrow(cand_models) - nrow(na.omit(cand_models)),
            " models because NAs were found")
  }
  cand_models <- na.omit(cand_models)

  #Subset models with significativa pROC
  if(verbose){
    message("Removing ", nrow(subset(cand_models,
                                     proc_pval.mean >= significance)),
            " models with non-significant values of pROC")
  }
  cand_models <- subset(cand_models, proc_pval.mean <= significance)

  #Subset models by omission rate
  cand_om <- subset(cand_models, cand_models[,om_thr] <= omrat_threshold/100)
  if(verbose){
    message(nrow(cand_om), " models were selected with omission rate below ",
            omrat_threshold, "%")
  }

  if(nrow(cand_om) == 0 & !allow_tolerance) {
  stop("There is no models with values of omission rate below than ",
       omrat_threshold, ".\nTry with allow_tolerance = T")
  }

  #If 0 models were selected and allow tolerance
  if(nrow(cand_om) == 0 & allow_tolerance) {
    min_thr <- min(cand_models[,om_thr])
    cand_om <- subset(cand_models, cand_models[,om_thr] <= min_thr + tolerance)
    if(verbose){
      message("Minimum value of omission rate (",  round(min_thr*100, 1), "%) is above the selected theshold (", (omrat_threshold),"%).\nApplying tolerance and selecting ", nrow(cand_om), " models with omission rate <", round(min_thr*100 + tolerance, 1), "%")
    }
  }

  #Calculate delta AIC
  cand_om$dAIC <- cand_om[, AIC] - min(cand_om[, AIC])
  #Select delta AIC
  cand_final <- subset(cand_om, cand_om$dAIC <= delta_aic)

  if(verbose){
    message("Selecting ", nrow(cand_final), " final model(s) with delta AIC <",
            delta_aic)
  }

  if(save_file == T){
    if(is.null(output_dir)){
      stop("outpur_dir need to be defined")
      }
  write.csv(cand_final, file.path(output_dir, "selected_models.csv"),
              row.names = F)
  if(verbose){
    message("selected_models.csv written in ", output_dir)
  }
    }

  return(cand_final)
  }

# #Test function
# #With minimum omission rate below the selected threshold
# bm <- sel_best_models(cand_models = cr,
#                        test_concave = TRUE,
#                        omrat = 5,
#                        omrat_threshold = 5, #5%
#                        allow_tolerance = T,
#                        tolerance = 0.01,
#                        AIC = "japones",
#                        significance = 0.05,
#                        verbose = TRUE,
#                        save_file = T,
#                        output_dir = dir,
#                        delta_aic = 2)
# #Save best model
# write.csv(bm, )
#
# #With minimum omission rate above the selected threshold, allowing tolerance
# bm2 <- sel_best_models(cand_models = cr,
#                       test_concave = TRUE,
#                       omrat = 5,
#                       omrat_threshold = 1, #1%
#                       allow_tolerance = T,
#                       tolerance = 0.01,
#                       AIC = "japones",
#                       significance = 0.05,
#                       verbose = TRUE,
#                       delta_aic = 2,
#                       save_file = F,
#                       output_dir = NUL)
#
# #With minimum omission rate above the selected threshold, allowing tolerance
# bm3 <- sel_best_models(cand_models = cr,
#                       test_concave = TRUE,
#                       omrat = 5,
#                       omrat_threshold = 1, #1%
#                       allow_tolerance = F,
#                       tolerance = 0.01,
#                       AIC = "japones",
#                       significance = 0.05,
#                       verbose = TRUE,
#                       delta_aic = 2,
#                       save_file = F,
#                       output_dir = NUL)
#
#
#
#
# ####Function to fit best models####
# fit_best <- function(data, #Data in SWD format
#                      pr_bg,
#                      spat_var,
#                      var_categorical,
#                      best_model,
#                      replicates = TRUE,
#                      summary_method = "median",
#                      n_replicates = 10,
#                      rep_type = c(method = "kfold", n = 4),
#                      parallel = TRUE,
#                      ncores = 1,
#                      progress_bar = TRUE,
#                      writeFiles = FALSE, #Write files?
#                      out_dir = NULL, #Name of the folder to write final models
#                      parallelType = "doSNOW",
#                      verbose = TRUE,
#                      to_export = c("aic_glmnetmx", "aic_maxnet", "eval_stats",
#                                    "get_formulas_maxnet",
#                                    "glmnet_mx", "kfold_part",
#                                    "maxnet.default.regularization",
#                                    "omrat_maxnet",
#                                    "predict.glmnet_mx", "empty_replicates",
#                                    "empty_summary",
#                                    "hinge", "hingeval", "thresholds",
#                                    "thresholdval", "categorical",
#                                    "categoricalval")) {
#
#   #If writefiles = TRUE, create directory
#   if(writeFiles){
#     if(!file.exists(out_dir))
#       dir.create(out_dir)
#   }
#
#   #Convert categorical variables to factor, if necessary
#   if(!is.null(var_categorical)) {
#     if(!is.factor(data[, var_categorical]))
#       data[, var_categorical] <- as.factor(data[, var_categorical])
#   }
#
#   #Remove NA from data
#   n_before <- nrow(data)
#   data <- na.omit(data)
#   n_after <- n_before - nrow(data)
#   if(verbose & n_after > 0){
#     message(n_after, " rows were excluded from database because NAs were found")
#   }
#
#   #Extracts IDs from models
#   m_ids <- bm$ID
#
#   #Fit models with no replicates
#   if(!replicates) {
#     best_maxnet <- lapply(1:nrow(best_model), function(i){
#       bm_i <- bm[i,]
#       best_formula <- bm_i$Formulas
#       best_regm <- bm_i$regm
#       m_i <- glmnet_mx(p = data[,pr_bg], data = data,
#                        f = as.formula(best_formula),
#                        regmult = best_regm, calculate_AIC = F)
#       return(m_i)
#     })
#     #Get predictions to current
#     best_pred <- rast(lapply(seq_along(best_maxnet), function(i){
#       bm_i <- best_maxnet[[i]]
#       best_maxnet_i <- terra::predict(spat_var, bm_i, na.rm = TRUE,
#                                       type = "cloglog")
#       return(best_maxnet_i)
#     }))
#     names(best_pred) <- paste0("Model_", m_ids)
#
#     # #Get consensus of differente models
#     # m_cons <- terra::app(best_pred, summary_method)
#   } #End of !replicates
#
#   #Prepare data (index) to replicates
#
#   if(replicates) {
#     if("kfold" %in% rep_type) {
#       n_fold <- as.numeric(rep_type[["n"]])
#       #Split in n kfolds
#       data <- kfold_part(data, n_folds = n_fold, pr_bg = pr_bg)
#       k_fold <- sort(unique(data[, "folds"]))
#       folds <- "folds"
#
#       best_pred <- lapply(1:nrow(best_model), function(i){
#         bm_i <- bm[i,]
#         best_formula <- bm_i$Formulas
#         best_regm <- bm_i$regm
#
#         pred_rep <- lapply(1:n_replicates, function(z){
#
#
#
#         pred_r <- rast(lapply(1:length(k_fold), function(x) {
#           data_x <- data[data[, folds] != x, ] #Select i k-fold
#           #Run model
#           mod_x <- glmnet_mx(p = data_x[,pr_bg], data = data_x,
#                            f = as.formula(best_formula),
#                            regmult = best_regm,
#                            calculate_AIC = FALSE)
#           #Predict
#           pred_x <- terra::predict(spat_var, mod_x, type = "cloglog",
#                                    na.rm = TRUE)
#           return(pred_x)
#           }))
#         pred_summary <- terra::app(pred_r, summary_method)
#         names(pred_summary) <- paste0("Model_", bm_i$ID)
#         return(pred_summary)
#         })
#       names(best_pred) <- paste0("Model_", m_ids)
#       }) #End of looping in replicates
#     } #End of kfold
#
#
#   } #End of replicates
#
#
#
#
#
#   } #End of function
#
