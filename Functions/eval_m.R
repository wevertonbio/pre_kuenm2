#Looping through candidate models using cores
eval_m <- function(data, #Data in SWD format
          pr_bg, #Column name with presence (1) or background (0)
          formula_grid, #Grid with formulas
          var_categorical = NULL,
          test_concave = TRUE, #Test concave curves in quadratic models?
          folds = 4, #Columns name with k_folds or vector indicating k_folds
          parallel = TRUE,
          ncores = 1,
          progress_bar = TRUE, #Show progress bar? Only works if parallelType = "doSNOW"
          writeFiles = FALSE, #Write candidate evaluations?
          out_dir = NULL, #Name of the folder to write candidate evaluations
          parallelType = "doSNOW",
          only_summary = TRUE,
          omrat_thr = c(5, 10),
          skip_existing_models = FALSE, #Only works if writeFiles = TRUE
          verbose = TRUE,
          #Check if it's necessary to export in the package
          to_export = c("aic_glmnetmx", "aic_maxnet", "eval_stats",
    "get_formulas_maxnet",
    "glmnet_mx", "kfold_part",
    "maxnet.default.regularization",
    "omrat_maxnet",
    "predict.glmnet_mx", "empty_replicates",
    "empty_summary",
    "hinge", "hingeval", "thresholds",
    "thresholdval", "categorical",
    "categoricalval")){
  #If writefiles = TRUE, create directory
  if(isTRUE(writeFiles)){
    if(!file.exists(out_dir))
      dir.create(out_dir)
  }

  #If skip_existing_models = TRUE, update grid
  if(skip_existing_models & writeFiles) {
    ready_models <- list.files(out_dir, full.names = T)
    ready_models <- do.call("rbind", lapply(seq_along(ready_models), function(i){
      read.csv( ready_models[i])
    }))
    run_models <- setdiff(formula_grid$ID, ready_models$ID)
    if(length(run_models) == 0) {
      stop(paste("All models completed. Check the folder:", out_dir))
    } else { #Update formula grid
      formula_grid <- formula_grid[formula_grid$ID %in% run_models, ]
    }
  }

  #Convert categorical variables to factor, if necessary
  if(!is.null(var_categorical)) {
    if(!is.factor(data[, var_categorical]))
      data[, var_categorical] <- as.factor(data[, var_categorical])
  }

  #Remove NA from data
  n_before <- nrow(data)
  data <- na.omit(data)
  n_after <- n_before - nrow(data)
  if(verbose & n_after > 0){
   message(n_after, " rows were excluded from database because NAs were found")
  }

  #Get k_folds
  if(is.numeric(folds)){
    data <- kfold_part(data, n_folds = folds, pr_bg = pr_bg)
    k_fold <- sort(unique(data[, "folds"]))
    folds <- "folds"
  } else {
    k_fold <- sort(unique(data[, folds]))
  }


  #Make cluster
  cl <- parallel::makeCluster(ncores)

  #If test_concave = TRUE
  if(test_concave){
    if(verbose){
      cat("\n
        Task 1 of 2: checking concave curves in quadratic models\n")
    }
    #Get only quadratic grids with higher regularization multiplier
    q_grids <- formula_grid[grepl("q", formula_grid$Features) &
                              formula_grid$regm == max(formula_grid$regm), ]

    #Set number of iteration
    n <- nrow(q_grids)
    #If n = 0, do not run
    if(n == 0) {
      warning("All quadratic models have been already tested")
    } else {
      #Show progress bar?
      if (isTRUE(progress_bar)) {
        pb <- txtProgressBar(0, n, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n) }

      if (parallelType == "doParallel") {
        doParallel::registerDoParallel(cl)
        opts <- NULL
      }

      if (parallelType == "doSNOW") {
        doSNOW::registerDoSNOW(cl)
        if (isTRUE(progress_bar))
 opts <- list(progress = progress)
        else opts <- NULL
      }
      results_concave <- foreach::foreach(
        x = 1:n, .options.snow = opts,
        .export = to_export) %dopar% {
 #Get grid x
 grid_x <- q_grids[x,] #Get i candidate model
 formula_x <- as.formula(grid_x$Formulas) #Get formula from grid x
 reg_x <- grid_x$regm #Get regularization multiplier from grid x

 #Complete model with AIC
 m_aic <- try(glmnet_mx(p = data[,pr_bg], data = data,
           f = formula_x, regmult = reg_x, calculate_AIC = T),
     silent = TRUE)
 if(any(class(m_aic) == "try-error")) {
      npar <- NA
      AICc <- NA
      is_c <- NA
      mods <- NA
      class(mods) <- "try-error"
      } else {

 #Get number of parameters
 npar <- length(m_aic$betas)

 #Calculate AIC from Warren
 vals <- predict.glmnet_mx(m_aic, data[data[,pr_bg] == 1, ],
       type = "exponential")
 LL <- sum(log(vals + .Machine$double.eps))
 #npar <- length(m_aic$betas)
 AICc <- ((2 * npar) - (2 * LL)) + (2 * npar * (npar +
   1)/(length(vals) - npar - 1))

 #Check if model has concave curves
 m_betas <- m_aic$betas
 # #Check plot
 # p_aic <- predict.glmnet_mx(m_aic, data, type = "cloglog")
 # plot_aic <- cbind(data, pred = p_aic)
 # plot_aic[,c(5,14)] %>% plot()

 #Select only quadratic betas
 q_betas <- m_betas[grepl("\\^2", names(m_betas))]
 #Check if is concave
 if(length(q_betas) == 0) {
   is_c <- FALSE} else {
     is_c <- any(q_betas > 0)
   }
      }

 #If is concave, write results and check another combination
 if(isTRUE(is_c) | is.na(is_c)){
   ####Save metrics in a dataframe
   all_reg <- unique(formula_grid$regm)
   grid_q <- do.call("rbind",
                     lapply(seq_along(all_reg), function(k){
                       grid_x_i <- grid_x
                       grid_x_i$regm <- all_reg[k]
                       #Check id
                       grid_x_i$ID<- formula_grid[formula_grid$Formulas == grid_x$Formulas & formula_grid$regm == all_reg[k], "ID"]
                       return(grid_x_i)
  }))

   df_eval_q <- empty_replicates(omrat_thr = omrat_thr,
    n_row = nrow(grid_q)*length(k_fold),
    replicates = k_fold, is_c = is_c)

   df_eval_q2 <- cbind(grid_q, df_eval_q)

   #z <- eval_stats(calib_results = df_eval_q2)

   #Summarize results?
   if(only_summary){
     eval_final_q <- empty_summary(omrat_thr = omrat_thr,
      is_c = is_c)
     eval_final_q <- cbind(grid_q, eval_final_q)
   } else {
     eval_final_q <- df_eval_q2 }
 } else {#If is not concave, keep calculating metrics
   mods <- lapply(1:length(k_fold), function(i) {
     data_i <- data[data[, folds] != i, ] #Select i k-fold
     #Run model
     mod_i <- glmnet_mx(p = data_i[,pr_bg], data = data_i,
    f = formula_x, regmult = reg_x,
    calculate_AIC = FALSE)

     #Predict model only to background
     pred_i <- as.numeric(predict(object = mod_i, newdata = data, clamp = FALSE,
     type = "cloglog"))
     pred_i <- cbind(data[,c(pr_bg, folds)], "suit" = pred_i)

     # #Predict model to spatraster, only to check
     # vars_r <- terra::rast("Models/Araucaria_angustifolia/PCA_variables.tiff")
     # pred_r <- terra::predict(vars_r, mod_i, type = "cloglog", na.rm = TRUE)
     # plot(pred_r)

     #Extract suitability in train and test points
     suit_val_cal <- na.omit(pred_i[pred_i[, folds] != i & pred_i[, pr_bg] == 1, "suit"])
     suit_val_eval <- na.omit(pred_i[pred_i[, folds] == i & pred_i[, pr_bg] == 1, "suit"])
     ####Calculate omission rate following kuenm####
     om_rate <- omrat_maxnet(threshold = omrat_thr,
         pred_train = suit_val_cal,
         pred_test = suit_val_eval)
     #Calculate pROC following enmpa
     proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
  prediction = pred_i$suit)


     # ####Calculate pROC following kuenm####
     # proc_i <- proc_maxnet(pred_test = suit_val_eval, model = pred_i$suit, threshold = 5,
     #     rand_percent = 50, iterations = 500, parallel = F)

     ####Save metrics in a dataframe
     df_eval_q <- data.frame(Replicate = i,
         t(data.frame(om_rate)),
         proc_auc_ratio = proc_i$pROC_summary[1],
         proc_pval = proc_i$pROC_summary[2],
         AIC = m_aic$AIC,
         AIC_Warren = AICc,
         npar = npar,
         is_concave = is_c,
         row.names = NULL)
     df_eval_q2 <- cbind(grid_x, df_eval_q)
     return(df_eval_q2)
   })
   #Return evaluaton final
   eval_final_q <- do.call("rbind", mods)

   #Summarize results?
   if(isTRUE(only_summary)){
     eval_final_q <- eval_stats(eval_final_q)
   }
 } #End of else

 #If writeFiles = T...
 if(isTRUE(writeFiles)){
   write.csv(eval_final_q, file.path(out_dir,
        paste0("cand_model_", grid_x$ID, ".csv")),
    row.names = F)
 }

 return(eval_final_q)
        } #End of if( n == 0)
    } #End of first foreach
  }  #End of If test_concave = TRUE

  #Update grid
  if(!test_concave) {n = 0}
  if(test_concave & n > 0){
    #Identify formulas with concave curves
    d_concave <- do.call("rbind", results_concave)
    d_concave <- d_concave[d_concave$is_concave == TRUE, ]
    f_concave <- d_concave$Formulas
    formula_grid2 <- formula_grid[!(formula_grid$ID %in% d_concave$ID), ]
  } else {
    formula_grid2 <- formula_grid
  }

  #Set number of iteration based on new grid
  n <- nrow(formula_grid2)
  #If n == 0, skip non-quadratic models
  if(n == 0) {
    d_res <- d_concave
    warning("All non-quadratic models have been already tested")
  } else {

    #Show progress bar? - Update
    if (isTRUE(progress_bar)) {
      pb <- txtProgressBar(0, n, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n) }

    if (parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    }

    if (parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if (isTRUE(progress_bar))
        opts <- list(progress = progress)
      else opts <- NULL
    }

    if(verbose) {
      if(test_concave) {
        cat("\nTask 2 of 2: calibrating non-quadratic models and quadratic models
        without concave curves\n
        ")
      } else {
        cat("
        Task 1 of 1: calibrating models\n")
      } }
    ####Results non-concave####
    results <- foreach::foreach(x = 1:n, .options.snow = opts,
   .export = to_export) %dopar% {
     #Get grid x
     grid_x <- formula_grid2[x,] #Get i candidate model
     formula_x <- as.formula(grid_x$Formulas) #Get formula from grid x
     reg_x <- grid_x$regm #Get regularization multiplier from grid x

     #Complete model with AIC
     m_aic <- try(glmnet_mx(p = data[,pr_bg], data = data,
          f = formula_x, regmult = reg_x, calculate_AIC = T),
         silent = TRUE)

     if(any(class(m_aic) == "try-error")) {
       npar <- NA
       AICc <- NA
       is_c <- NA
       mods <- NA
       class(mods) <- "try-error"
     } else {

       #Get number of parameters
       npar <- length(m_aic$betas)

       #Calculate AIC from Warren
       vals <- predict.glmnet_mx(m_aic, data[data[,pr_bg] == 1, ],
    type = "exponential")
       LL <- sum(log(vals + .Machine$double.eps))
       #npar <- length(m_aic$betas)
       AICc <- ((2 * npar) - (2 * LL)) + (2 * npar * (npar +
         1)/(length(vals) - npar - 1))

       #Check if model has concave curves
       m_betas <- m_aic$betas
       #Select only quadratic betas
       q_betas <- m_betas[grepl("\\^2", names(m_betas))]
       #Check if is concave
       if(length(q_betas) == 0) {
         is_c <- FALSE} else {
  is_c <- any(q_betas > 0)
         }

       mods <- try(lapply(1:length(k_fold), function(i) {
         data_i <- data[data[, folds] != i, ] #Select i k-fold
         #Run model
         mod_i <- glmnet_mx(p = data_i[,pr_bg], data = data_i,
          f = formula_x, regmult = reg_x,
          calculate_AIC = FALSE)

         #Predict model only to background
         pred_i <- as.numeric(predict(object = mod_i, newdata = data, clamp = FALSE,
         type = "cloglog"))
         pred_i <- cbind(data[,c(pr_bg, folds)], "suit" = pred_i)

         # #Predict model to spatraster, only to check
         # vars_r <- terra::rast("Models/Araucaria_angustifolia/PCA_variables.tiff")
         # pred_r <- terra::predict(vars_r, mod_i, type = "cloglog", na.rm = TRUE)
         # plot(pred_r)

         #Extract suitability in train and test points
         #Extract suitability in train and test points
         suit_val_cal <- na.omit(pred_i[pred_i[, folds] != i &
                                          pred_i[, pr_bg] == 1, "suit"])
         suit_val_eval <- na.omit(pred_i[pred_i[, folds] == i &
                                           pred_i[, pr_bg] == 1, "suit"])
         ####Calculate omission rate following kuenm####
         om_rate <- omrat_maxnet(threshold = omrat_thr,
    pred_train = suit_val_cal,
    pred_test = suit_val_eval)
         #Calculate pROC following enmpa
         proc_i <- enmpa::proc_enm(test_prediction = suit_val_eval,
      prediction = pred_i$suit)


         # ####Calculate pROC following kuenm####
         # proc_i <- proc_maxnet(pred_test = suit_val_eval, model = pred_i$suit, threshold = 5,
         #     rand_percent = 50, iterations = 500, parallel = F)

         ####Save metrics in a dataframe
         df_eval <- data.frame(Replicate = i,
  t(data.frame(om_rate)),
  proc_auc_ratio = proc_i$pROC_summary[1],
  proc_pval = proc_i$pROC_summary[2],
  AIC = m_aic$AIC,
  AIC_Warren = AICc,
  npar = npar,
  is_concave = is_c,
  row.names = NULL)
         df_eval2 <- cbind(grid_x, df_eval)
         return(df_eval2)
       }), silent = TRUE)}

     #####Create empty dataframe if mods is an error####
     if(class(mods) == "try-error") {
       eval_final <- cbind(grid_x, empty_replicates(omrat_thr = omrat_thr,
     n_row = length(k_fold),
     replicates = k_fold, is_c = is_c))
     } else{
       #Return evaluation final
       eval_final <- do.call("rbind", mods) }

     #Summarize results?
     if(isTRUE(only_summary)){
       if(class(mods) == "try-error") {
         eval_final <- cbind(grid_x,
           empty_summary(omrat_thr = omrat_thr, is_c = is_c))
       } else {
         eval_final <- eval_stats(eval_final) }
     }

     #If writeFiles = T...
     if(isTRUE(writeFiles)){
       write.csv(eval_final, file.path(out_dir,
 paste0("cand_model_", grid_x$ID, ".csv")),
        row.names = F)
     }
     return(eval_final) } #End of foreach

    #Stop cluster
    parallel::stopCluster(cl)

    #Convert to dataframe and join results with results concave, if it exists
    d_res <- do.call("rbind", results)
    if(test_concave) {
      d_res <- rbind(d_res, d_concave)
    }

  } #End of if(n == 0)
  return(d_res)
} #End of function
