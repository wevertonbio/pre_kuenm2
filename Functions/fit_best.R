#Load packages
library(terra)
library(dplyr)
library(pbapply)
library(pbapply)
library(foreach)
library(parallel)

#Load functions
source("Functions/Metrics_Functions.R")
source("Functions/eval_m.R")
source("Functions/part_data.R")

#Prediction to dataframe or raster
#In df, each replicate in a column

#Select 3 best models from candidate results
cand_res <- read.csv("Models/Piper_fuligineum/candidate_results.csv")
selected_models <- cand_res %>%
  filter(is_concave == F, proc_pval.mean < 0.05,
         Omission_rate_at_5.mean <= 0.05) %>% slice_min(AIC, n =3)

#Import data to create the functions
data <- read.csv("Models/Piper_fuligineum/occ_bg.csv")
pr_bg <- "pr_bg"
var_categorical = NULL
replicates = TRUE
n_replicates <- 5
rep_type = "subsample"
train_portion = 0.7
write_models = TRUE
out_dir = "Models/Piper_fuligineum/Best_models/" #Name of the folder to write final models
parallel = TRUE
ncores = 1
progress_bar = TRUE
parallelType = "doSNOW"
verbose = TRUE
to_export = c("aic_glmnetmx", "aic_maxnet", "eval_stats",
              "get_formulas_maxnet",
              "glmnet_mx", "kfold_part",
              "maxnet.default.regularization",
              "omrat_maxnet",
              "predict.glmnet_mx", "empty_replicates",
              "empty_summary",
              "hinge", "hingeval", "thresholds",
              "thresholdval", "categorical",
              "categoricalval")




####Function to fit best models####
fit_best <- function(data, #Data in SWD format
                     pr_bg,
                     var_categorical = NULL,
                     selected_models,
                     replicates = TRUE,
                     n_replicates = 10,
                     rep_type = "kfold",
                     train_portion = 0.7,
                     write_models = FALSE, #Write files?
                     out_dir = NULL, #Name of the folder to write final models
                     parallel = TRUE,
                     ncores = 1,
                     parallelType = "doSNOW",
                     progress_bar = TRUE,
                     verbose = TRUE,
                     to_export = c("aic_glmnetmx", "aic_maxnet", "eval_stats",
                                   "get_formulas_maxnet",
                                   "glmnet_mx", "kfold_part",
                                   "maxnet.default.regularization",
                                   "omrat_maxnet",
                                   "predict.glmnet_mx", "empty_replicates",
                                   "empty_summary",
                                   "hinge", "hingeval", "thresholds",
                                   "thresholdval", "categorical",
                                   "categoricalval")) {

  #If writefiles = TRUE, create directory
  if(write_models){
    if(!file.exists(out_dir))
      dir.create(out_dir, recursive = T)
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

  #Extracts IDs from models
  m_ids <- selected_models$ID

  #If replicates = FALSE, set n_replicates = 1
  if(!replicates) {
    n_replicates <- 1
  }

  #Create grid of fitted models
  dfgrid <- expand.grid(models = m_ids, replicates = 1:n_replicates)
  n <-  nrow(dfgrid)

  #Prepare data (index) to replicates
  if(replicates) {
    #Partitioning data
    rep_data <- part_data(data = data, pr_bg = pr_bg,
                          train_portion = train_portion,
                          n_replicates = n_replicates,
                          method = rep_type) }


  # #Fit models with no replicates
  # if(!replicates) {
  #   best_models <- lapply(1:nrow(selected_models), function(i){
  #     best_models_i <- selected_models[i,]
  #     best_formula <- best_models_i$Formulas
  #     best_regm <- best_models_i$regm
  #     m_i <- glmnet_mx(p = data[,pr_bg], data = data,
  #                      f = as.formula(best_formula),
  #                      regmult = best_regm, calculate_AIC = F)
  #     return(m_i)
  #   })
  #   names(best_models) <- paste0("Model_", m_ids)
  # } #End of !replicates



  #Set parallelization
  if(parallel) {
    #Make cluster
    cl <- parallel::makeCluster(ncores)
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
  } else {opts <- NULL}

  #Fit models with replicates
  best_models <- foreach::foreach(x = 1:n, .options.snow = opts,
                                  .export = to_export) %dopar% {
    #Get grid
    grid_x <- dfgrid[x,]
    m_id <- grid_x$models
    rep_x <- grid_x$replicates

    #Get best model
    best_models_i <- selected_models[which(selected_models$ID == m_id),]
    #best_models_i <- selected_models[i,]
    best_formula <- best_models_i$Formulas
    best_regm <- best_models_i$regm

    #Get replicate, if necessary
    if(replicates){
    rep_i <- rep_data[[rep_x]]
    data_x <- data[rep_i, ] } else { #Select i k-fold
      data_x <- data }
    #Run model
    mod_x <- glmnet_mx(p = data_x[,pr_bg], data = data_x,
                       f = as.formula(best_formula),
                       regmult = best_regm,
                       calculate_AIC = FALSE)
    #Only to test

    mod_x$checkModel <- m_id
    mod_x$checkreplicate <- rep_x
    # #Predict
    # pred_x <- terra::predict(spat_var, mod_x, type = "cloglog",
    #                          na.rm = TRUE)
    return(mod_x) } #End of foreach

  #End cluster
  parallel::stopCluster(cl)

  #Split list
  # Crie um vetor com o número de réplicas para cada modelo
  num_repl <- tapply(dfgrid$replicates, dfgrid$models, FUN = length)
  num_repl <- num_repl[match(dfgrid$models, names(num_repl))]

  # Split list
  best_models <- split(best_models, dfgrid$models)
  best_models <- lapply(best_models, function(sublist) {
    names(sublist) <- paste0("Rep_", seq_along(sublist))
    return(sublist)
  })

  #Rename models
  names(best_models) <- paste0("Model_", names(best_models))

  #Write models?
  if(write_models){
    saveRDS(best_models, file.path(out_dir, "Best_models.RDS"))
    if(verbose){
      message("Best_models.RDS saved in ", out_dir)
    }
  }

  return(best_models)
} #End of function

#Test function
#Load functions
source("Functions/Metrics_Functions.R")
source("Functions/eval_m.R")
source("Functions/part_data.R")

#Prediction to dataframe or raster
#In df, each replicate in a column

#Select 3 best models from candidate results
cand_res <- read.csv("Models/Piper_fuligineum/candidate_results.csv")
selected_models <- cand_res %>%
  filter(is_concave == F, proc_pval.mean < 0.05,
         Omission_rate_at_5.mean <= 0.05) %>% slice_min(AIC, n =3)
occ_bg <- read.csv("Models/Piper_fuligineum/occ_bg.csv")

#With kfold
best_model <- fit_best(data = occ_bg,
                       pr_bg = "pr_bg",
                       selected_models = selected_models,
                       var_categorical = NULL,
                       replicates = TRUE,
                       n_replicates = 4,
                       train_portion = 0.7,
                       rep_type = "kfold",
                       parallel = TRUE,
                       ncores = 5,
                       progress_bar = TRUE,
                       write_models = FALSE, #Write files?
                       out_dir = NULL, #Name of the folder to write final models
                       parallelType = "doSNOW",
                       verbose = TRUE)

#With subsample
best_model2 <- fit_best(data = occ_bg,
                       pr_bg = "pr_bg",
                       selected_models = selected_models,
                       var_categorical = NULL,
                       replicates = TRUE,
                       n_replicates = 5,
                       train_portion = 0.7,
                       rep_type = "subsample",
                       parallel = TRUE,
                       ncores = 5,
                       progress_bar = TRUE,
                       write_models = FALSE, #Write files?
                       out_dir = NULL, #Name of the folder to write final models
                       parallelType = "doSNOW",
                       verbose = TRUE)
#With bootstrap
best_model3 <- fit_best(data = occ_bg,
                        pr_bg = "pr_bg",
                        selected_models = selected_models,
                        var_categorical = NULL,
                        replicates = TRUE,
                        n_replicates = 10,
                        train_portion = 0.7,
                        rep_type = "bootstrap",
                        parallel = TRUE,
                        ncores = 5,
                        progress_bar = TRUE,
                        write_models = FALSE, #Write files?
                        out_dir = NULL, #Name of the folder to write final models
                        parallelType = "doSNOW",
                        verbose = TRUE)
#With kfold, no parallel and writing results
best_model4 <- fit_best(data = occ_bg,
                        pr_bg = "pr_bg",
                        selected_models = selected_models,
                        var_categorical = NULL,
                        replicates = TRUE,
                        n_replicates = 4,
                        train_portion = 0.7,
                        rep_type = "kfold",
                        parallel = T,
                        ncores = 5,
                        progress_bar = TRUE,
                        write_models = TRUE, #Write files?
                        out_dir = "Models/Piper_fuligineum/Best_models/", #Name of the folder to write final models
                        parallelType = "doSNOW",
                        verbose = TRUE)
