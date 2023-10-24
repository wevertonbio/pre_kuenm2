source("Functions/Metrics_Functions.R")

models <- readRDS("Models/Myrcia_hatschbachii/Best_models.RDS")
spat_var <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")
consensus_per_model = TRUE
consensus_general = TRUE
consensus = c("median", "range", "mean", "stdev") #mean, weighted mean, variance
write_files = TRUE
write_replicates = FALSE
type = "cloglog"
out_dir = "Models/Myrcia_hatschbachii/Predictions"
overwrite = TRUE

#If predict to other scenarios, keep replicates?

predict_models <- function(models,
                           spat_var,
                           write_files = FALSE,
                           write_replicates = FALSE,
                           out_dir = NULL,
                           consensus_per_model = TRUE,
                           consensus_general = TRUE,
                           consensus = c("median", "range", "mean", "stdev"), #weighted mean
                           type = "cloglog",
                           overwrite = TRUE) {
  #Get models names
  nm <- names(models)

  #Get predictions for each replicate
  p_models <- lapply(models, function(i){
    terra::rast(lapply(i, function(x) {
      terra::predict(spat_var, x, na.rm = TRUE,
                                 type = type)
    }))
  })
  names(p_models) <- nm

  #Start to store results
  res <- list(Replicates = p_models)


  #Get consensus by model
  if(consensus_per_model) {
  if("median" %in% consensus) {
    res$Consensus_per_model$median <- rast(lapply(p_models, median))
  }
  if("mean" %in% consensus) {
    res$Consensus_per_model$mean <- rast(lapply(p_models, mean))
  }
  if("stdev" %in% consensus) {
    res$Consensus_per_model$stdev <- rast(lapply(p_models, stdev))
  }
  if("range" %in% consensus) {
    res$Consensus_per_model$range <- rast(lapply(p_models, function(r) {
      diff(range(r))
    }))
    }
      }

  #Get general consensus
  if(consensus_general & length(p_models) == 1 & consensus_per_model){

    if("median" %in% consensus) {
      res$Consensus_general$median <- res$Consensus_per_model$median
    }
    if("mean" %in% consensus) {
      res$Consensus_general$mean <- res$Consensus_per_model$mean
    }
    if("stdev" %in% consensus) {
      res$Consensus_general$stdev <- res$Consensus_per_model$stdev
    }
    if("range" %in% consensus) {
      res$Consensus_general$range <- res$Consensus_per_model$range
    }
  } else {
    if(consensus_general & length(p_models) == 1 & !consensus_per_model) {
      all_rep <- rast(p_models)

      if("median" %in% consensus) {
        res$Consensus_general$median <- median(all_rep)
      }
      if("mean" %in% consensus) {
        res$Consensus_general$mean <- mean(all_rep)
      }
      if("stdev" %in% consensus) {
        res$Consensus_general$stdev <- stdev(all_rep)
      }
      if("range" %in% consensus) {
        res$Consensus_general$range <- diff(range(all_rep))
      }
    }

    if(consensus_general  & length(p_models) > 1){
      all_rep <- rast(p_models)

      if("median" %in% consensus) {
        res$Consensus_general$median <- median(all_rep)
      }
      if("mean" %in% consensus) {
        res$Consensus_general$mean <- mean(all_rep)
      }
      if("stdev" %in% consensus) {
        res$Consensus_general$stdev <- stdev(all_rep)
      }
      if("range" %in% consensus) {
        res$Consensus_general$range <- diff(range(all_rep))
      }
    }
  }


  # #Get general consensus
  # if(consensus_general){
  #   all_rep <- rast(p_models)
  #
  #   if("median" %in% consensus) {
  #     res$Consensus_general$median <- median(all_rep)
  #   }
  #   if("mean" %in% consensus) {
  #     res$Consensus_general$mean <- mean(all_rep)
  #   }
  #   if("stdev" %in% consensus) {
  #     res$Consensus_general$stdev <- stdev(all_rep)
  #   }
  #   if("range" %in% consensus) {
  #     res$Consensus_general$range <- diff(range(all_rep))
  #     }
  #   }


  #Write files?
  if(write_files) {
    if(!file.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
  #Get names in results
  nmres <- setdiff(names(res), "Replicates")

  #Create folders
  sapply(nmres, function(i){
    dir.create(file.path(out_dir, i))
  })

  #Save files
  sapply(nmres, function(i){
    r_i <- res[[i]]
    sapply(consensus, function(x){
      writeRaster(r_i[[x]],
                  paste0(out_dir, "/", i, "/", x, ".tiff"),
                  overwrite = overwrite) })
  })

  }

  #Write replicates
  if(write_replicates) {
    dir.create(file.path(out_dir, "Replicates"))
    sapply(nm, function(i){
      r <- p_models[[i]]
      writeRaster(r,
                  paste0(out_dir, "/Replicates/", i, ".tiff"),
                  overwrite = overwrite)
    })
    }

  return(res)
  }#End of function

#Predict simple model
source("Functions/Metrics_Functions.R")

models <- readRDS("Models/Piper_fuligineum/Best_models/Best_models.RDS")
spat_var <- rast("Models/Piper_fuligineum/PCA_variables.tiff")

#Predict with consensus per model and consensus general, without write files
pm1 <- predict_models(models = models, spat_var = spat_var,
                     write_files = FALSE,
                     write_replicates = FALSE,
                     out_dir = NULL,
                     consensus_per_model = TRUE,
                     consensus_general = TRUE,
                     consensus = c("median", "range", "mean", "stdev"), #weighted mean
                     type = "cloglog",
                     overwrite = TRUE)

#Predict with only consensus general, without write files
pm2 <- predict_models(models = models, spat_var = spat_var,
                      write_files = FALSE,
                      write_replicates = FALSE,
                      out_dir = NULL,
                      consensus_per_model = FALSE,
                      consensus_general = TRUE,
                      consensus = c("median", "range", "mean", "stdev"), #weighted mean
                      type = "cloglog",
                      overwrite = TRUE)

#Predict with only consensus per model, without write files
pm3 <- predict_models(models = models, spat_var = spat_var,
                      write_files = FALSE,
                      write_replicates = FALSE,
                      out_dir = NULL,
                      consensus_per_model = TRUE,
                      consensus_general = FALSE,
                      consensus = c("median", "range", "mean", "stdev"), #weighted mean
                      type = "cloglog",
                      overwrite = TRUE)
#Predict with consensus per model and consensus general, writing files, but not replicates
pm4 <- predict_models(models = models, spat_var = spat_var,
                      write_files = TRUE,
                      write_replicates = FALSE,
                      out_dir = "Models/Piper_fuligineum/Predictions",
                      consensus_per_model = TRUE,
                      consensus_general = TRUE,
                      consensus = c("median", "range", "mean", "stdev"), #weighted mean
                      type = "cloglog",
                      overwrite = TRUE)
#Predict with consensus per model and consensus general, writing consensus and replicates
pm5 <- predict_models(models = models, spat_var = spat_var,
                      write_files = TRUE,
                      write_replicates = TRUE,
                      out_dir = "Models/Piper_fuligineum/Predictions",
                      consensus_per_model = TRUE,
                      consensus_general = TRUE,
                      consensus = c("median", "range", "mean", "stdev"), #weighted mean
                      type = "cloglog",
                      overwrite = TRUE)
