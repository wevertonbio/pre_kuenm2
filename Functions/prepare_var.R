#Load packages
#library(terra)
#library(dplyr)

#Example data
# var <- list.files("C:/Users/wever/Desktop/KU_Models/Current_Neotropic/",
#                   full.names = T) %>% rast()
# names(var)
# var <- var[[setdiff(names(var), c("Altitude", "Nitrogen",
#                            "OrganicCArbon", "ph", "Ruggedness",
#                            "TopoIndex"))]]
# m <- vect("C:/Users/wever/Desktop/KU_Models/M_simulations/New_M/Araucaria_angustifolia.gpkg")
# spat_var = var
# mask = m
# do_PCA = TRUE
# exclude_from_PCA = names(var)[!grepl("Bio", names(var))]
# writeFiles = F
# out_dir = "Models/Araucaria_angustifolia/"
# var_portion = 95
# write_PCA_files = T

prepare_var <- function(spat_var,
                        mask = NULL,
                        do_PCA = F,
                        exclude_from_PCA = NULL,
                        var_portion = 95,
                        writeFiles = F,
                        write_PCA_files = F,
                        out_dir = NULL,
                        overwrite = T){
  if(!is.null(mask)){
    spat_var <- terra::crop(spat_var, mask, mask = TRUE)
  }

  if(do_PCA){
    if(!is.null(exclude_from_PCA)) {
      spat_var_pca <- spat_var[[setdiff(names(spat_var),exclude_from_PCA)]]
    } else {
      spat_var_pca <- spat_var}
    #Do PCA
    #PCA
    p_env <- stats::prcomp(as.data.frame(spat_var_pca),
                           retx = TRUE,
                           scale. = TRUE,
                           center = TRUE)

    #Predict to current
    vars_now <- terra::predict(spat_var_pca, p_env)
    #Select n-variables based on portion of explained variance
    var_exp <- t(summary(p_env)$importance)
    var_sel <- which(var_exp[,3] > var_portion/100)[1]
    vars_now <- vars_now[[1:var_sel]]

    if(!is.null(exclude_from_PCA)) {
      spat_final <- c(vars_now, spat_var[[exclude_from_PCA]])
    } else {spat_final <- vars_now }

    #If writeFiles
    if(writeFiles){
      writeRaster(spat_final, file.path(out_dir, "PCA_var.tiff"),
                  overwrite = overwrite)
    }

    if(write_PCA_files){
      dir.create(file.path(out_dir, "PCA_results"))
      #Get coefficients
      cof_env <- p_env$rotation %>% as.data.frame() %>%
        mutate(Variable = row.names(.), .before = 1)
      write.csv(cof_env, file.path(out_dir, "PCA_results/PCA_rotation.csv"))
      saveRDS(p_env, file.path(out_dir, "PCA_results/PCA_model.RDS"))
      write.csv(var_exp, file.path(out_dir, "PCA_results/PCA_importance.csv"))

    }
  } else {
    spat_final <- spat_var
  }

  return(spat_final)
}

# #Test function
# #Example data
# var <- list.files("C:/Users/wever/Desktop/KU_Models/Current_Neotropic/",
#                   full.names = T) %>% rast()
# names(var)
# var <- var[[setdiff(names(var), c("Altitude", "Nitrogen",
#                                   "OrganicCArbon", "ph", "Ruggedness",
#                                   "TopoIndex"))]]
# m <- vect("C:/Users/wever/Desktop/KU_Models/M_simulations/New_M/Araucaria_angustifolia.gpkg")
#
# pca_var <- prepare_var(spat_var = var, #SpatRaster of variables
#                        mask = m, #Spatvector of the mask/M
#                        do_PCA = T, #Do PCA?
#                        #Var to exclude from PCA:
#                        exclude_from_PCA = names(var)[!grepl("Bio",
#                                                             names(var))],
#                        var_portion = 95, #Select axis that explained x% of the variance
#                        writeFiles = T, #Write files?
#                        write_PCA_files = F, #Write PCA files (model, importance and rotation?)
#                        out_dir = NULL, #Output directory
#                        overwrite = T) #Overwrite files
