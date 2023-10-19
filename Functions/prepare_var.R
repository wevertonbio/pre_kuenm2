# #Load packages
# library(terra)
# library(dplyr)
#
# #Example data
# var <- list.files("C:/Users/wever/Desktop/KU_Models/Current_Neotropic/",
#                   full.names = T) %>% rast()
# names(var)
# var <- var[[setdiff(names(var), c("Altitude", "Nitrogen",
#                                   "OrganicCArbon", "ph", "Ruggedness",
#                                   "TopoIndex"))]]
# m <- vect("C:/Users/wever/Desktop/KU_Models/M_simulations/New_M/Araucaria_angustifolia.gpkg")
#
# #Future layers
# fvar <- list(rast("C:/Users/wever/Desktop/KU_Models/Future_AF/SSP245_2041-2060_ACCESS-CM2.tiff"),
#              rast("C:/Users/wever/Desktop/KU_Models/Future_AF/SSP245_2041-2060_IPSL-CM6A-LR.tiff"))
# #Rename
# sources_fvar <- sapply(fvar, terra::sources)
# names(fvar) <- sub(".*/(.*?)\\.tif*", "\\1", sources_fvar)
#
#
#
# spat_var = var
# mask = m
# do_PCA = TRUE
# exclude_from_PCA = names(var)[!grepl("Bio", names(var))]
# writeFiles = TRUE
# out_dir = "Models/Araucaria_angustifolia/"
# var_portion = 95
# write_PCA_files = T
# overwrite = T
# project_PCA = TRUE
# proj_spat = fvar
# proj_folder = "C:/Users/wever/Desktop/KU_Models/Future_AF/"
# verbose = TRUE

prepare_var <- function(spat_var,
                        mask = NULL,
                        do_PCA = F,
                        exclude_from_PCA = NULL,
                        var_portion = 95,
                        writeFiles = F,
                        write_PCA_files = F,
                        out_dir = NULL,
                        overwrite = TRUE,
                        project_PCA = FALSE,
                        proj_spat = NULL, #Need to be a list with SpatRast. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
                        proj_folder = NULL, #Variables need to be stacked in a single tiff file by scenario, and file need to be the name of the scenario. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
                        return_proj = TRUE,
                        verbose = TRUE) {
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
      write.csv(var_exp, file.path(out_dir, "PCA_results/PCA_importance.csv"))}
  } else {
    spat_final <- spat_var
  }

  #If project_PCA = TRUE...
  if(project_PCA) {
    if(writeFiles)
      dir.create(file.path(out_dir, "Variable_projections"))
    if(verbose & writeFiles)
      message("PCA projections will be saved in: ",
              file.path(out_dir, "Variable_projections"))

    #If spatrasters are in folders
    if(!is.null(proj_folder)) {
      #List tiff files in directory
      list_proj <- list.files(path = proj_folder, pattern = "\\.tiff$|\\.tif$",
                              full.names = TRUE, recursive = TRUE)
      #Rasterize files and save in a list
      proj_spat <- lapply(list_proj, terra::rast)
      #Rename scenarios
      names(proj_spat) <- sub(".*/(.*?)\\.tif*", "\\1",
                              sapply(proj_spat, sources))
    }

    #If spatrasters are in a list
    if(!is.null(proj_spat)) {
      r_proj <- lapply(names(proj_spat), function(x){
        #Subset variables
        name_x <- x
        psx <- proj_spat[[name_x]][[names(spat_var_pca)]]
        #Project PCA
        psr <- terra::predict(psx, p_env, na.rm = TRUE)
        #Select axix
        psr <- psr[[1:var_sel]]

        #If there is variable to exclude from PCA
        if(!is.null(exclude_from_PCA)) {
          var_proj_exclude <- proj_spat[[name_x]][[exclude_from_PCA]]
          psr <- c(psr, var_proj_exclude)
        }
        #Save
        if(writeFiles)
          writeRaster(psr,
                      paste0(file.path(out_dir, "Variable_projections"),
                             "/", name_x, ".tiff"), overwrite = overwrite)
        if(return_proj) {
          return(psr) } else {
            return(NULL)}
      })
      if(return_proj)
        names(r_proj) <- names(proj_spat)

    } #End of if(!is.null(proj_spat))

  } #End of if(project_PCA)

  #See what return
  if(project_PCA & return_proj) {
    result <- list(Current = spat_final,
                   Projections = r_proj)
  } else {return(spat_final)}
} #End of function

# #Test function
# #Example data
# var <- list.files("C:/Users/wever/Desktop/KU_Models/Current_Neotropic/",
#                   full.names = T) %>% rast()
# names(var)
# var <- var[[setdiff(names(var), c("Altitude", "Nitrogen",
#                                   "OrganicCArbon", "ph", "Ruggedness",
#                                   "TopoIndex"))]]
# m <- vect("C:/Users/wever/Desktop/KU_Models/M_simulations/New_M/Araucaria_angustifolia.gpkg")
# exclude_from_PCA = names(var)[!grepl("Bio", names(var))]
#
#
# ####Test without prediction and not write
# pca_var <- prepare_var(spat_var = var,
#                        mask = m,
#                        do_PCA = TRUE,
#                        exclude_from_PCA = exclude_from_PCA,
#                        var_portion = 95,
#                        writeFiles = FALSE,
#                        write_PCA_files = T,
#                        out_dir = "Models/Araucaria_angustifolia/",
#                        overwrite = TRUE,
#                        project_PCA = FALSE,
#                        proj_spat = NULL, #Need to be a list with SpatRast. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                        proj_folder = NULL, #Variables need to be stacked in a single tiff file by scenario, and file need to be the name of the scenario. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                        return_proj = TRUE,
#                        verbose = TRUE) #Overwrite files
#
# ####Test with predictions on list, but not write data
# #Future layers
# fvar <- list(rast("C:/Users/wever/Desktop/KU_Models/Future_AF/SSP245_2041-2060_ACCESS-CM2.tiff"),
#              rast("C:/Users/wever/Desktop/KU_Models/Future_AF/SSP245_2041-2060_IPSL-CM6A-LR.tiff"))
# #Rename
# sources_fvar <- sapply(fvar, terra::sources)
# names(fvar) <- sub(".*/(.*?)\\.tif*", "\\1", sources_fvar)
#
# pca_var2 <- prepare_var(spat_var = var,
#                         mask = m,
#                         do_PCA = TRUE,
#                         exclude_from_PCA = exclude_from_PCA,
#                         var_portion = 95,
#                         writeFiles = FALSE,
#                         write_PCA_files = T,
#                         out_dir = "Models/Araucaria_angustifolia/",
#                         overwrite = TRUE,
#                         project_PCA = TRUE, #If project_PCA is TRUE, need to define out_dir
#                         proj_spat = fvar, #Need to be a list with SpatRast. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                         proj_folder = NULL, #Variables need to be stacked in a single tiff file by scenario, and file need to be the name of the scenario. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                         return_proj = TRUE,
#                         verbose = TRUE) #Overwrite files
#
# ####Test with predictions on folders, and write data
# #Future layers
# fvar <- "C:/Users/wever/Desktop/KU_Models/Future_AF/"
#
# pca_var3 <- prepare_var(spat_var = var,
#                         mask = m,
#                         do_PCA = TRUE,
#                         exclude_from_PCA = exclude_from_PCA,
#                         var_portion = 95,
#                         writeFiles = TRUE,
#                         write_PCA_files = T,
#                         out_dir = "Models/Araucaria_angustifolia/",
#                         overwrite = TRUE,
#                         project_PCA = TRUE, #If project_PCA is TRUE, need to define out_dir
#                         proj_spat = NULL, #Need to be a list with SpatRast. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                         proj_folder = fvar, #Variables need to be stacked in a single tiff file by scenario, and file need to be the name of the scenario. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                         return_proj = TRUE,
#                         verbose = TRUE) #Overwrite files
#
# ####Test with predictions on folders, write data and do not return projections
# #Future layers
# fvar <- "C:/Users/wever/Desktop/KU_Models/Future_AF/"
#
# pca_var4 <- prepare_var(spat_var = var,
#                         mask = m,
#                         do_PCA = TRUE,
#                         exclude_from_PCA = exclude_from_PCA,
#                         var_portion = 95,
#                         writeFiles = TRUE,
#                         write_PCA_files = T,
#                         out_dir = "Models/Araucaria_angustifolia/",
#                         overwrite = TRUE,
#                         project_PCA = TRUE, #If project_PCA is TRUE, need to define out_dir
#                         proj_spat = NULL, #Need to be a list with SpatRast. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                         proj_folder = fvar, #Variables need to be stacked in a single tiff file by scenario, and file need to be the name of the scenario. If there are variables to exclude from PCA (for example, soil type), this variable must be included in the tiff file
#                         return_proj = FALSE,
#                         verbose = TRUE) #Overwrite files
