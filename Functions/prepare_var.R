#' Prepare environmental variables
#'
#' @description
#' Prepare variables to run models, masking and/or generating PCA variables
#'
#' @param variables
#' @param mask
#' @param do_PCA
#' @param exclude_from_PCA
#' @param var_portion
#' @param writeFiles
#' @param write_PCA_files
#' @param out_dir
#' @param overwrite
#' @param project_PCA
#' @param proj_spat
#' @param proj_folder
#' @param return_proj
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
prepare_var <- function(variables,
                        calib_area = NULL,
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
  if(!is.null(calib_area)){
    variables <- terra::crop(variables, calib_area, mask = TRUE)
  }

  if(do_PCA){
    if(!is.null(exclude_from_PCA)) {
      variables_pca <- variables[[setdiff(names(variables),exclude_from_PCA)]]
    } else {
      variables_pca <- variables}
    #Do PCA
    #PCA
    p_env <- stats::prcomp(as.data.frame(variables_pca),
                           retx = TRUE,
                           scale. = TRUE,
                           center = TRUE)

    #Predict to current
    vars_now <- terra::predict(variables_pca, p_env)
    #Select n-variables based on portion of explained variance
    var_exp <- t(summary(p_env)$importance)
    var_sel <- which(var_exp[,3] > var_portion/100)[1]
    vars_now <- vars_now[[1:var_sel]]

    if(!is.null(exclude_from_PCA)) {
      spat_final <- c(vars_now, variables[[exclude_from_PCA]])
    } else {spat_final <- vars_now }

    #If writeFiles
    if(writeFiles){
      writeRaster(spat_final, file.path(out_dir, "PCA_var.tiff"),
                  overwrite = overwrite)
    }

    if(write_PCA_files){
      dir.create(file.path(out_dir, "PCA_results"), recursive = TRUE)
      #Get coefficients
      rotation_matrix <- as.data.frame(p_env$rotation)

      # Add the "Variable" column at the beginning of the rotation matrix
      cof_env <- cbind(Variable = row.names(rotation_matrix), rotation_matrix)

      write.csv(cof_env, file.path(out_dir, "PCA_results/PCA_rotation.csv"))
      saveRDS(p_env, file.path(out_dir, "PCA_results/PCA_model.RDS"))
      write.csv(var_exp, file.path(out_dir, "PCA_results/PCA_importance.csv"))}
  } else {
    spat_final <- variables
  }

  #If project_PCA = TRUE...
  if(project_PCA) {
    if(writeFiles)
      dir.create(file.path(out_dir, "Variable_projections"), recursive = TRUE)
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
        psx <- proj_spat[[name_x]][[names(variables_pca)]]
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
