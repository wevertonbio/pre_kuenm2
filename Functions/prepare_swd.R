#### Prepare SWD ####
#library(terra)
prepare_swd <- function(occ,
                        species = NULL,
                        x,
                        y,
                        spat_var,
                        nbg = 10000,
                        kfolds = 4,
                        include_xy = TRUE,
                        writeFiles = F,
                        out_dir = NULL){
  #Extract variables to occurrences
  xy <- as.matrix(occ[,c(x,y)])
  colnames(xy) <- c("x", "y")
  occ_var <- cbind(xy, terra::extract(x = spat_var, y = xy))
  occ_var$pr_bg <- 1
  #Get background points
  cell_samp <- terra::as.data.frame(spat_var[[1]], na.rm = TRUE,
                                    cells = TRUE)[, "cell"]
  cell_samp <- sample(cell_samp, size = nbg, replace = FALSE,
                      prob = NULL)
  bg_var <- terra::extract(x = spat_var, y = cell_samp, xy = TRUE)
  bg_var$pr_bg <- 0
  #Join occ and background
  occ_bg <- rbind(occ_var, bg_var)
  #Reorder columns
  occ_bg <- occ_bg[,c(1, 2, which(names(occ_bg) == "pr_bg"),
                      (3:(ncol(occ_bg)-1))
                      [-which(names(occ_bg) == "pr_bg")])]
  #Include xy?
  if(!include_xy){
    occ_bg <- subset(occ_bg, select = c(-x, -y))
  } else {
    colnames(occ_bg)[colnames(occ_bg) %in% c("x", "y")] <- c(x, y)
  }

  #Split in kfolds?
  if(is.numeric(kfolds)) {
    occ_bg <- kfold_part(data = occ_bg, n_folds = kfolds, pr_bg = "pr_bg")
    #Reorder columns
    pos_pr_bg <- which(names(occ_bg) == "pr_bg")
    occ_bg <- occ_bg[, c(1:(pos_pr_bg), ncol(occ_bg), (pos_pr_bg + 1):(ncol(occ_bg)-1))]
  }

  #Append species name?
  if(!is.null(species)){
    occ_bg$species <- species
    #Reorder columns
    occ_bg <- occ_bg[, c("species",
                         colnames(occ_bg)[colnames(occ_bg) != "species"])]
  }
  #Save results?
  if(writeFiles) {
    if(is.null(out_dir)){
      out_dir <- getwd()
    } else{
      if(!file.exists(out_dir))
        dir.create(out_dir, recursive = TRUE)
    }
    write.csv(occ_bg, file.path(out_dir, "occ_bg.csv"),
              row.names = FALSE)
  }
  return(occ_bg)
}
