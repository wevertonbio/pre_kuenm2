####Basic workflow to run pre-KUENM2####
#Load packages
library(terra)
library(dplyr)
library(pbapply)
library(pbapply)
library(foreach)


#### 1. Prepare SWD ####
#Load functions
source("Functions/Metrics_Functions.R")

#List species
spp <- list.dirs("Models/", recursive = F, full.names = F)

#Start looping
#At the first time, run the model specie by species, replacing i by 1, 2, 3...
i <- 4

pbapply(seq_along(spp), function(i){
  #Get i specie
  sp <- spp[i]
  #Set directory of the specie
  sp_dir <- file.path("Models", sp)
  #Get occurrence
  occ_sp <- read.csv(file.path(sp_dir, "Occurrences.csv"))
  #Get variables
  var_sp <- rast(file.path(sp_dir, "PCA_variables.tiff"))
  #Prepare swd
  sp_swd <- prepare_swd(occ = occ_sp, species = sp, x = "x", y = "y",
                        spat_var = var_sp, nbg = 10000, kfolds = 4, include_xy = T,
                        writeFiles = T, out_dir = "Models")
})


#### 2. Fit candidate models ####
#Load functions
source("Functions/Metrics_Functions.R")
source("Functions/eval_m.R")

#List species
spp <- list.dirs("Models/", recursive = F, full.names = F)

#Start looping
#At the first time, run the model specie by species, replacing i by 1, 2, 3...
#At the first time, start run in the line 22
i <- 4

pbapply(seq_along(spp), function(i){
  #Get i specie
  sp <- spp[i]
  #Set directory of the specie
  sp_dir <- file.path("Models", sp)
  #Get swd
  occ_bg <- read.csv(file.path(sp_dir, "occ_bg.csv"))
  #Get formulas
  f_grid <- create_grid(#swd = occ_bg,
                        var_names = c("PC1", "PC2", "PC3",
                        "PC4", "PC5"),
                        x = "x", y = "y",
                        min.number = 2,
                        categorical_var = NULL,
                        features = c("l", "lq"),
                        regm = c(0.1, 1, 3, 5))
  #Fit candidate models
  #USe 70% of the available cores
  ncores <- round(parallel::detectCores()* 0.7, 0)
  m <- eval_m(data = occ_bg, pr_bg = "pr_bg", formula_grid = f_grid,
              test_concave = T, folds = "folds",
              parallel = T, ncores = ncores, progress_bar = T, writeFiles = F,
              only_summary = T, omrat_thr = c(5, 10))
  #Save candidate models
  write.csv(m, file.path(sp_dir, "candidate_results.csv"), row.names = F)

})

#### 3. Select best models ####
#Load functions
source("Functions/Metrics_Functions.R")
source("Functions/eval_m.R")
source("Functions/select_best_models.R")

#List species
spp <- list.dirs("Models/", recursive = F, full.names = F)

#Start looping
#At the first time, run the model specie by species, replacing i by 1, 2, 3...
#At the first time, start run in the line 22
i <- 4

pbapply(seq_along(spp), function(i){
  #Get i specie
  sp <- spp[i]
  #Set directory of the specie
  sp_dir <- file.path("Models", sp)
  #Get candidate results
  cand_res <- read.csv(file.path(sp_dir, "candidate_results.csv"))
  #Select best model
  bm <- sel_best_models(cand_models = cand_res,
                        test_concave = T, omrat = 5, omrat_threshold = 5,
                        allow_tolerance = T, tolerance = 0.01, AIC = "AIC_Warren",
                        significance = 0.05,
                        verbose = TRUE,
                        save_file = F,
                        output_dir = NULL)

  #Save candidate models
  write.csv(bm, file.path(sp_dir, "selected_models.csv"), row.names = F)
  })
