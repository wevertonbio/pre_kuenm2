#### General workflow to work with (pre)KUENM2 ####
#To start, you need 3 things (4 if you want to project to others scenarios):
 #Dataframe with occurrences of species (with longitude and latitude)
 #Spatraster with the variables to fit the model
 #SpatVector with the M of the specie
 #If deal with projections, tiff files storing the variables to project.
#Each tiff file must correspond to one scenario and has the name of the scenario
#For example: "SSP245_2041-2060_ACCESS-CM2.tiff" if you want to project to the future or "USA.tiff" if you want to project the niche to USA

#Let's start preparing the variables. We will cut the variables using a M, run a PCA to get the x axis that explain 95% of the variance and project this PCAs also to other scenarios (future)

#Before start, load the packages and the necessary functions
#For now, you find the functions in: https://github.com/wevertonbio/pre_kuenm2
library(terra)



####Prepare variables####
source("Functions/prepare_var.R") #Function to prepare variables


#Load the variables
var <- rast(list.files("Current_Neotropic/", full.names = TRUE))
names(var)
#Delete variables that mix temperature and precipitation:
var <- var[[setdiff(names(var), c("Bio08", "Bio09", "Bio18", "Bio18"))]]
names(var)
#Lets run the pca only with the worldclima variables
wc_var <- var[[grepl("Bio", names(var))]]
names(wc_var)
#And let's identify the raw version of the others variables
others_vars <- names(var)[!grepl("Bio", names(var))]

#Get the M of the specie
m <- vect("M_simulations/Myrcia_hatschbachii.gpkg")

#Now, let's run the function prepare_var
prepare_var(spat_var = var, #Variables to fit the model
            mask = m, #M or another mask to crop the variables
            do_PCA = T, #Do PCA?
            exclude_from_PCA = others_vars, #Variables to exclude from PCA
            var_portion = 95, #Get the x axis that explain 95% of the variance
            writeFiles = T, #Write rasters?
            write_PCA_files = T, #Write PCA files? (Model, rotation and importance)
            out_dir = "Models/Myrcia_hatschbachii", #Outuput directory
            overwrite = T, #Overwrite when rasters exists?
            project_PCA = TRUE, #Project PCA?
            proj_folder = "Future_AF/", #Folder with variables to project PCA
            return_proj = FALSE, #Return projections as a object in the environment?
            verbose = TRUE) #Show messages?

#Check the folder set in out_dir :)
#Clean the environment
rm(list = ls())


####Prepare swd####
source("Functions/prepare_data.R") #FUnction to prepare swd

#Now, let's prepare the swd file, which is a file with the values of the variables extracted to each occurrence record of the species, and to n (10.000) background points
#Now, we need the occurrences of the species
occ <- read.csv("Models/Myrcia_hatschbachii/Occurrences.csv")

#Check the column names
colnames(occ)

#Load the PCA variables
pca_var <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")
plot(pca_var) #Check plot

#Prepare SWD
prepare_data(occ = occ,
             species = "species",
             x = "x",
             y = "y",
             spat_variables = pca_var,
             categorical_variables = "SoilType",
             nbg = 10000,
             kfolds = 4,
             include_xy = TRUE,
             write_files = T,
             file_name = "Models/Myrcia_hatschbachii/Data_to_calib",
             seed = 42,
             verbose = TRUE) #Outuput directory
#Check the file occ_bg in the out_dir :)
#Clean the environment
rm(list = ls())

####Fit candidate models####
library(foreach)
source("Functions/calibration_glmnetmx.R") #FUnction to fit candidate models
source("Functions/helpers_calibration_glmnetmx.R") #
source("Functions/calibration_grid_glmnetmx.R") #
source("Functions/glmnet_mx.R")
source("Functions/helpers_glmnetmx.R")
source("Functions/omrat_glmnetmx.R")
source("Functions/select_best_models.R")

#Get swd file
data <- readRDS("Models/Myrcia_hatschbachii/Data.RDS")

#Let's create a grid of formulas, combining differents sets of variables, features and regularization multipliers
f_grid <- calibration_grid_glmnetmx(swd = data$calibration_data, #Swd dataframe
                      x = "x", y = "y", #If swd dataframe cointains x and y, set here; if not, use NULL
                      min.number = 2, #Minimum number of variables in each set of variavles
                      categorical_var = "SoilType", #Identify categorical variables, if exists
                      features = c("l", "lq", "lqp"), #Features
                      regm = c(0.1, 1, 5)) #Regularization multipliers
#See how many candidate models will be testes
nrow(f_grid)

#Fit candidate models
#USe 70% of the available cores
ncores <- round(parallel::detectCores()* 0.7, 0)

m <- calibration_glmnetmx(data = data, #Data in **CLASS??** format
          formula_grid = f_grid[1:50,], #Grid with formulas
          test_convex = TRUE, #Test concave curves in quadratic models?
          parallel = TRUE,
          ncores = 7,
          progress_bar = TRUE, #Show progress bar? Only works if parallelType = "doSNOW"
          write_summary = FALSE, #Write candidate evaluations?
          out_dir = NULL, #Name of the folder to write candidate evaluations
          parallel_type = "doSNOW",
          return_replicate = TRUE,
          omrat_thr = c(5, 10),
          omrat_threshold = 10,
          skip_existing_models = FALSE, #Only works if writeFiles = TRUE
          verbose = TRUE)

#Save candidate models
saveRDS(m, "Models/Myrcia_hatschbachii/candidate_results.RDS")
#Clean the environment
rm(list = ls())


# ####Select the best models#### - Inside Fit candidate models
# #Now, let's select the best models:
# source("Functions/select_best_models.R")
#
# #Import results from candidate models
# m <- readRDS("Models/Myrcia_hatschbachii/candidate_results.RDS")
# bm <- sel_best_models(cand_models = m$Summary, #dataframe with Candidate results
#                       test_convex = T, #Remove models with concave curves?
#                       omrat_threshold = 10, #Omission rate (test points) used to select the best models
#                       allow_tolerance = T, #If omission rate is higher than set, select the model with minimum omission rate
#                       tolerance = 0.01, #If allow tollerance, select the model with minimum omission rate + tolerance
#                       AIC = "nk", #Which AIC? japones (nk) or Warrien (ws?
#                       significance = 0.05, #Significante to select models based on pROC
#                       verbose = TRUE, #Show messages?
#                       delta_aic = 2, #Delta AIC to select best models
#                       save_file = T, #Save file with best models?
#                       file_name = "Models/Myrcia_hatschbachii/best_model") #Outuput directory to save file
# #The result will be a subset of the candidate models with the best models
# #Clean the environment
# rm(list = ls())

#### Fit best models ####
library(terra)
library(foreach)
source("Functions/calibration_glmnetmx.R") #FUnction to fit candidate models
source("Functions/helpers_calibration_glmnetmx.R") #
source("Functions/calibration_grid_glmnetmx.R") #
source("Functions/glmnet_mx.R")
source("Functions/helpers_glmnetmx.R")
source("Functions/omrat_glmnetmx.R")
source("Functions/fit_selected_glmnetmx.R")
source("Functions/part_data.R")

#Import data to fit best model
bm <- readRDS("Models/Myrcia_hatschbachii/candidate_results.RDS") #Best models
#Fit best models: the result will be a RDS file with the replicates of the model
res_kfold <- fit_selected_glmnetmx(calibration_results = bm,
                             #selected_models = bm,
                             # replicates = TRUE,
                             n_replicates = 10,
                             rep_type = "kfold",
                             train_portion = 0.7,
                             write_models = TRUE, #Write files?
                             file_name = "Models/Myrcia_hatschbachii/best_models_kfold", #Name of the folder to write final models
                             parallel = TRUE,
                             ncores = 8,
                             parallelType = "doSNOW",
                             progress_bar = TRUE,
                             verbose = TRUE) #Show messages?
res_subsample <- fit_selected_glmnetmx(calibration_results = bm,
                                   #selected_models = bm,
                                   # replicates = TRUE,
                                   n_replicates = 10,
                                   rep_type = "subsample",
                                   train_portion = 0.7,
                                   write_models = TRUE, #Write files?
                                   file_name = "Models/Myrcia_hatschbachii/best_models_kfold", #Name of the folder to write final models
                                   parallel = TRUE,
                                   ncores = 8,
                                   parallelType = "doSNOW",
                                   progress_bar = TRUE,
                                   verbose = TRUE) #Show messages?
#Clean the environment
rm(list = ls())

####Predict####
library(terra)
source("Functions/predict_selected_glmnetmx.R")
source("Functions/helpers_glmnetmx.R")

res <- readRDS("Models/Myrcia_hatschbachii/best_models_kfold.RDS")
#Test with 2 models
res2 <- c(res, res)
names(res2) <- c("Model_21", "Model_42")

var <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")

p <- predict_selected_glmnetmx(models = res2,
                               spat_var = var,
                    write_files = TRUE,
                    write_replicates = FALSE,
                    out_dir = "Models/Myrcia_hatschbachii/Predictions/",
                    consensus_per_model = TRUE,
                    consensus_general = TRUE,
                    consensus = c("median", "range", "mean", "stdev"), #weighted mean
                    type = "cloglog",
                    overwrite = TRUE)
p
plot(p$Replicates[[1]])
plot(p$Consensus_per_model$median)
plot(p$Consensus_general$median)
