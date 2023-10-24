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
source("Functions/prepare_swd.R") #FUnction to prepare swd
source("Functions/Metrics_Functions.R") #Helper functions

#Now, let's prepare the swd file, which is a file with the values of the variables extracted to each occurrence record of the species, and to n (10.000) background points
#Now, we need the occurrences of the species
occ <- read.csv("Models/Myrcia_hatschbachii/Occurrences.csv")

#Check the column names
colnames(occ)

#Load the PCA variables
pca_var <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")
plot(pca_var) #Check plot

#Prepare SWD
prepare_swd(occ = occ, #Dataframe with the occurrences
            x = "x", y = "y", #Column names with longitude (x) and latitude (y)
            species = "Myrcia_hatschbachii", #Name of the specie (optional)
            spat_var = pca_var, #Variables to fit the model (here, PCA variable)
            nbg = 10000, #Number of background points
            kfolds = 4, #Number of k-folds to validate the candidate models
            include_xy = TRUE, #Include coordinates in the swd file?
            writeFiles = TRUE, #Write swd file?
            out_dir = "Models/Myrcia_hatschbachii") #Outuput directory
#Check the file occ_bg in the out_dir :)
#Clean the environment
rm(list = ls())

####Fit candidate models####
library(kuenm)
library(foreach)
source("Functions/eval_m.R") #FUnction to fit candidate models
source("Functions/Metrics_Functions.R") #Helper functions

#Get swd file
occ_bg <- read.csv("Models/Myrcia_hatschbachii/occ_bg.csv")
#Check column names
colnames(occ_bg)
#Let's create a grid of formulas, combining differents sets of variables, features and regularization multipliers
f_grid <- create_grid(swd = occ_bg, #Swd dataframe
                      x = "x", y = "y", #If swd dataframe cointains x and y, set here; if not, use NULL
                      min.number = 2, #Minimum number of variables in each set of variavles
                      categorical_var = "SoilType", #Identify categorical variables, if exists
                      features = c("l", "lq", "lqp"), #Features
                      regm = c(0.1, 1, 3, 5)) #Regularization multipliers
#See how many candidate models will be testes
nrow(f_grid)

#Fit candidate models
#USe 70% of the available cores
ncores <- round(parallel::detectCores()* 0.7, 0)

m <- eval_m(data = occ_bg, #SWD dataframe
            pr_bg = "pr_bg", #Column name in data indicating presence/background
            formula_grid = f_grid, #Formula grid
            var_categorical = "SoilType", #Identify categorical variables, if exists
            test_concave = T, #Test concave curves?
            folds = "folds", #Column name in data indicating folds to evaluate the model
            parallel = T, #Run in parallel?
            ncores = ncores, #Number of cores
            progress_bar = T, #Show progress bar?
            writeFiles = F, #Write files?
            only_summary = T, #Return results for all folds or only summary?
            omrat_thr = c(5, 10)) #Omission rate to evaluate the models

#Save candidate models
write.csv(m, "Models/Myrcia_hatschbachii/candidate_results.csv", row.names = F)
#Clean the environment
rm(list = ls())


####Select the best models####
#Now, let's select the best models:
source("Functions/select_best_models.R")

#Import results from candidate models
m <- read.csv("Models/Myrcia_hatschbachii/candidate_results.csv")
bm <- sel_best_models(cand_models = m, #dataframe with Candidate results
                      test_concave = T, #Remove models with concave curves?
                      omrat = 10, #Omission rate (train points) used to select the best models
                      omrat_threshold = 10, #Omission rate (test points) used to select the best models
                      allow_tolerance = T, #If omission rate is higher than set, select the model with minimum omission rate
                      tolerance = 0.01, #If allow tollerance, select the model with minimum omission rate + tolerance
                      AIC = "japones", #Which AIC? japones or Warrien?
                      significance = 0.05, #Significante to select models based on pROC
                      verbose = TRUE, #Show messages?
                      delta_aic = 2, #Delta AIC to select best models
                      save_file = T, #Save file with best models?
                      output_dir = "Models/Myrcia_hatschbachii") #Outuput directory to save file
#The result will be a subset of the candidate models with the best models
#Clean the environment
rm(list = ls())

#### Fit best models ####
library(terra)
library(foreach)
source("Functions/fit_best.R")
source("Functions/Metrics_Functions.R")
source("Functions/part_data.R")

#Import data to fit best model
occ_bg <- read.csv("Models/Myrcia_hatschbachii/occ_bg.csv") #SWD file
bm <- read.csv("Models/Myrcia_hatschbachii/selected_models.csv") #Best models
#Fit best models: the result will be a RDS file with the replicates of the model
res <- fit_best(data = occ_bg, #SWD dataframe
                       pr_bg = "pr_bg", #Column name in data indicating presence/background
                       selected_models = bm, #Best models (output of sel_best_models)
                       var_categorical = "SoilType", #Identify categorical variables, if exists
                       replicates = TRUE, #Run replicates?
                       n_replicates = 10, #Number of replicates
                       rep_type = "subsample", #Type of partition: subsample, bootstrap or kfold
                       train_portion = 0.7, #If subsample or bootstrap, portion of train points
                       parallel = TRUE, #Run in parallel?
                       ncores = 8, #Number of cores
                       progress_bar = TRUE, #Show progress bar?
                       write_models = TRUE, #Write files?
                       out_dir = "Models/Myrcia_hatschbachii", #Name of the folder to write final models
                       parallelType = "doSNOW", #Paralell type (doSnow or doParallel)
                       verbose = TRUE) #Show messages?
#Clean the environment
rm(list = ls())

####Predict####

var <- rast("Models/Myrcia_hatschbachii/PCA_var.tiff")
p <- predict_models(models = res, spat_var = var,
                    write_files = FALSE,
                    write_replicates = FALSE,
                    out_dir = NULL,
                    consensus_per_model = FALSE,
                    consensus_general = TRUE,
                    consensus = c("median", "range", "mean", "stdev"), #weighted mean
                    type = "cloglog",
                    overwrite = TRUE)
p
plot(p$Replicates$Model_771)
plot(p$Consensus_per_model$median)
plot(p$Consensus_general$median)
