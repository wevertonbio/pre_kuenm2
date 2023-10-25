####Ver como funciona categorical####





####Metrics functions####









# #Test function
# #Example
# #Get species list
# sp <- list.dirs(path = "Models", recursive = F, full.names = F)
# sp <- sp[1]
# #Get directory of specie
# dir.sp <- file.path("Models/", sp)
# ##Get occurrences
# occ_sp <- read.csv(file.path(dir.sp, "Occurrences.csv"))
# ##Get spatvar
# #Get variables
# vars_now <- rast(file.path(dir.sp, "PCA_variables.tiff"))
# swd <- prepare_swd(occ = occ_sp, x = "x", y = "y",
#                    spat_var = vars_now, nbg = 10000, kfolds = 4, include_xy = F)






# ####PROC####
# #Not necessary, use enmpa::proc_enm
# proc_maxnet <- function(pred_test, model, threshold = 5,
#                         rand_percent = 50,
#                         iterations = 500, parallel = FALSE) {
#   min_pred <- min(model, na.rm = T)
#   max_pred <- max(model, na.rm = T)
#   vals <- na.omit(model)
#   nvals <- length(vals)
#   vals <- c(vals, pred_test)
#   vals <- as.numeric(cut(vals, 500))
#   pred_test <- vals[(nvals + 1):length(vals)]
#   vals <- vals[1:nvals]
#   classpixels <- as.data.frame(table(vals), stringsAsFactors = FALSE)
#   colnames(classpixels) <- c("value", "count")
#   if (min_pred == max_pred) {
#     warning("\nmodel has no variability, pROC will return NA.\n")
#     p_roc <- rep(NA, 2)
#     names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", threshold,
#                              "%"), "pval_pROC")
#     auc_ratios <- rep(NA, 3)
#     names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC",
#                            "AUC_ratio")
#     p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
#   } else {
#     classpixels <- classpixels %>% dplyr::mutate(value = rev(value),
#                                                  count = rev(count), totpixperclass = cumsum(count),
#                                                  percentpixels = totpixperclass/sum(count)) %>% dplyr::arrange(value)
#     error_sens <- 1 - (threshold/100)
#     prediction_errors <- classpixels[, "value"]
#     fractional_area <- classpixels[, "percentpixels"]
#     n_data <- length(pred_test)
#     n_samp <- ceiling((rand_percent/100) * n_data)
#     big_classpixels <- matrix(rep(prediction_errors, each = n_samp),
#                               ncol = length(prediction_errors))
#     partial_AUC <- 1:iterations %>% purrr::map_df(~kuenm:::calc_aucDF(big_classpixels,
#                                                                       fractional_area, pred_test, n_data, n_samp, error_sens))
#     naID <- !is.na(partial_AUC$auc_ratio)
#     nona_valproc <- partial_AUC$auc_ratio[naID]
#     mauc <- mean(nona_valproc)
#     proc <- sum(nona_valproc <= 1)/length(nona_valproc)
#     p_roc <- c(mauc, proc)
#     names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", threshold,
#                              "%"), "pval_pROC")
#     auc_ratios <- partial_AUC
#     names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC",
#                            "AUC_ratio")
#     p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
#   }
# }
# #t <- proc_maxnet(pred_test = suit_val_eval, model = pred_i)

