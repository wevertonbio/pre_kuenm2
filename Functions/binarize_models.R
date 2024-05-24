binarize_models <- function(candidate_results, pred,
                            metric_consensus = "median",
                            omission_rate = 10) {
  m <- pred[[metric_consensus]]
  xy <- candidate_results$data_xy[which(
    candidate_results$calibration_data$pr_bg == 1),]
  v <- terra::extract(m, xy, ID = FALSE)[[1]]
  thr <- sort(v)[round(nrow(xy) * (omission_rate/100)) + 1]
  m_bin <- app(m, function(x) ifelse(x >= thr, 1, 0))
  return(m_bin)
}

# #Example
# library(terra)
# candidate_results <- readRDS("candidate_results.RDS")
# pred <- rast("Predictions/General_consensus.tiff")
# bn <- binarize_models(candidate_results, pred, omission_rate = 10)
# plot(bn)
# bn5 <- binarize_models(candidate_results, pred, omission_rate = 5)
# plot(bn5)
