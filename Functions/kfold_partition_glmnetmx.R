#### Get k-folds ####
kfold_partition_glmnetmx <- function(data, n_folds = 4, pr_bg, seed = 42) {
  pre <- which(data[, pr_bg] == 1)
  aus <- which(data[, pr_bg] == 0)
  set.seed(seed)
  foldp <- sample(cut(seq(1, length(pre)), breaks = n_folds, labels = FALSE))
  folda <- sample(cut(seq(1, length(aus)), breaks = n_folds, labels = FALSE))
  #Create columns to save pr_bg
  data$folds <- NA
  data$folds[which(data[,pr_bg] == 1)] <- foldp
  data$folds[which(data[,pr_bg] == 0)] <- folda
  return(data)
}
