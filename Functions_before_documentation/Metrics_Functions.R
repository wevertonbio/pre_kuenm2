####Ver como funciona categorical####

####Functions to fit hinge features####
hinge <- function(x, nknots=50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length=nknots)
  lh <- outer(x, utils::head(k,-1), function(w,h) hingeval(w, h, max))
  rh <- outer(x, k[-1], function(w,h) hingeval(w, min, h))
  colnames(rh) <- paste("", min, k[-1], sep=":")
  cbind(lh, rh)
}

hingeval <- function (x, min, max)
{
  pmin(1, pmax(0, (x - min)/(max - min)))
}

####Functions to fit thresholds features####
thresholds <- function (x, nknots = 50) {
  min <- min(x)
  max <- max(x)
  k <- seq(min, max, length = nknots + 2)[2:nknots + 1]
  f <- outer(x, k, function(w, t) ifelse(w >= t, 1, 0))
  colnames(f) <- paste("", k, sep = ":")
  f
}

thresholdval <- function (x, knot) {
  ifelse(x >= knot, 1, 0)
}

####Functions to fit categorical features####
categorical <- function (x) {
  f <- outer(x, levels(x), function(w, f) ifelse(w == f, 1,
                                                 0))
  colnames(f) <- paste("", levels(x), sep = ":")
  f
}
categoricalval <- function (x, category) {
  ifelse(x == category, 1, 0)
}



#### Get k-folds ####
kfold_part <- function(data, n_folds = 4, pr_bg) {
  pre <- which(data[, pr_bg] == 1)
  aus <- which(data[, pr_bg] == 0)
  foldp <- sample(cut(seq(1, length(pre)), breaks = n_folds, labels = FALSE))
  folda <- sample(cut(seq(1, length(aus)), breaks = n_folds, labels = FALSE))
  #Create columns to save pr_bg
  data$folds <- NA
  data$folds[which(data[,pr_bg] == 1)] <- foldp
  data$folds[which(data[,pr_bg] == 0)] <- folda
  return(data)
}

####Metrics functions####

####Omission Rate####
omrat_maxnet <- function(threshold = 5, pred_train, pred_test) {
  om_rate <- vector("numeric", length = length(threshold))
  for (i in 1:length(threshold)) {
    val <- ceiling(length(pred_train) * threshold[i]/100) + 1
    omi_val_suit <- sort(pred_train)[val]
    om_rate[i] <- as.numeric(length(pred_test[pred_test <
                                                    omi_val_suit])/length(pred_test))
  }
  names(om_rate) <- paste("Omission_rate_at_", threshold, sep = "")
  return(om_rate)
}
# t2 <- omrat_maxnet(threshold = 5,
#                    pred_train =  suit_val_cal,
#                    pred_test = suit_val_eval)


####AIC####
aic_maxnet <- function(exp_pred_occ, ncoefs, pred_aic = NULL) {
  #Get sum of whole area
  p.sum <- sum(values(pred_aic, na.rm=TRUE))
  exp_pred_occ <- exp_pred_occ/p.sum

  #Calculate AIC following ENMeval::aic.maxent()
  n.occs <- length(exp_pred_occ)
  AIC.valid <- ncoefs < n.occs
  for (i in 1:length(AIC.valid)) {
    if (AIC.valid[i] == FALSE) {
      message(paste("Warning: model", names(AIC.valid)[i],
                    "has more non-zero coefficients (ncoef) than occurrence records for training, so AIC cannot be calculated."))
    }
  }
  LL <- sum(log(exp_pred_occ), na.rm = TRUE)
  AICc <- (2 * ncoefs - 2 * LL) + (2 * (ncoefs) * (ncoefs +
                                                     1)/(n.occs - ncoefs - 1))
  AICc <- sapply(1:length(AICc), function(x)
    ifelse(AIC.valid[x] ==FALSE | is.infinite(AICc[x]), NA, AICc[x]))
  return(AICc)
  }


####Get formulas####
get_formulas_maxnet <- function (independent, type = "lqpht",
                                 categorical_var = NULL,
                                 minvar = 1, maxvar = NULL)
{
  if (is.character(type)) {
    if (!all(unlist(strsplit(type, "")) %in% c("l", "p",
                                               "q", "h", "t"))) {
      stop("'type' must be: 'l', 'p', 'q', 'h', 't', or a combination of those three.")
    }
  }  else {
    stop("'type' must be of class character.")
  }

  #Check if categorical variable is in independente variables. If not, set null
  if(!is.null(categorical_var)) {
    if(!(categorical_var %in% independent)) {
      categorical_var <- NULL
    }  }

  if(!is.null(categorical_var)) {
    independent <- setdiff(independent, categorical_var)
  }

  predictors <- independent
  npred <- length(predictors)
  aux <- " "
  #Linear
  if (grepl("l", type)) {
    aux <- paste(aux, paste(predictors, collapse = " + "),
                 sep = " + ")
  }
  #Quadratic
  if (grepl("q", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("I(", predictors[v], "^2)"),
                   sep = " + ")
    }
  }
  #Product
  if (grepl("p", type)) {
    if (npred > 1) {
      inter_tab <- utils::combn(predictors, 2)
      aux_inter <- paste0(" + ", paste(apply(inter_tab,
                                             2, paste, collapse = ":"), collapse = " + "))
      if (length(aux_inter) > 0) {
        aux <- paste0(aux, aux_inter)
      }
    }
    else {
      if (grepl("l", type) | grepl("q", type)) {
        message("'p' is is only possible with 2 or more independent variables.",
                "\nReturning other combinations.")
      }
      else {
        stop("'p' is is only possible with 2 or more independent variables.",
             "\nTry other combinations of type.")
      }
    }
  }
  #Hinge
  if (grepl("h", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("hinge(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  #Threshold
  if (grepl("t", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("thresholds(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  out <- paste0("~",
                gsub(pattern = "  \\+ ", x = aux, replacement = ""))

  if(!is.null(categorical_var)) {
    out <- paste0(out, " + categorical(", categorical_var, ")")
  }
  return(out)
}
#my_formula <- get_formulas(independent = c("PC1", "PC2"), type = "lq")

#### New functions for maxnet ####
glmnet_mx <- function(p, data, f, regmult = 1.0,
                      regfun = maxnet.default.regularization,
                      addsamplestobackground = TRUE, weights_1_0 = c(1, 100),
                      calculate_AIC = FALSE, ...) {
  if (anyNA(data)) {
    stop("NA values in data table. Please remove them and rerun.")
  }
  if (!is.vector(p)) {
    stop("p must be a vector.")
  }

  if (addsamplestobackground) {
    pdata <- data[p == 1, ]
    ndata <- data[p == 0, ]

    # add to background any presence data that isn't there already
    wadd <- !do.call(paste, pdata) %in% do.call(paste, ndata)
    p <- c(p, rep(0, sum(wadd)))
    toadd <- pdata[wadd, ]
    data <- rbind(data, toadd)
  }

  mm <- model.matrix(f, data)
  reg <- regfun(p, mm) * regmult
  weights <- ifelse(p == 1, weights_1_0[1], weights_1_0[2])
  lambdas <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) *
    sum(p) / sum(weights)

  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",
                          standardize = FALSE, penalty.factor = reg,
                          lambda = lambdas, weights = weights, ...)

  bb <- coef(model)[, 200]

  class(model) <- c("glmnet_mx", class(model))
  if (length(model$lambda) < 200) {
    msg <- "glmnet failed to complete regularization path. Model may be infeasible."
    if (!addsamplestobackground) {
      msg <- paste(msg, "\n\tTry re-running with addsamplestobackground = TRUE.")
    }
    stop(msg)
  }

  # AIC calculation
  filter <- bb[-1] != 0
  bb <- c(bb[1], bb[-1][filter])

  if (calculate_AIC & sum(filter) != 0) {
    model$AIC <- aic_glmnetmx(x = as.matrix(mm[, filter]), y = p, beta = bb)
  } else {
    model$AIC <- NA
  }

  # returning other elements
  model$betas <- bb[-1]
  model$alpha <- 0
  rr <- predict.glmnet_mx(model, data[p == 0, , drop = FALSE],
                          type = "exponent", clamp = FALSE)
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw))
  model$alpha <- -log(sum(rr))

  model$penalty.factor <- reg
  model$featuremins <- apply(mm, 2, min)
  model$featuremaxs <- apply(mm, 2, max)

  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2, max)

  means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data)[!vv], function(n) {
    which.max(table(data[p == 1, n, drop = FALSE]))
  })
  names(majorities) <- names(data)[!vv]
  model$samplemeans <- unlist(c(means, majorities))

  model$levels <- lapply(data, levels)

  return(model)
}



#### maxnet default regularization ####
maxnet.default.regularization <- function (p, m) {
  isproduct <- function(x) {grepl(":", x) & !grepl("\\(", x)}
  isquadratic <- function(x) {grepl("^I\\(.*\\^2\\)", x)}
  ishinge <- function(x) {grepl("^hinge\\(", x)}
  isthreshold <- function(x) {grepl("^thresholds\\(", x)}
  iscategorical <- function(x) {grepl("^categorical\\(", x)}

  regtable <- function(name, default) {
    if (ishinge(name)) {
      return(list(c(0, 1), c(0.5, 0.5)))
    }
    if (iscategorical(name)) {
      return(list(c(0, 10, 17), c(0.65, 0.5, 0.25)))
    }
    if (isthreshold(name)) {
      return(list(c(0, 100), c(2, 1)))
    }
    default
  }

  lregtable <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
  qregtable <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5, 0.25, 0.05))
  pregtable <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9, 0.55, 0.05))

  mm <- m[p == 1, ]
  np <- nrow(mm)
  lqpreg <- lregtable

  if (sum(isquadratic(colnames(mm)))) {
    lqpreg <- qregtable
  }
  if (sum(isproduct(colnames(mm)))) {
    lqpreg <- pregtable
  }

  classregularization <- sapply(colnames(mm), function(n) {
    t <- regtable(n, lqpreg)
    approx(t[[1]], t[[2]], np, rule = 2)$y
  })
  classregularization <- classregularization / sqrt(np)


  ishinge <- grepl("^hinge\\(", colnames(mm))

  hmindev <- sapply(1:ncol(mm), function(i) {
    if (!ishinge[i]) {
      return(0)
    }
    avg <- mean(mm[, i])
    std <- max(sd(mm[, i]), 1 / sqrt(np))
    std * 0.5 / sqrt(np)
  })

  tmindev <- sapply(1:ncol(mm), function(i) {
    ifelse(isthreshold(colnames(mm)[i]) && (sum(mm[, i]) ==
                                              0 || sum(mm[, i]) == nrow(mm)), 1, 0)
  })

  pmax(0.001 * (apply(m, 2, max) - apply(m, 2, min)), hmindev,
       tmindev, apply(as.matrix(mm), 2, sd) * classregularization)
}


#### aic for glmnet like maxent ####
aic_glmnetmx <- function(x, y, beta) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.")
  }
  if (mode(x) != "numeric") {
    stop("x must be numeric.")
  }

  x_glm <- glm(y ~ x, family = "binomial")
  predict_glm <- predict(x_glm)

  pr0 <- 1 - 1 / (1 + exp(predict_glm))
  pr <- 1 - 1 / (1 + exp(cbind(1, x) %*% beta))
  cc <- as.vector(beta)
  jc <- abs(cc) > 1e-05
  xj <- cbind(1, x)[, jc]
  Pi <- diag(as.vector(pr * (1 - pr)))
  Pi0 <- diag(as.vector(pr0 * (1 - pr0)))
  j22 <- t(xj) %*% Pi %*% xj
  i22 <- t(xj) %*% Pi0 %*% xj

  aic <- -2 * sum(y * log(pr) + (1 - y) * log(1 - pr)) +
    2 * sum(diag(solve(j22) %*% i22))

  return(aic)
}


#### predict gmlm like maxent to new data ####
predict.glmnet_mx <- function (object, newdata, clamp = FALSE,
                               type = c("link", "exponential",
                                        "cloglog", "logistic")) {
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v]),
                           object$varmax[v])
    }
  }

  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
               names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
               terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
               terms)
  f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))

  mm <- model.matrix(f, data.frame(newdata))

  if (clamp) {
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                 object$featuremaxs[names(object$betas)]))
  }

  link <- (mm %*% object$betas) + object$alpha
  type <- match.arg(type)

  # return prediction
  if (type == "link") {
    return(link)
  }
  if (type == "exponential") {
    return(exp(link))
  }
  if (type == "cloglog") {
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic") {
    return(1/(1 + exp(-object$entropy - link)))
  }
}

####Summarize results####
#Function to summarize data
eval_stats <- function(calib_results, toagg = NULL, agg_by = NULL,
                       to_keep = NULL){
  if(is.null(toagg)){
  toagg <- c("Omission_rate_at_5", "Omission_rate_at_10",
             "proc_auc_ratio", "proc_pval") }
  if(is.null(agg_by)){
    agg_by <- c("Formulas", "regm", "Features")
  }
  if(is.null(to_keep)){
    to_keep <- c("ID", "Formulas", "regm", "Features", "AIC",
                 "AIC_Warren", "npar", "is_concave")
  }
  agg_formula <- paste("~", paste(agg_by, collapse = " + "))

  #Get summary
  xy <- lapply(toagg, function(x) {
    do.call(data.frame, stats::aggregate(as.formula(paste(x, agg_formula)),
                                  data = calib_results, FUN = function(y) c(mean = round(mean(y), 4), sd = round(sd(y), 4)), na.action=NULL))
  })

  #Summarize stats
  stats <- Reduce(function(x, y) merge(x, y,
                                       by = agg_by),
                  xy)

  stats_AICS <- calib_results[!duplicated(calib_results[,to_keep]),][,to_keep]
  stats_final <- merge(stats, stats_AICS, by = agg_by)
  return(stats_final)
}


####Create empty dataframes####
empty_replicates <- function(omrat_thr = c(5, 10),
                             n_row = 4, replicates = 1:4,
                             is_c = NA) {
  column_names <- c("Replicate", paste0("Omission_rate_at_", omrat_thr),
                    "proc_auc_ratio", "proc_pval", "AIC", "AIC_Warren", "npar", "is_concave")
  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names
  df_eval_q$Replicate <- replicates
  df_eval_q$is_concave = is_c
  return(df_eval_q)
}

empty_summary <- function(omrat_thr, is_c){
  om_means <- paste0("Omission_rate_at_", omrat_thr, ".mean")
  om_sd <- paste0("Omission_rate_at_", omrat_thr, ".sd")
  column_names <- c(om_means, om_sd,
                    "proc_auc_ratio.mean", "proc_auc_ratio.sd", "proc_pval.mean",
                    "proc_pval.sd", "AIC", "AIC_Warren", "npar", "is_concave")
  eval_final_q  <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$is_concave = is_c
  return(eval_final_q)
}

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

####Get grid of formulas####
create_grid <- function(var_names = NULL, swd = NULL, x = NULL, y = NULL,
                        fold_column = "folds",
                        min.number = 2,
                        categorical_var = NULL,
                        features = c("l", "q", "lq", "lqp", "p"),
                        min_continuous = NULL,
                        regm = c(0.1, 1, 2, 3, 5)){
  if(is.null(var_names) & !is.null(swd)) {
    ignore_columns <- intersect(c(x, y, "pr_bg", "species", fold_column),
                                colnames(swd))
    var_names <- colnames(swd[, !names(swd) %in% ignore_columns])}

  #Get var combinations
  var_comb <- kuenm::all_var_comb(var.names = var_names,
                                  min.number = min.number)

  #Split features
  formula_x <- unlist(lapply(seq_along(features), function(x){
    f_x <- features[x]
    if(grepl("p", f_x) & !is.null(categorical_var)) {
      var_comb_new <- var_comb[sapply(var_comb,
                                      function(x)
                                        sum(!grepl(categorical_var, x))) >= 2]
    } else {
      var_comb_new <- var_comb}
    ff_x <- unlist(lapply(seq_along(var_comb_new), function(i){
      #If type = p, get only combinations with 2 or more combinations

      f_l <- get_formulas_maxnet(independent = var_comb_new[[i]], type = f_x,
                                 categorical_var = categorical_var)
      names(f_l) <- f_x
      return(f_l) }))
  }))
  #Create dataframe with formulas
  formula_d <- data.frame(formula = paste(formula_x, -1),
                          features = names(formula_x))
  #Expand grid
  f_grid <- unique(expand.grid(formula_d$formula, regm))
  colnames(f_grid) <- c("Formulas", "regm")
  f_grid <- merge(f_grid, formula_d, by.x = "Formulas", by.y = "formula")
  f_grid$ID <- 1:nrow(f_grid)
  colnames(f_grid) <- c("Formulas", "regm", "Features", "ID")
  #Reorder columns
  f_grid <- f_grid[, c("ID", "Formulas", "regm", "Features")]
  #Class formulas
  f_grid$Formulas <- as.character(f_grid$Formulas)
  return(f_grid)
}


# #Test
# my_grid <- create_grid(var_names = c("bio_1", "bio_7", "bio_12", "bio_15", "soil"),
#                       min.number = 2, categorical_var = "soil",
#                       features =  c("l", "q", "lq", "lqp", "p"),
#                       regm = c(0.1, 1, 2, 3, 5))



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

