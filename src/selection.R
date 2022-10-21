###########################################################################
# Fit estimators using all data / sample splitting / cross-validation
###########################################################################

fit_cross_validation_fixed_range <- function(folds, f, data_train, data_val) {

  # Assumes that the range of gamma is always the same

  measurements_per_fold <- vector("list", folds)
  max_measure <- 0

  for (i in 1:folds) {

    fsol <- f(data = data_train[[i]])
    gamma <- fsol$gamma

    measurements_per_fold[[i]] <- measure(
      beta = fsol$beta,
      data = data_val[[i]],
      metric = abs_difference)

    max_measure <- max(max_measure, max(measurements_per_fold[[i]]))

  }

  avg_measurement <- do.call(rbind, measurements_per_fold)
  avg_measurement <- colMeans(avg_measurement)

  idx_selected <- which.min(avg_measurement)
  gamma_selected <- gamma[idx_selected]

  ########################################
  # plot cross-validation selection
  ########################################

  # build data-frame
  df <- data.frame(x = gamma, y = avg_measurement,
                   fold = 0, type = 'Average')
  for (i in 1:folds) {
    df <- rbind(df, data.frame(
      x = gamma,
      y = measurements_per_fold[[i]],
      fold = i, type = 'Fold'))
  }
  df <- rbind(df, data.frame(
    x = c(gamma_selected, gamma_selected),
    y = c(min(df$y), max(df$y)),
    fold = -1, type = 'Selected'))

  p <- ggplot(df) +
    geom_line(mapping =
                aes(x = x, y = y, group = fold, color = type, linetype = type)) +
    scale_color_manual(values = c('black', 'black', 'red')) +
    labs(color = "", linetype = "") +
    theme(legend.title = element_blank())

  ticks <- c(0, 0.25, 0.5, 0.75, 1, gamma_selected)
  label_ticks <- expression(
    hat(beta)[CD], '0.25', '0.5', '0.75', hat(beta)[OLS], hat(beta)[CR])
  idx <- order(ticks)

  p <- p +
    scale_x_continuous(breaks = ticks[idx], labels = label_ticks[idx]) +
    ylab('Absolute risk difference') +
    xlab(expression('Regularization path' ~ hat(beta)[gamma]))

  ########################################

  return(list(opt = gamma_selected, plot = p))
}

fit_cross_validation <- function(
  folds, f, data_train, data_val) {

  return(fit_cross_validation_fixed_range(
    folds = folds,
    f = f,
    data_train = data_train,
    data_val = data_val
  ))

}

subset_data <- function(data, sele = NULL, selo = NULL) {

  Xo <- data$Xo
  yo <- data$yo

  Xe <- data$Xe
  ye <- data$ye
  
  if (!is.null(selo)) {
    Xo <- as.matrix(Xo[selo,])
    yo <- as.matrix(yo[selo,])
  }

  if (!is.null(sele)) {
    Xe <- as.matrix(Xe[sele,])
    ye <- as.matrix(ye[sele,])
  }

  return(list(Xe = Xe, ye = ye, Xo = Xo, yo = yo))

}

split_idx <- function(folds, n) {

  #equal sized splits
  g <- sample(rep(1:folds, floor(n / folds)))

  res <- n %% folds
  if (res > 0) {
    g <- c(g, sample(1:folds, res))
  }

  return(g)

}

cross_validation <- function(folds, data, estimators) {
  # The functions assumes that all provided estimators have fixed range

  go <- split_idx(folds, n = dim(data$Xo)[1])
  ge <- split_idx(folds, n = dim(data$Xe)[1])

  data_train <- vector("list", folds)
  data_val <- vector("list", folds)

  for (i in 1:folds) {

    data_train[[i]] <- subset_data(data,
                                   sele = (ge != i),
                                   selo = (go != i))

    data_val[[i]] <- subset_data(data,
                                 sele = (ge == i),
                                 selo = (go == i))

  }

  estimands <- vector("list", length(estimators))

  for (k in 1:length(estimators)) {

    estimands[[k]] <- fit_cross_validation(folds = folds,
                                           f = estimators[[k]],
                                           data_train = data_train,
                                           data_val = data_val)

  }

  return(estimands)
}

