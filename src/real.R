# assuming it's run from the experiment folder
source("./src/datautil.R") # load data utilities
source("./src/selection.R")
source("./src/measurements.R") # load metrics
source("./src/plots.R")
source("./src/plotutils.R")
source("./src/cd.R")
source("./src/ols.R")

plot_datasets <- function(d, de_train_idx, do_train_idx, d_test_idx) {

  d_ <- d
  d_$g <- "l"
  d_[de_train_idx,]$g <- paste0("In-Sample (e) (", sum(de_train_idx), " observations)")
  d_[do_train_idx,]$g <- paste0("In-Sample (o) (", sum(do_train_idx), " observations)")
  d_[d_test_idx,]$g <- paste0("Out-of-Sample (", sum(d_test_idx), " observations)")
  p <- ggplot(d_) + geom_density(aes(y, group = g, fill = g), alpha = 0.5)
  plot(p)

}

select <- function(d, idx, names) {
  d_ <- d[idx,]

  X <- as.matrix(d_[, (colnames(d_) %in% names)])
  y <- t(t(d_$y))

  return(list(X = X, y = y, names = names))
}

center_train_data <- function(data) {

  data_ <- data
  mean_ <- (mean(data_$ye) + mean(data_$yo)) / 2
  data_$ye <- data_$ye - mean_
  data_$yo <- data_$yo - mean_

  meanX_ <- NULL

  return(list(data = data_, mean = mean_, meanX = meanX_))

}


# compute
resample_train <- function(d, names,
                           de_train_idx, do_train_idx, d_test_idx,
                           method,
                           compute_estimator,
                           folds,
                           resamples,
                           debug = FALSE,
                           centerX = FALSE) {

  plot_datasets(d,
                de_train_idx = de_train_idx,
                do_train_idx = do_train_idx,
                d_test_idx = d_test_idx)

  rols <- vector(mode = "numeric", length = resamples)
  rcd <- vector(mode = "numeric", length = resamples)
  rcr <- vector(mode = "numeric", length = resamples)
  min_ <- vector(mode = "numeric", length = resamples)
  plots <- vector(mode = "list", length = resamples)

  d_test <- select(d, d_test_idx, names)

  for (i in 1:resamples) {
    writeLines(paste("Resample", i))

    do_train_idx_B <- sample(which(do_train_idx), replace = TRUE)
    do_train <- select(d, do_train_idx_B, names)

    de_train_idx_B <- sample(which(de_train_idx), replace = TRUE)
    de_train <- select(d, de_train_idx_B, names)

    data <- list(Xe = de_train$X,
                 ye = de_train$y,
                 Xo = do_train$X,
                 yo = do_train$y)

    centering <- center_train_data(data, centerX = centerX)
    data <- centering$data

    # choose gamma by cross-validation
    estimands <- cross_validation(folds = folds, data = data, estimators = c(method))
    l <- estimands[[1]]$opt

    # plot cv loss and coefficient path
    # plot(estimands[[1]]$plot_cv)
    plots[[i]] <- estimands[[1]]$plot
    plot(plots[[i]])

    # s <- method(data)
    # p <- plot_path(gamma = s$gamma, beta = s$beta, names = de_train$names)
    # p <- p + geom_vline(xintercept = l, color = "red", linetype = "dotted")
    # plot(p)

    m <- moments(data)
    cr <- compute_estimator(m, l = l)
    ols <- compute_ols(m)
    cd <- compute_cd(m)

    # test set
    # print(centering$mean)
    # d_test$y <- d_test$y - centering$mean

    rcr[[i]] <- risk(d_test$X, d_test$y, beta = cr)
    rcd[[i]] <- risk(d_test$X, d_test$y, beta = cd)
    rols[[i]] <- risk(d_test$X, d_test$y, beta = ols)

    # smallest mean square error per experiment
    risks <- c(rcr[[i]], rcd[[i]])
    order_ <- order(risks, decreasing = FALSE)
    min_[[i]] <- order_[1]

  }

  return(list(min = min_,
              plots = plots,
              cr = rcr,
              ols = rols,
              cd = rcd))

}

resample_test <- function(d, names,
                          de_train_idx, do_train_idx, d_test_idx,
                          method,
                          compute_estimator,
                          folds,
                          resamples,
                          debug = FALSE) {

  de_train <- select(d, de_train_idx, names)
  do_train <- select(d, do_train_idx, names)

  data <- list(Xe = de_train$X, ye = de_train$y, Xo = do_train$X, yo = do_train$y)

  centering <- center_train_data(data)
  data <- centering$data

  estimands <- cross_validation(folds = folds, data = data, estimators = c(method))
  l <- estimands[[1]]$opt

  writeLines(paste("Choosen regularization", l))

  m <- moments(data)
  cr <- compute_estimator(m = m, l = l)
  cd <- compute_cd(m = m)

  if (debug) {
    # plot cv loss and coefficient path
    plot(estimands[[1]]$plot)

    s <- method(data)
    p <- plot_path(gamma = s$gamma, beta = s$beta, names = de_train$names)
    p <- p + geom_vline(xintercept = l, color = "red", linetype = "dotted")
    p <- p + geom_hline(yintercept = cd, linetype = "dashed", color = "black")
    p <- p + geom_hline(yintercept = ols, linetype = "dashed", color = "black")
    plot(p)

    writeLines(paste("OLS", ols))
    writeLines(paste("CR", cr))
    writeLines(paste("CD", cd))
  }

  # test sets
  rcd <- vector(mode = "numeric", length = resamples)
  rcr <- vector(mode = "numeric", length = resamples)

  for (i in 1:resamples) {

    idx <- sample(which(d_test_idx), replace = TRUE)

    d_test_ <- select(d, idx, names)
    d_test_$y <- d_test_$y - centering$mean

    rcr[[i]] <- risk(d_test_$X, d_test_$y, beta = cr)
    rcd[[i]] <- risk(d_test_$X, d_test_$y, beta = cd)

  }

  return(list(cv = estimands[[1]],
              cr = rcr,
              cd = rcd))

}