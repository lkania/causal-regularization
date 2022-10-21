####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/data.R")
source("./src/datautil.R") # load data utilities
source("./src/measurements.R") # load metrics
source("./src/plotutils.R")
path <- './experiments/out_of_sample_risk/'

library(doParallel)
registerDoParallel(cores = 12)

out_of_sample <- function(target, B,
                          beta,
                          oshifts, # out-of-sample shifts
                          on = 1e6,
                          normalize = TRUE) { # out-of-sample sample size

  nbeta <- dim(beta)[1]
  noshifts <- length(oshifts)

  os <- matrix(0, nrow = noshifts, ncol = nbeta)

  for (j in 1:noshifts) {

    shift <- oshifts[[j]]

    odata <- gen(target = target, B = B, n = on, sd = shift) # out-of-sample data

    for (i in 1:nbeta) {

      b <- t(t(beta[i,]))

      os[j, i] <- risk(X = odata$X, y = odata$y, beta = b)

      if (normalize) {
        os[j, i] <- os[j, i] / shift
      }

    }

  }

  return(os)

}

homogenize <- function(x, n, shifts, gamma, fun) {

  summary_ <- array(unlist(x), c(dim(x[[1]]), length(x)))
  summary_ <- t(apply(summary_, 1:2, FUN = fun))

  summary_ <- data.frame(summary_)
  names(summary_) <- shifts
  summary_$gamma <- gamma
  summary_ <- melt(summary_,
                   id.vars = 'gamma',
                   value.name = 'risk',
                   variable.name = 'shift')
  summary_$n <- n

  return(summary_)

}

do_trial <- function(target, B, n, init_shift,
                     oshifts, compute_estimator,
                     gamma, normalize = TRUE) {
  # in-sample fit
  data <- data_gen(target = target, B = B, n = n, sd = init_shift)

  # out-of-sample performance
  return(out_of_sample(target = target,
                       B = B,
                       beta = compute_estimator(data = data, gamma = gamma)$beta,
                       oshifts = oshifts, normalize = normalize))
}

out_of_sample_risk <- function(target,
                               B,
                               ns,
                               oshifts,
                               init_shift,
                               trials,
                               compute_estimator,
                               fun,
                               gamma, normalize = TRUE) {

  os_per_n <- vector("list", length(ns))

  for (k in 1:length(ns)) {

    n <- ns[[k]]
    writeLines(paste("Sample size", n))

    oss <- foreach(i = 1:trials) %dopar% { do_trial(target = target,
                                                    B = B,
                                                    n = n,
                                                    gamma = gamma,
                                                    compute_estimator = compute_estimator,
                                                    init_shift = init_shift,
                                                    oshifts = oshifts,
                                                    normalize = normalize) }

    os_per_n[[k]] <- oss

    writeLines(paste("n =", n, "done"))

  }

  ffun <- (function(x, n) homogenize(x = x, n = n,
                                     shifts = oshifts,
                                     gamma = gamma,
                                     fun = fun))
  summary <- mapply(ffun, os_per_n, ns, SIMPLIFY = FALSE)

  summary <- data.table::rbindlist(summary, use.names = TRUE)

  return(summary)

}


target <- 4
B <- matrix(c(0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0,
              1, 1, 0, 0, 0, 0, 0,
              0, 1, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0), nrow = 7, byrow = TRUE)

source("./src/scr.R")
gamma_ <- c(seq(0, 0.4, 0.01), seq(0.45, 1, 0.1), 1)
compute_estimator <- scr

set.seed(123)
data_median <- out_of_sample_risk(target = target,
                                  B = B,
                                  ns = c(250, 1000),
                                  oshifts = c(10, 100, 1000),
                                  init_shift = 1,
                                  trials = 100,
                                  compute_estimator = compute_estimator,
                                  gamma = gamma_,
                                  fun = median)
save(data_median, file = paste0(path, 'data_median.Rda'))

produce_plot <- function(data_) {
  main_text_size <- 12
  pdata <- data_
  pdata$n <- as.factor(pdata$n)
  pdata$shift <- as.factor(pdata$shift)
  levels(pdata$n) <- paste("Sample size =", levels(pdata$n), sep = " ")
  p <- ggplot(pdata, aes(x = gamma, y = risk, color = shift, linetype = shift)) +
    geom_line(lwd = 1) +
    facet_wrap(. ~ n, ncol = 2, scales = 'free') +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = expression(hat(beta)[CD],
                                           '0.25',
                                           '0.5',
                                           '0.75',
                                           hat(beta)[OLS])) +
    ylab('Median normalized out-of-sample risk') +
    xlab(expression('Regularization path' ~ hat(beta)[gamma])) +
    theme(legend.position = c(0.85, 0.8)) +
    labs(color = "Out-of-sample shift", linetype = "Out-of-sample shift") +
    theme(legend.background =
            element_rect(colour = 'black', fill = 'white', linetype = 'solid')) +
    theme(strip.text.x = element_text(size = main_text_size),
          legend.text = element_text(size = main_text_size),
          legend.title = element_text(size = main_text_size),
          axis.text = element_text(size = main_text_size),
          axis.title = element_text(size = main_text_size)) +
    guides(color = guide_legend(nrow = 1))
  return(p)
}

p <- produce_plot(data_median)
save_plot(p,
          filename = paste0(path, 'out_of_sample_risk.pdf'),
          dims = list(height = 3.5, width = 8))