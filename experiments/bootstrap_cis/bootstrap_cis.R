library(doParallel)
registerDoParallel(cores = 12)
library(boot)

####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/data.R") # load data generation methods
source("./src/datautil.R") # load data utilities
source("./src/measurements.R") # load metrics
source("./src/utils.R")
source("./src/scr.R")
path <- './experiments/bootstrap_cis/'

bootrisk <- function(data, idxs, l) {
  sdata <- data[idxs,]
  ps <- dim(sdata)[2]
  p <- ps / 2 - 1
  Xe <- sdata[, 1:p]
  ye <- t(t(sdata[, p + 1]))
  Xo <- sdata[, (p + 2):(ps - 1)]
  yo <- t(t(sdata[, ps]))

  datab <- list(Xe = Xe, ye = ye, Xo = Xo, yo = yo)

  m <- moments(datab)
  b <- compute_scr(m = m, l = l)
  return(difference(datab, beta = b))

}

bootstrap_cis <- function(target,
                          B,
                          init_shift,
                          resamples,
                          trials,
                          l,
                          ns,
                          alpha) {

  ##################################################
  # population risk
  ##################################################

  # population quantities
  betapa <- t(t(B[target, -target]))
  betach <- t(t(B[-target, target]))
  Bx <- B[-target, -target]
  p <- dim(Bx)[1]
  Ip <- diag(p)

  iM <- Ip - Bx - betach %*% t(betapa)
  M <- solve(iM)
  h <- t(M) %*% betapa
  m <- h + solve(betach %*% t(betach) + (3 / 2) * Ip, betach)

  # beta computation
  CX_A <- M %*% (betach %*% t(betach) + 2 * Ip) %*% t(M)
  CX_0 <- M %*% (betach %*% t(betach) + Ip) %*% t(M)

  EXY_A <- CX_A %*% betapa + M %*% betach
  EXY_0 <- CX_0 %*% betapa + M %*% betach

  weight <- 2 * l - 1
  beta <- solve(weight * CX_0 + CX_A, weight * EXY_0 + EXY_A)

  # risk computation
  b <- t(M) %*% beta
  population_risk_difference <- sum((b - h)^2)

  writeLines(paste('population risk diff',
                   population_risk_difference))

  ##################################################
  # checking that the population risk 
  # difference computation was correct
  # using large n=1M
  ##################################################

  pop_data <- data_gen(target = target, B = B,
                       n = 1e6, sd = init_shift,
                       confounding = FALSE)

  numerical_risk_difference <- abs_difference(data = pop_data,
                                              beta = compute_scr(
                                                m = moments(pop_data),
                                                l = l))
  writeLines(paste('numerical error',
                   abs(population_risk_difference - numerical_risk_difference)))

  ##################################################
  # conf. intervals
  ##################################################

  stats <- vector("list", length = length(ns))
  width <- vector("list", trials)
  coverage <- vector("list", trials)
  point_estimates <- vector("list", trials)

  k <- 1
  for (n in ns) {

    for (trial in 1:trials) {

      ##################################################
      # data generation
      ##################################################

      data <- data_gen(target, B, n, sd = init_shift, confounding = FALSE)

      ##################################################
      # confidence interval via bootstrap 
      # for \hat{beta}_\lambda
      ##################################################

      boot_ <- boot(cbind(data$Xe, data$ye, data$Xo, data$yo),
                    (function(data, idx) bootrisk(data = data, idxs = idx, l = l)),
                    R = resamples,
                    parallel = "multicore",
                    ncpus = 10)

      ##################################################
      # pivotal boostrap confidence interval
      ##################################################
      q <- as.numeric(unlist(
        quantile(unlist(boot_$t, use.names = FALSE), c(alpha / 2, 1 - alpha / 2))))
      point_estimate <- boot_$t0
      high <- 2 * point_estimate - q[1]
      low <- 2 * point_estimate - q[2]

      # since we know that the parameter of interest is non-negative
      low <- max(0, low)

      # compute statistics
      point_estimates[[trial]] <- point_estimate
      width[[trial]] <- high - low
      below_high <- (high >= population_risk_difference)
      above_low <- (population_risk_difference >= low)
      coverage[[trial]] <- below_high && above_low

    }

    width <- unlist(width)
    coverage <- unlist(coverage)
    point_estimates <- unlist(point_estimates)

    average_point_estimate <- mean(point_estimates)
    average_width <- mean(width)
    median_width <- median(width)
    average_coverage <- mean(coverage)

    stats[[k]] <- list(average_width = average_width,
                       median_width = median_width,
                       average_coverage = average_coverage,
                       average_point_estimate = average_point_estimate)
    k <- k + 1

    writeLines(paste('n =', n,
                     'avg_coverage = ', round(average_coverage, 3),
                     'avg_width = ', round(average_width, 3),
                     'avg_pe = ', round(average_point_estimate, 6),
                     sep = '\t'))

  }

  stats <- do.call(rbind, stats)
  stats <- data.frame(stats)
  stats$n <- ns
  stats <- sapply(stats, as.numeric)
  stats <- as.data.frame(stats)

  return(stats)
}

target <- 4
B <- matrix(c(0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0,
              1, 1, 0, 0, 0, 0, 0,
              0, 1, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0), nrow = 7, byrow = TRUE)

set.seed(123)
s <- bootstrap_cis(
  target = target,
  B = B,
  init_shift = 1,
  resamples = 1000,
  trials = 500,
  l = 0.3,
  ns = seq(100, 2000, 100),
  alpha = 0.05)

save(s, file = paste0(path, 'bootstrap_cis.Rda'))

source("./src/plotutils.R") # load plot utilities

main_text_size <- 12
sample_size_axis <- c(200, 500, 1000, 1500, 2000)
names(s) <- c("Average", "Median", "coverage", "point_estimate", "n")
melted <- melt(subset(s[-1,], select = c("Average", "n")), 'n')
p1 <- ggplot(melted, aes(x = n, y = value)) +
  geom_line() +
  xlab('Sample size') +
  ylab('') +
  ggtitle('Average width') +
  scale_x_continuous(breaks = sample_size_axis) +
  theme(strip.text.x = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size),
        legend.title = element_text(size = main_text_size),
        axis.text = element_text(size = main_text_size),
        axis.title = element_text(size = main_text_size))
melted <- subset(s[-1,], select = c(n, coverage))
p2 <- ggplot(melted, aes(x = n, y = coverage)) +
  geom_line() +
  xlab('Sample size') +
  ylab('') +
  ggtitle('Empirical coverage') +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  geom_hline(yintercept = min(s$coverage), linetype = "dashed", color = "black", alpha = 0.3) +
  geom_hline(yintercept = max(s$coverage), linetype = "dashed", color = "black", alpha = 0.3) +
  scale_y_continuous(breaks = c(0.9,
                                round(min(s$coverage), 3),
                                0.95,
                                round(max(s$coverage), 3),
                                1),
                     limits = c(0.9, 1)) +
  scale_x_continuous(breaks = sample_size_axis) +
  theme(strip.text.x = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size),
        legend.title = element_text(size = main_text_size),
        axis.text = element_text(size = main_text_size),
        axis.title = element_text(size = main_text_size))
p <- ggarrange(p1, p2, ncol = 2)
save_plot(p,
          filename = paste0(path, 'bootstrap_ci.pdf'),
          dims = list(height = 4, width = 8))
