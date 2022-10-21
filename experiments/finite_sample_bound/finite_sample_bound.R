####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/data.R") # load data generation methods
source("./src/datautil.R") # load data utilities
source("./src/measurements.R") # load metrics
source("./src/plotutils.R")
path <- './experiments/finite_sample_bound/'

finite_sample_bound <- function(target, B, ns, trials) {

  #############################################
  # Assumptions: 
  # 1. The noise and shift variables are centered
  # 2. The covariance of the noise variable is the identity matrix
  # 3. The covariance of the shift variable for the covariates, i.e. A_X^e,
  #    is the identity matrix
  #############################################
  init_shift <- 1

  #############################################
  # Compute population quantities
  #############################################

  # population quantities
  betapa <- t(t(B[target, -target]))
  betach <- t(t(B[-target, target]))
  Bx <- B[-target, -target]
  p <- dim(Bx)[1]
  Ip <- diag(p)

  iM <- Ip - Bx - betach %*% t(betapa)
  M <- solve(iM)

  # betas computation
  CX_A <- M %*% (betach %*% t(betach) + 2 * Ip) %*% t(M)
  CX_0 <- M %*% (betach %*% t(betach) + Ip) %*% t(M)

  EXY_A <- CX_A %*% betapa + M %*% betach
  EXY_0 <- CX_0 %*% betapa + M %*% betach

  beta <- betapa

  G_delta <- CX_A - CX_0
  G_plus <- CX_0 + CX_A
  betaols <- solve(G_plus, EXY_0 + EXY_A)

  # risk computation
  population_pool_risk <- t(betaols - beta) %*% solve(G_plus, betaols - beta)
  population_risk_difference <- t(betapa - beta) %*% solve(G_delta, betapa - beta)
  population_worst_risk <- 0.5 * population_pool_risk + (0.5 + init_shift) * population_risk_difference

  writeLines(paste('population risk difference', population_risk_difference))
  writeLines(paste('population pool risk', population_pool_risk))
  writeLines(paste('population worst risk', population_worst_risk))

  t <- data.frame(matrix(0, nrow = length(ns) + 1, ncol = 6))
  names(t) <- c('n',
                'sample risk difference',
                'erd',
                'sample pool risk',
                'epd',
                'ner')

  idx <- length(ns) + 1
  t[idx, 1] <- 8
  t[idx, 2] <- population_risk_difference
  t[idx, 3] <- 0
  t[idx, 4] <- population_pool_risk
  t[idx, 5] <- 0
  t[idx, 6] <- 0

  #############################################
  # Empirical averages
  #############################################

  normalization_constant <- (1 + init_shift) * (sum(abs(beta))^2 + 1)

  k <- 1
  for (n in ns) {

    avg_sample_risk_difference <- 0
    avg_sample_pool_risk <- 0

    for (trial in 1:trials) {

      data <- data_gen(target = target, B = B, n = n, sd = init_shift, confounding = FALSE)

      sample_risk_difference <- abs_difference(data = data, beta = beta)
      avg_sample_risk_difference <- avg_sample_risk_difference + sample_risk_difference

      m <- moments(data)
      empirical_betaols <- solve(m$Gplus, m$Zplus)
      sample_pool_risk <- t(empirical_betaols - beta) %*% solve(m$Gplus, empirical_betaols - beta)
      avg_sample_pool_risk <- avg_sample_pool_risk + sample_pool_risk

    }

    avg_sample_risk_difference <- avg_sample_risk_difference / trials
    avg_excess_risk_difference <- (avg_sample_risk_difference - population_risk_difference) / normalization_constant

    avg_sample_pool_risk <- avg_sample_pool_risk / trials
    avg_excess_pool_risk <- (avg_sample_pool_risk - population_pool_risk) / normalization_constant

    avg_normalized_excess_sup_risk <- 0.5 * avg_excess_pool_risk + (0.5 + init_shift) * avg_excess_risk_difference

    writeLines(paste('n =', n))

    writeLines(paste('\tsample risk difference', avg_sample_risk_difference))
    writeLines(paste('\texcess risk difference', avg_excess_risk_difference, '\n'))
    writeLines(paste('\tsample pool risk', avg_sample_pool_risk))
    writeLines(paste('\texcess pool risk', avg_excess_pool_risk, '\n'))
    writeLines(paste('\tnormalized excess risk', avg_normalized_excess_sup_risk))

    t[k, 1] <- n
    t[k, 2] <- avg_sample_risk_difference
    t[k, 3] <- avg_excess_risk_difference
    t[k, 4] <- avg_sample_pool_risk
    t[k, 5] <- avg_excess_pool_risk
    t[k, 6] <- avg_normalized_excess_sup_risk

    k <- k + 1
  }

  return(t)

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
ns <- c(1e1, 1e2, 1e3, 1e4, 1e5)
result <- finite_sample_bound(target = target,
                              B = B,
                              ns = ns,
                              trials = 1000)

save(result, file = paste0(path, 'result.Rda'))

# plot

data <- result[-(length(ns) + 1),]

empirical <- abs(data[, 6])
predicted <- 1 / sqrt(ns)
df <- rbind(data.frame(x = data[, 1], y = empirical, type = 'empirical'),
            data.frame(x = data[, 1], y = predicted, type = 'predicted'))
main_text_size <- 12
scientific10 <- function(x) parse(text = paste0("10^", as.character(log10(x))))
p <- ggplot(df, aes(x = x, y = y, linetype = type)) +
  geom_line() +
  scale_linetype_discrete(labels = c(
    'predicted' = expression(paste('Predicted ', O[p](1 / sqrt(n)))),
    'empirical' = "Empirical average")) +
  theme(legend.title = element_blank()) +
  theme(aspect.ratio = 1,
        legend.text.align = 0,
        legend.margin = margin(t = 0, r = 0.05, b = 0, l = 0.05, unit = "inch")) +
  xlab('Sample size = n (log scale)') +
  ylab('Absolute excess risk (log scale)') +
  scale_y_continuous(trans = 'log10', breaks = 1 / sqrt(c(ns, 1e6)), labels = scientific10) +
  scale_x_continuous(trans = 'log10', breaks = ns, labels = scientific10) +
  theme(axis.text = element_text(size = main_text_size),
        axis.title = element_text(size = main_text_size)) +
  theme(legend.position = c(0.75, 0.9)) +
  theme(legend.background =
          element_rect(colour = 'black', fill = 'white', linetype = 'solid'))
plot(p)
save_plot(p,
          filename = paste0(path, 'finite-sample-bound.pdf'),
          dims = list(width = 4, height = 4))


