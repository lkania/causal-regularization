###########################################################################
# Data generation
###########################################################################

data_gen <- function(target, B, n, mean = 0, sd = 0, confounding = TRUE) {

  do <- gen(target, B, n, mean = mean, sd = 0, confounding = confounding)
  de <- gen(target, B, n, mean = mean, sd = sd, confounding = confounding)
  return(list(Xe = de$X, ye = de$y, Xo = do$X, yo = do$y))

}

gen <- function(target, B, n, mean = 0, sd = 0, confounding = TRUE) {
  p <- dim(B)[1]

  shift <- 0
  if (sd > 0) {

    shift <- sqrt(sd) * matrix(rnorm(p * n, mean = mean, sd = 1), nrow = p)
    shift[target,] <- 0

  }

  conf <- 0
  if (confounding) {
    conf <- t(matrix(1, n, p) * rnorm(n, 0, 1))
  }

  noise <- matrix(rnorm(p * n, 0, 1), nrow = p)

  all <- t(solve(diag(p) - B, shift + noise + conf))
  X <- all[, -target]
  y <- as.matrix(all[, target])
  return(list(X = X, y = y))
}
