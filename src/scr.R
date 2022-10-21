###########################################################################
# Squared Causal Regularization (analytic implementation)
###########################################################################

# assuming it's run from the experiments folder
source("./src/datautil.R") # load data utilities

compute_scr <- function(m, l) {
  weight <- 2 * l - 1
  Gl <- weight * m$XXo + m$XXe
  Zl <- weight * m$XYo + m$XYe
  return(limSolve::Solve(A = Gl, B = Zl, tol = 1e-6))

}

scr <- function(data, length = 200, gamma = NULL) {

  m <- moments(data)

  if (is.null(gamma)) {
    gamma <- seq(0, 1, length = length)
  }

  beta <- NULL

  for (l in gamma) {

    beta <- rbind(beta, t(compute_scr(m, l)))

  }

  return(list(gamma = gamma, beta = beta, name = "SCR"))

}