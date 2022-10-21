###########################################################################
# Causal Dantzig
###########################################################################

library(CVXR)

compute_cd <- function(m) {

  beta <- Variable(rows = m$p)
  objective <- Minimize(sum_squares(m$G %*% beta - m$Z))
  problem <- Problem(objective)
  sol <- solve(problem)

  return(sol$getValue(beta))

}
