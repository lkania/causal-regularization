quant <- function(x,q=c(0.05,0.95)){
  
  return(quantile(unlist(x,use.names = FALSE),q))
  
}