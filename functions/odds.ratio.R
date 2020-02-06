# Odds ratio
odds.ratio <- function(a, b, total){
  a.b <- length(intersect(a, b))
  a <- length(a) 
  b <- length(b) 
  a.nonb <- a - a.b
  nona.b <- b - a.b
  nona.nonb <- total - a.b - a.nonb - nona.b
  v <- c(a.b, a.nonb, nona.b, nona.nonb)
  x <- matrix(v, 2, byrow = TRUE)
  or <- OddsRatio(x) # Odds of b
  or
}