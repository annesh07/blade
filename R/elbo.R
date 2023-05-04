elbo <- function(X, Y, N0, N00, N000){
  V <- X[1:N0]
  W <- X[(N0 + 1):(N0 + N00)]
  B <- X[(N0 + N00 + 1):(N0 + N00 + N000)]

  a1 <- elbo1(V, W)
  a2 <- elbo2(V, W, B, Y)
  a3 <- elbo3(B)
  a4 <- elbo4(W)
  a5 <- elbo5(B)

  e <- a1 + a2 + a3 - a4 - a5

  return(e)
}
