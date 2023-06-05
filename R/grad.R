grad <- function(X, Y, N0, N00, N000, mu, a0, a, Ngene, Nsample, Ncell){
  V <- X[1:N0]
  W <- X[(N0 + 1):(N0 + N00)]
  B <- X[(N0 + N00 + 1):(N0 + N00 + N000)]

  g1 <- grad1(V, W, B, mu, a0, Ngene, Ncell, Nsample)
  g2 <- grad2(V, W, B, Y, mu, a0, Ngene, Ncell, Nsample)
  g3 <- grad3(V, W, B, Y, a, Ngene, Ncell, Nsample)

  g <- c(g1,g2,g3)

  return(e)
}
