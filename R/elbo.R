elbo <- function(X, Y, N0, N00, N000, mu, a, a0, Ngene, Nsample, Ncell){
  V <- X[1:N0]
  W <- X[(N0 + 1):(N0 + N00)]
  B <- X[(N0 + N00 + 1):(N0 + N00 + N000)]

  a1 <- elbo1(V, W, mu, a0, Ngene, Ncell, Nsample)
  a2 <- elbo2(V, W, B, Y, Ngene, Nsample, Ncell)
  a3 <- elbo3(B, a, Nsample, Ncell)
  a4 <- elbo4(W, Ngene, Ncell)
  a5 <- elbo5(B, Nsample, Ncell)

  e <- a1 + a2 + a3 - a4 - a5

  return(e)
}
