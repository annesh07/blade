#' BLADE main function wrapper
#'
#' call elbo1 to elbo5 functions to perform variational inference
#'
#' @param Y description to add
#' @param mu description to add
#' @param a description to add
#' @param k0 description to add
#' @param a0 description to add
#'
#' @export
#'
blade <- function(Y, mu, sigma, a, k0, a0){

  Nsample <- dim(Y)[1]
  Ngene <- dim(Y)[2]
  Ncell <- dim(mu)[2]
  #mu0=mu
  #b0=a0*sigma
  #s0=c()
  #for(i in 1:Ngene){
  #  s0[i]=var(log(Y[,i]))
  #}
  #s=t(replicate(Nsample,s0))
  #s=matrix(1,nrow=Nsample,ncol=Ngene)

  N0 <- Nsample*Ngene*Ncell
  N00 <- Ngene*Ncell
  N000 <- Nsample*Ncell

  r <- optim(par = c(array(0, dim = c(Nsample, Ngene, Ncell)),
            matrix(0, Ngene, Ncell), matrix(0.01, Nsample, Ncell)),
          fn = elbo, gr = grad,
          Y = Y, N0 = N0, N00 = N00, N000 = N000, mu=mu, a0=a0,
          Ngene = Ngene, Nsample=Nsample, Ncell= Ncell, a=a,
          method = "L-BFGS-B",
          lower = c(array(0.01, dim = c(Nsample, Ngene, Ncell)),
                    matrix(0.01, Ngene, Ncell),
                    matrix(0.01, Nsample, Ncell)),
          upper = c(array(Inf, dim = c(Nsample, Ngene, Ncell)),
                    matrix(Inf, Ngene, Ncell), matrix(1, Nsample,Ncell)),
          control = list("trace"=1, "REPORT"=1, maxit = 2)
  )

  return(r)
}
