blade=function(Y, mu, sigma, a, k0, a0){
  Nsample=dim(Y)[1]
  Ngene=dim(Y)[2]
  Ncell=dim(mu)[2]
  #mu0=mu
  #b0=a0*sigma
  #s0=c()
  #for(i in 1:Ngene){
  #  s0[i]=var(log(Y[,i]))
  #}
  #s=t(replicate(Nsample,s0))
  #s=matrix(1,nrow=Nsample,ncol=Ngene)
  
  source("elbo1.R")
  source("elbo2.R")
  source("elbo3.R")
  source("elbo4.R")
  source("elbo5.R")
  
  N0=Nsample*Ngene*Ncell
  N00=Ngene*Ncell
  N000=Nsample*Ncell
  
  elbo=function(X){
    V=X[1:N0]
    W=X[(N0+1):(N0+N00)]
    B=X[(N0+N00+1):(N0+N00+N000)]
    a1=elbo1(V,W)
    a2=elbo2(V,W,B)
    a3=elbo3(B)
    a4=elbo4(W)
    a5=elbo5(B)
    e=a1+a2+a3-a4-a5
    return(e)
  }
  r=optim(c(array(0,dim=c(Nsample,Ngene,Ncell)),matrix(0,Ngene,Ncell),matrix(0.01,Nsample,Ncell)),elbo,method="L-BFGS-B",lower=c(array(0.01,dim=c(Nsample,Ngene,Ncell)),matrix(0.01,Ngene,Ncell),matrix(0.01,Nsample,Ncell),upper=c(array(Inf,dim=c(Nsample,Ngene,Ncell)),matrix(Inf,Ngene,Ncell),matrix(1,Nsample,Ncell))))
}
