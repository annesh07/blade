#you might want to assign the location of the working directory to that folder where the codes will be downloaded
setwd("D:/Bordeaux/Codes") 
library('rBeta2009')

Nsample=20
Ngene=200
Ncell=3
Noise=0.5

F0=rdirichlet(Nsample,rep(5,Ncell))
mu=matrix(rnorm(Ngene*Ncell,mean=0,sd=2),nrow=Ngene)
sigma=matrix(Noise,nrow=Ngene,ncol=Ncell)
x=array(rnorm(Nsample,mean=mu,sd=sigma),dim=c(Ngene,Ncell,Nsample))

Y=matrix(0,nrow=Nsample,ncol=Ngene)
for(i in 1:Nsample){
  Y[i,]=exp(x[,,i])%*%F0[i,]
}

a=matrix(5, nrow=Nsample, ncol=Ncell)
a0=matrix(0.1, nrow=Ngene, ncol=Ncell)
k0=matrix(1, nrow=Ngene, ncol=Ncell)

source("blade.R")
blade(Y,mu,sigma,a,k0,a0)
