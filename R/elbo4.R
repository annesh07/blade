elbo4=function(W, Ngene, Ncell){
  W=matrix(W, nrow=Ngene, ncol=Ncell)

  e4=(-Nsample/2)*sum(log(2*pi*(W)^2)-1)
  return(e4)
}
