elbo1=function(V,W){
  V=array(V,dim=c(Ngene,Ncell,Nsample))
  W=matrix(W, nrow=Ngene, ncol=Ncell)
  mu0=mu
  b0=a0*sigma
  
  Vjt=rowSums(V,dims=2)/Nsample
  Vijt=sweep(V,1:2,Vjt,"-")
  V0=rowSums((Vijt^2),dims=2)
  e1=-sum((a0+Nsample/2)*log(b0+((Nsample-1)/2)*W^2+0.5*V0+(Nsample*k0)*((W^2)/2+(Vjt-mu0)^2)/(2*(k0+Nsample))))
  return(e1)
}