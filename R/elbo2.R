elbo2=function(V, W, B, Y, Ngene, Nsample, Ncell){

  V=array(V,dim=c(Nsample,Ngene,Ncell))
  W=matrix(W, nrow=Ngene, ncol=Ncell)
  B=matrix(B, nrow=Nsample, ncol=Ncell)

  # s0=c()
  # for(i in 1:Ngene){
  #   s0[i]=var(log(Y[,i]))
  # }
  s0 <- matrixStats::colVars(log(Y))
  s <- t(replicate(Nsample, s0))

  B0=rowSums(B)
  Bt=sweep(B,1,B0,'/')

  t1=exp(sweep(V,2:3,(W^2)/2,'+'))
  eq=rowSums(sweep(t1,c(1,3),Bt,'*'),dims=2)

  t21=exp(sweep(2*V,2:3,2*(W^2),'+'))
  t22=sweep(Bt*(1-Bt),1,B0+1,'/')+Bt^2
  t20=sweep(t21,c(1,3),t22,'*')
  t23=exp(sweep(2*V,2:3,(W^2),'+'))
  t200=sweep(t23,c(1,3),(Bt^2),'*')
  t2=rowSums(t20-t200,dims=2)

  t3_list_i <- list()
  for(i in 1:Ncell){
    t3_list_j <- list()
    for(j in c(1:Ncell)[-i]){
      V1=V[,,-i]
      V2=V[,,-j]
      W1=W[,-i]
      W2=W[,-j]
      Bt1=Bt[,-i]
      Bt2=Bt[,-j]
      t31=exp(sweep(V1+V2,2:3,(W1^2+W2^2)/2,'+'))
      t32=sweep(-1*Bt1*Bt2,1,(B0+1),'/')
      t33=sweep(t31,c(1,3),t32,'*')
      t3_list_j[[j]] <- rowSums(t33, dims=2)
    }
    t3_list_i[[i]] <- do.call(sum, t3_list_j)
  }
  t3 <- do.call(sum, t3_list_i)

  vq <- t2+t3

  Eq <- ((Y-log(eq)) - vq / (2*eq^2))^2
  Vq <- vq/(eq^2)

  e21 <- Eq + Vq
  e22 <- 2*(s^2)

  e2 <- sum(-1 * e21 / e22)

  return(e2)
}
