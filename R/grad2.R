grad2=function(V, W, B, Y, mu, a0, Ngene, Ncell, Nsample){
  V <- array(V,dim=c(Ngene,Ncell,Nsample))
  W <- matrix(W, nrow=Ngene, ncol=Ncell)
  mu0 <- mu
  b0 <- a0*sigma

  Vjt <- rowSums(V,dims=2)/Nsample
  Vijt <- sweep(V,1:2,Vjt,"-")
  V0 <- rowSums((Vijt^2),dims=2)
  e1 <- (b0+((Nsample-1)/2)*W^2+0.5*V0+(Nsample*k0)*((W^2)/2+(Vjt-mu0)^2)/(2*(k0+Nsample)))

  e11 <- (Nsample-1)*W+(k0*W)/(k0+Nsample)

  d1 <- sum(-(a0+Nsample/2)*e11/e1)

  V=array(V,dim=c(Nsample,Ngene,Ncell))
  B=matrix(B, nrow=Nsample, ncol=Ncell)

  s0 <- matrixStats::colVars(log(Y))
  s <- t(replicate(Nsample, s0))

  B0 <- rowSums(B)
  Bt <- sweep(B,1,B0,'/')

  t1 <- exp(sweep(V,2:3,(W^2)/2,'+'))
  eq <- rowSums(sweep(t1,c(1,3),Bt,'*'),dims=2)
  eq1 <- sweep(t1,c(1,3),Bt,'*')
  deq <- rowSums(sweep(eq1,2:3,W,"*"),dims=2)

  t21 <- exp(sweep(2*V,2:3,2*(W^2),'+'))
  t22 <- sweep(Bt*(1-Bt),1,B0+1,'/')+Bt^2
  t20 <- sweep(t21,c(1,3),t22,'*')
  dt20 <- 4*sweep(t20,2:3,W,"*")
  t23 <- exp(sweep(2*V,2:3,(W^2),'+'))
  t200 <- sweep(t23,c(1,3),(Bt^2),'*')
  dt200 <- 2*sweep(t200,2:3,W,"*")
  t2 <- rowSums(t20-t200,dims=2)
  dt2 <- rowSums(dt20-dt200,dims=2)

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
  vq <- t2 + t3

  dt3_list_i <- list()
  for(i in 1:Ncell){
    dt3_list_j <- list()
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
      t333 <- sweep(t33,2:3,W1,"*")
      dt3_list_j[[j]] <- rowSums(t333, dims=2)
    }
    dt3_list_i[[i]] <- do.call(sum, dt3_list_j)
  }

  dt3 <- do.call(sum, dt3_list_i)
  dvq <- dt2+dt3

  da <- (dvq*(eq^2)-2*deq*vq)/(eq^3)
  db <- -(Y-log(eq)-vq/(2*eq^2))*(2*deq/eq+da)

  e22 <- 2*(s^2)

  d2 <- sum((-1)*(da+db)/e22)

  d4 <- sum((-1)*(Nsample)/W)

  dw <- d1+d2-d4
  return(dw)
}
