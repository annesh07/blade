grad3=function(V, W, B, Y, a, Ngene, Ncell, Nsample){
  V <- array(V,dim=c(Ngene,Ncell,Nsample))
  W <- matrix(W, nrow=Ngene, ncol=Ncell)
  B <- matrix(B, nrow=Nsample, ncol=Ncell)

  s0 <- matrixStats::colVars(log(Y))
  s <- t(replicate(Nsample, s0))

  B0 <- rowSums(B)
  Bt <- sweep(B,1,B0,'/')
  Bt2 <- sweep(B,1,B0^2,'/')

  t1 <- exp(sweep(V,2:3,(W^2)/2,'+'))
  eq <- rowSums(sweep(t1,c(1,3),Bt,'*'),dims=2)
  deq <- rowSums(sweep(t1,c(1,3),(1/B0),'*'),dims=2)-rowSums(sweep(t1,c(1,3),Bt2,'*'),dims=2)

  t21 <- exp(sweep(2*V,2:3,2*(W^2),'+'))
  t22 <- sweep(Bt*(1-Bt),1,B0+1,'/')+Bt^2
  t20 <- sweep(t21,c(1,3),t22,'*')
  t23 <- exp(sweep(2*V,2:3,(W^2),'+'))
  t200 <- sweep(t23,c(1,3),(Bt^2),'*')
  t2 <- rowSums(t20-t200,dims=2)

  B01 <- sweep(B,1,(B0^2),'/')
  B02 <- sweep(-2*B,1,B0,'+')
  B03 <- sweep(B02,1,((B0^2)*(B0+1)),'/')+2*B01
  dt201 <- sweep(t21,c(1,3),B03,'*')
  dt2001 <- sweep(t23,c(1,3),(2*B01),'*')
  dt21 <- rowSums(dt201-dt2001,dims=2)

  B04 <- sweep(B^2,1,(B0^3),'/')
  B05 <- sweep(-1*B,1,B0,'+')*B
  B06 <- sweep(B,1,((B0^2)*(B0+1)),'/')-sweep(B05,1,(3*B0+2)/((B0^3)*(B0+1)^2),'*')-2*B04
  dt202 <- sweep(t21,c(1,3),B06,'*')
  dt2002 <- sweep(t23,c(1,3),(2*B04),'*')
  dt22 <- rowSums(dt201+dt2001,dims=2)

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

  dt31_list_i <- list()
  for(i in 1:Ncell){
    dt31_list_j <- list()
    for(j in c(1:Ncell)[-i]){
      V1=V[,,-i]
      V2=V[,,-j]
      W1=W[,-i]
      W2=W[,-j]
      B1=B[,-i]
      B2=B[,-j]
      t31=exp(sweep(V1+V2,2:3,(W1^2+W2^2)/2,'+'))
      t32=sweep(B1*B2,1,(3*B0+2)/((B0^3)*(B0+1)^2),'*')
      t33=sweep(t31,c(1,3),t32,'*')
      dt31_list_j[[j]] <- rowSums(t33, dims=2)
    }
    dt31_list_i[[i]] <- do.call(sum, dt31_list_j)
  }
  dt31 <- do.call(sum, dt31_list_i)

  dt32_list_i <- list()
  for(i in 1:Ncell){
    dt32_list_j <- list()
    for(j in c(1:Ncell)[-i]){
      V1=V[,,-i]
      V2=V[,,-j]
      W1=W[,-i]
      W2=W[,-j]
      B1=B[,-i]
      B2=B[,-j]
      t31=exp(sweep(V1+V2,2:3,(W1^2+W2^2)/2,'+'))
      t32=sweep(B1,1,((B0^2)*(B0+1)),'/')
      t33=sweep(t31,c(1,3),t32,'*')
      dt32_list_j[[j]] <- rowSums(t33, dims=2)
    }
    dt32_list_i[[i]] <- do.call(sum, dt32_list_j)
  }
  dt32 <- do.call(sum, dt32_list_i)

  vq <- t2 + t3
  dvq <- dt21 + dt22 + dt31 - dt32

  da <- (dvq*(eq^2)-2*deq*vq)/(eq^3)
  db <- -(Y-log(eq)-vq/(2*eq^2))*(2*deq/eq+da)

  e22 <- 2*(s^2)

  d2 <- sum((-1)*(da+db)/e22)

  f <- (a-1)*trigamma(B)-sweep((a-1),1,trigamma(rowSums(B)),'*')
  f0 <- rowSums(f)
  d3 <- sum(f0)

  g <- (B-1)*trigamma(B)-sweep((B-1),1,trigamma(rowSums(B)),'*')
  g0 <- rowSums(g)
  d5 <- sum(g0)

  db <- d2+d3-d5
  return(db)
}
