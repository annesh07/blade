elbo5=function(B, Nsample, Ncell){
  B=matrix(B, nrow=Nsample, ncol=Ncell)

  f=(B-1)*sweep(digamma(B),1,digamma(rowSums(B)),"-")
  f0=rowSums(f)
  f00=(apply(gamma(B),1,prod))/gamma(rowSums(B))
  e5=sum(-log(f00)+f0)
  return(e5)
}
