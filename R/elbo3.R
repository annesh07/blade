elbo3=function(B){
  B=matrix(B, nrow=Nsample, ncol=Ncell)

  f=(a-1)*sweep(digamma(B),1,digamma(rowSums(B)),"-")
  f0=rowSums(f)
  f00=(apply(gamma(a),1,prod))/gamma(rowSums(a))
  e3=sum(-log(f00)+f0)
  return(e3)
}
