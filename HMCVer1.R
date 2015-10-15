library(MASS)
HMC1= function(TotalSamples, density, q,M=diag(length(q)),epsilon,L, densitydiff="NaN", burnin) {
  oldq=q; #Initialisation of q
  d=length(q); #Dimensions
  N=TotalSamples+burnin; #How many samples will be needed
  p=mvnorm(N,rep(0,d),M); 
  #Generate multivariate normal, this is an N by d quantity. I generate all at once here.
  for (k in 1:N){
    oldp=p[k,]; #Take one momentum vector.
    for (t in 1:L){
      ptand
    }
  }
}