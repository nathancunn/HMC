library(MASS)
HMC1= function(TotalSamples, density, q,M=diag(length(q)),epsilon,L, densitydiff="NaN", burnin) {
  oldq=q; #Initialisation of q
  d=length(q); #Dimensions
  MassVector=colSums(M); #Mass Of particles
  N=TotalSamples+burnin; #How many samples will be needed
  Samples=matrix(0,N,2*d);
  #For allocation of memory
  p=mvnorm(N,rep(0,d),M); 
  #Generate multivariate normal, this is an N by d quantity. I generate all at once here.
  H=runif(N);
  for (k in 1:N){
    oldp=p[k,]; #Take one momentum vector.
    for (t in 1:L){
      #LeapFrog L times with step size epsilon
      ptplusepsilonhalf= oldp- (epsilon/2)*densitydiff(oldq);
      qtandepsilon =  oldq + ((epsilon/MassVector)*ptplusepsilonhalf);
      #Assume independent in p vector hence use column sums
      ptplusepsilon= ptplusepsilonhalf -(epsilon/2)*densitydiff(qtandepsilon);
    }
    OldUenergy= -log(density(oldq));
    #Energy of target density at old state
    OldKenergy= (oldp/MassVector)%*%(t(oldp));
    #Energy of momentum at Old state
    OldHamiltonian= OldUenergy+OldKenergy;
    #Old Hamiltonian
    ProposeUenergy= -log(density(qtandepsilon));
    ProposeKenergy=(ptplusepsilon/(2*MassVector))%*%(t(ptplusepsilon));
    #Propose Energy
    ProposeHamiltonian=ProposeUenergy+ProposeKenergy;
    #Propose Hamiltonian
    if (H[k]<=exp(-OldHamiltonian+ProposeHamiltonian)) {
      oldp=ptplusepsilon;
      oldq=qtandepsilon;
    #Metropolis Hastings, if accepted move to new state
    }
    Samples[k,]=rbind(oldp,oldq)
    #Samples
    HMC1=Samples[(burnin+1:N),]
    #forget Intialisation
  }
}