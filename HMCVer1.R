library(MASS)
HMCVer1 <- function(TotalSamples, densityA,M, q,epsilon,L, densitydiff, burnin) {
  oldq <-qtandepsilon <- q #Initialisation of q and qtandepsilon
  d <- length(q) #Dimensions
  if(length(M)>1){MassVector <- colSums(M)} else MassVector <- M
  N <- TotalSamples+burnin; #How many samples will be needed
  Samples <- matrix(0,N,d);
  #For allocation of memory
  p <- mvrnorm(N,rep(0,d),diag(d)); 
  #Generate multivariate normal, this is an N by d quantity. I generate all at once here.
  for (k in 1:N){
    oldp <-ptplusepsilon <- p[k,] 
    #Take one momentum vector.
    for (t in 1:L){
      #LeapFrog L times with step size epsilon
      ptplusepsilonhalf= ptplusepsilon - (epsilon/2)*densitydiff(qtandepsilon)
      qtandepsilon =  qtandepsilon + ((epsilon/MassVector)*ptplusepsilonhalf)
      #Assume independent in p vector hence use column sums
      ptplusepsilon= ptplusepsilonhalf - (epsilon/2)*densitydiff(qtandepsilon)
    }
    OldUenergy= -log(densityA(oldq));
    #Energy of target densityA at old state
    OldKenergy= sum((oldp^2/MassVector))/2;
    #Energy of momentum at Old state
    OldHamiltonian= OldUenergy+OldKenergy;
    #Old Hamiltonian
    ProposeUenergy= -log(densityA(qtandepsilon));
    ProposeKenergy=sum((ptplusepsilon^2/MassVector))/2;
    #return(c(ProposeKenergy))
    #Propose Energy
    ProposeHamiltonian=ProposeUenergy+ProposeKenergy;
  
    #Propose Hamiltonian
    if (runif(1)< exp(OldHamiltonian-ProposeHamiltonian)) {
      oldp=ptplusepsilon;
      oldq=qtandepsilon;
    #Metropolis Hastings, if accepted move to new state
    }
    Samples[k]=oldq
    #Samples
    #forget Intialisation
  }
  return(Samples[(burnin+1:N),])
}


JJ <- HMCVer1(TotalSamples = 10000,densityA = dnorm,M=1,q = 0,epsilon = 0.05,L = 20,densitydiff = function(x) x,burnin = 0)
hist(JJ,freq=F)
lines(seq(-4,4,0.01),dnorm(seq(-4,4,0.01)))

# Multivariate example - needs work
JJ <- HMCVer1(TotalSamples = 1000,densityA = function(l) dmvnorm(l,c(0,0),matrix(c(1,0.7,0.7,1),2,2)),q = c(0,0),M=diag(length(q)),epsilon = 0.05,L = 60,densitydiff = function(x) t(t(x)%*%solve(matrix(c(1,0.7,0.7,1),2,2))), burnin = 0)



