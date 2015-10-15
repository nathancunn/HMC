library(MASS)
HMCVer1= function(TotalSamples, densityA , q,M,epsilon,L, densitydiff, burnin) {
  oldq=q; #Initialisation of q
  d=length(q); #Dimensions
  #MassVector=colSums(M); #Mass Of particles
  #MassVector=1
  if(length(M)>1){MassVector <- colSums(M)} else MassVector <- M
  N=TotalSamples+burnin; #How many samples will be needed
  Samples=matrix(0,N,2*d);
  #For allocation of memory
  p=mvrnorm(N,rep(0,d),M); 
  #Generate multivariate normal, this is an N by d quantity. I generate all at once here.
  Uni=runif(N);
  for (k in 1:N){
    oldp=p[k,]; #Take one momentum vector.
    for (t in 1:L){
      #LeapFrog L times with step size epsilon
      ptplusepsilonhalf= oldp- (epsilon/2)*densitydiff(oldq);
      qtandepsilon =  oldq + ((epsilon/MassVector)*ptplusepsilonhalf);
      #Assume independent in p vector hence use column sums
      ptplusepsilon= ptplusepsilonhalf -(epsilon/2)*densitydiff(qtandepsilon);
    }
    OldUenergy= -log(densityA(oldq));
    #Energy of target densityA at old state
    OldKenergy= (t(oldp)/(2*MassVector))%*%((oldp));
    #Energy of momentum at Old state
    OldHamiltonian= OldUenergy+OldKenergy;
    #Old Hamiltonian
    ProposeUenergy= -log(densityA(qtandepsilon));
    ProposeKenergy=(((ptplusepsilon))/(2*MassVector))%*%t(ptplusepsilon);
    #return(c(ProposeKenergy))
    #Propose Energy
    ProposeHamiltonian=ProposeUenergy+ProposeKenergy;
  
    #Propose Hamiltonian
    if (Uni[k]<=exp(-OldHamiltonian+ProposeHamiltonian)) {
      oldp=ptplusepsilon;
      oldq=qtandepsilon;
    #Metropolis Hastings, if accepted move to new state
    }
    Samples[k,]=rbind(oldp,oldq)
    #Samples
    #forget Intialisation
  }
  return(Samples[(burnin+1:N),])
}


JJ <- HMCVer1(TotalSamples = 10000,densityA = dnorm,q = 0,M=1,epsilon = 0.2,L = 20,densitydiff = function(x) x,burnin = 0)
TotalSamples = 100
density = dmvnorm
q=c(0,0)
M=diag(length(q))
epsilon=0.2
L=20
densitydiff = function(x) x
burnin=0



plot(JJ[,1])
lag.plot(JJ[,1])
lag(JJ[,1],1)
JJ[,1]
hist(JJ[,1])

plot(JJ[,1],type="b")



JJ <- HMCVer1(TotalSamples = 100,densityA = dmvnorm,q = c(0,0),M=diag(length(q)),epsilon = 0.2,L = 20,densitydiff = function(x) x,burnin = 0)



JJ <- HMCVer1(TotalSamples = 1000,densityA = function(l) dmvnorm(l,c(0,0),matrix(c(1,0.7,0.7,1),2,2)),q = c(0,0),M=diag(length(q)),epsilon = 0.05,L = 60,densitydiff = function(x) t(t(x)%*%solve(matrix(c(1,0.7,0.7,1),2,2))), burnin = 0)



