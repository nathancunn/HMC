libs <- list("MASS","mvtnorm")
lapply(libs,library,character.only = T)
HMC <- function(TotalSamples, densityA, M, q, epsilon, L, densitydiff, burnin) {
  # Arguments
  # TotalSamples - the number of simulations you'd like output
  # densityA
  # M - the mass matrix
  # q - starting values
  # epsilon
  # L - the number of leapfrog iterations
  # densitydiff - the derivative of the distribution
  # burnin - the burn-in
  epsilonFix <- epsilon
  rejections <- 0
  oldq <-q.epsilon <- q # Initialisation of q and q.epsilon
  d <- length(q) # Dimensions
  # The elements of p are assumed independent, hence use column sums
  if(length(M)>1){mass.vector <- colSums(M)} else mass.vector <- M
  #Use in leapfrog
  N <- TotalSamples+burnin; #How many samples will be needed
  pb <- txtProgressBar(min = 0, max = N, style = 3);
  OldHamiltonian=rep(0,N);
  ProposeHamiltonian=rep(0,N);
  Samples <- matrix(0,N,d);
  #For allocation of memory
  p <- mvrnorm(N,rep(0,d),M); 
  #Generate momentum for multivariate normal, this is an N by d quantity. I generate all at once here.
  for (k in 1:N){
    if (N==1){
      oldp <-p.epsilon <- p
    }
    else {
      oldp <-p.epsilon <- p[k,] 
    }
    #Take one momentum vector
    
    #Varies epsilon to avoid periodic behaviour and also 'encourage' ergodicity
    epsilon <- runif(1,epsilonFix*0.8,epsilonFix*1.2);
    #Varies epsilon to avoid periodic behaviour and also 'encourage' ergodicity
    
    
    # Leapfrog iterated 'L' times with stepsize epsilon
    for (t in 1:L){
      p.epsilon.half <- p.epsilon - (epsilon/2) * as.vector(densitydiff(q.epsilon))
      q.epsilon <-  q.epsilon + ((epsilon / mass.vector) * p.epsilon.half)
      p.epsilon <- p.epsilon.half - (epsilon/2) * as.vector(densitydiff(q.epsilon))
    }
    
    
    OldUenergy= -log(densityA(oldq));
    #Energy of target densityA at old state
    OldKenergy= sum((oldp^2/mass.vector))/2;
    #Energy of momentum at Old state
    OldHamiltonian[k]= OldUenergy+OldKenergy;
    #Old Hamiltonian
    ProposeUenergy= -log(densityA(q.epsilon));
    ProposeKenergy=sum((p.epsilon^2)/mass.vector)/2;
    #Propose Energy 
    ProposeHamiltonian[k]=ProposeUenergy+ProposeKenergy;
    #Propose Hamiltonian
    if (runif(1)< exp(OldHamiltonian[k]-ProposeHamiltonian[k])) {
      oldp=p.epsilon;
      oldq=q.epsilon;
      #Metropolis Hastings, if accepted move to new state
    }
    else {q.epsilon=oldq;
    rejections= rejections+1;
    }
    Samples[k,] <- oldq;
    # if (k==1){
    #plot(Samples[k,1],Samples[k,2],pch=26,col=26,,xlab="x",ylab="y",xlim=c(-2,2), ylim=c(-2,2))}
    #else {
    #lines(Samples[k,1],Samples[k,2],col=26,,xlab="x",ylab="y")
    #points(Samples[k,1],Samples[k,2],pch=4,col=26,,xlab="x",ylab="y")
    # }
    # Sys.sleep(0.03)
    #In order to plot points sequentially
    #scatterplot3d(Samples[k,1],Samples[k,2],Samples[k,3],xlim=c(-4,4),ylim=c(-4,4),zlim=c(-4,4))
    #par(new=TRUE)
    #Samples
    #forget Intialisation
    setTxtProgressBar(pb, k)
  }
  #  cat("Total Samples after Burn in:",TotalSamples,"\nStep Size used was (+-10%)",epsilonFix,"\nLeapFrog iterations were",L,"\nMean Of Samples=",apply(Samples,2,mean),"\nCovariance=",cov(Samples), "\nRejectionRate=",rejections/N)
  #hist(Samples[((burnin+1):N),],freq=F)
  #points(Samples[,1],Samples[,2])
  #plot(Samples[,1],Samples[,2],col=2)
  #lines(seq(-4,4,0.01),dnorm(seq(-4,4,0.01)))
  #acf(Samples)
  plot(1:N,OldHamiltonian-ProposeHamiltonian)
  # plot(Samples[,1])
  #plot(Samples[,1],Samples[,2],type='l',col=2,xlab="X",ylab="Y")
  #points(Samples[,1],Samples[,2],pch=4)
  after=Sys.time()-now;
  print(after)
  close(pb)
  (Samples[((burnin+1):N),])
}



HMC <- function(TotalSamples, densityA, M, q, epsilon, L, densitydiff, burnin) {
  # Arguments
  # TotalSamples - the number of simulations you'd like output
  # densityA
  # M - the mass matrix
  # q - starting values
  # epsilon
  # L
  # densitydiff - the derivative of the distribution
  # burnin - the burn-in
  epsilonFix <- epsilon
  
  #now <- Sys.time()
  #For display at the end.
  #quartz()
  #Bring up window for plot
  #plot.new()
  #  plot.window(xlim=c(-3,3), ylim=c(-3,3))
  #  axis(side=1, pos=-3)
  #  axis(side=2,pos=-3)
  #  title(main="Bivariate Normal Simulation")
  rejections <- 0
  #Counts Rejection
  oldq <-q.epsilon <- q # Initialisation of q and q.epsilon
  d <- length(q) # Dimensions
  if(length(M)>1){mass.vector <- colSums(M)} else mass.vector <- M
  #Use in leapfrog
  N <- TotalSamples+burnin; #How many samples will be needed
  pb <- txtProgressBar(min = 0, max = N, style = 3);
  OldHamiltonian=rep(0,N);
  ProposeHamiltonian=rep(0,N);
  Samples <- matrix(0,N,d);
  #For allocation of memory
  p <- mvrnorm(N,rep(0,d),M); 
  #Generate momentum for multivariate normal, this is an N by d quantity. I generate all at once here.
  for (k in 1:N){
    if (N==1){
      oldp <-p.epsilon <- p
    }
    else {
      oldp <-p.epsilon <- p[k,] 
    }
    #Take one momentum vector.
    epsilon= runif(1,epsilonFix*0.8,epsilonFix*1.2);
    #Varies epsilon to avoid periodic behaviour and also 'encourage' ergodicity
    for (t in 1:L){
      #LeapFrog L times with step size epsilon
      p.epsilon.half= p.epsilon - (epsilon/2)*as.vector(densitydiff(q.epsilon));
      q.epsilon =  q.epsilon + ((epsilon/mass.vector)*p.epsilon.half);
      #Assume independent in p vector hence use column sums
      p.epsilon= p.epsilon.half - (epsilon/2)*as.vector(densitydiff(q.epsilon));
    }
    OldUenergy= -log(densityA(oldq));
    #Energy of target densityA at old state
    OldKenergy= sum((oldp^2/mass.vector))/2;
    #Energy of momentum at Old state
    OldHamiltonian[k]= OldUenergy+OldKenergy;
    #Old Hamiltonian
    ProposeUenergy= -log(densityA(q.epsilon));
    ProposeKenergy=sum((p.epsilon^2)/mass.vector)/2;
    #Propose Energy 
    ProposeHamiltonian[k]=ProposeUenergy+ProposeKenergy;
    #Propose Hamiltonian
    if (runif(1)< exp(OldHamiltonian[k]-ProposeHamiltonian[k])) {
      oldp=p.epsilon;
      oldq=q.epsilon;
      #Metropolis Hastings, if accepted move to new state
    }
    else {q.epsilon=oldq;
    rejections= rejections+1;
    }
    Samples[k,] <- oldq;
    # if (k==1){
    #plot(Samples[k,1],Samples[k,2],pch=26,col=26,,xlab="x",ylab="y",xlim=c(-2,2), ylim=c(-2,2))}
    #else {
    #lines(Samples[k,1],Samples[k,2],col=26,,xlab="x",ylab="y")
    #points(Samples[k,1],Samples[k,2],pch=4,col=26,,xlab="x",ylab="y")
    # }
    # Sys.sleep(0.03)
    #In order to plot points sequentially
    #scatterplot3d(Samples[k,1],Samples[k,2],Samples[k,3],xlim=c(-4,4),ylim=c(-4,4),zlim=c(-4,4))
    #par(new=TRUE)
    #Samples
    #forget Intialisation
    setTxtProgressBar(pb, k)
  }
  #  cat("Total Samples after Burn in:",TotalSamples,"\nStep Size used was (+-10%)",epsilonFix,"\nLeapFrog iterations were",L,"\nMean Of Samples=",apply(Samples,2,mean),"\nCovariance=",cov(Samples), "\nRejectionRate=",rejections/N)
  #hist(Samples[((burnin+1):N),],freq=F)
  #points(Samples[,1],Samples[,2])
  #plot(Samples[,1],Samples[,2],col=2)
  #lines(seq(-4,4,0.01),dnorm(seq(-4,4,0.01)))
  #acf(Samples)
  plot(1:N,OldHamiltonian-ProposeHamiltonian)
  # plot(Samples[,1])
  #plot(Samples[,1],Samples[,2],type='l',col=2,xlab="X",ylab="Y")
  #points(Samples[,1],Samples[,2],pch=4)
  after=Sys.time()-now;
  print(after)
  close(pb)
  (Samples[((burnin+1):N),])
}
JJ <- HMCVer1(TotalSamples = 1,densityA = dnorm,M=1,q = 0,epsilon = 0.05,L = 20,densitydiff = function(x) x,burnin = 100)

#Bivariate Normal
JJ <- HMCVer1(TotalSamples = 25,densityA = function(l) dmvnorm(l,c(0,0),matrix(c(1,0.95,0.95,1),2,2)),q = c(-2,-2),M=diag(2),epsilon = 0.18,L = 20,densitydiff = function(x) {solve(matrix(c(1,0.95,0.95,1),2,2))%*%as.matrix(x)}, burnin = 0)

JJ <- HMCVer1(TotalSamples = 1000,densityA = function(l) dmvnorm(l,c(0,0,0),matrix(c(1,0.9,0.8,0.9,1,0.7,0.8,0.7,1),3,3)),M=diag(3),q = c(10,10,-10),epsilon = 0.05,L = 60,densitydiff = function(x) {solve(matrix(c(1,0.9,0.8,0.9,1,0.7,0.8,0.7,1),3,3))%*%(as.matrix(x)-as.matrix(c(20,4,-10)))}, burnin = 0)
#50 dimensions
JJ <- HMCVer1(TotalSamples = 5000,densityA = function(l) dmvnorm(l,rep(0,150),diag(seq(from=0.02,to=1,length=150)^2)),q = rep(0,150),M=diag(150),epsilon = 0.014,L = 100,densitydiff = function(x) {solve(diag(seq(from=0.02,to=1,length=150)^2))%*%as.matrix(x)}, burnin = 0)
#Make matrix 50 by 50
MatrixMake=seq(from=0.02,to=1,length=150)^2
plot(MatrixMake[1],var(JJ[,1]),xlim=c(0,1),ylim=c(0,1),xlab="Real Variance of coordinate",ylab="Sample Variance",main="Hamiltonian Monte Carlo")
for (t in 2:150) {
  points(MatrixMake[t],var(JJ[,t]))
}
lines(MatrixMake,MatrixMake)
plot(MatrixMake,apply(JJ,2,mean),xlab="Real Variance",ylab="Mean",main="Hamiltonian Monte Carlo")