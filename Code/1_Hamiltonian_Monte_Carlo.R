libs <- list("MASS","mvtnorm")
lapply(libs,library,character.only = T)
HMC <- function(total.samples, q.density, M, q, epsilon, L, diff.density, burnin, output = 1) {
  # Arguments
  # total.samples - the number of simulations you'd like output
  # q.density - the density of interest to sample from
  # M - the mass matrix
  # q - starting values
  # epsilon
  # L - the number of leapfrog iterations
  # diff.density - the derivative of the density to sample from
  # burnin - the burn-in
  # Output - whether you want text output
  init.epsilon <- epsilon
  rejections <- 0
  init.q <-q.epsilon <- q # Initialisation of q and q.epsilon
  d <- length(q) # Dimensions
  # The elements of p are assumed independent, hence use column sums
  if(length(M)>1){mass.vector <- colSums(M)} else mass.vector <- M

  N <- total.samples+burnin
  # Set a progress bar
  if(output==1) {pb <- txtProgressBar(min = 0, max = N, style = 3)}
  old.h=rep(0,N)
  proposal.h=rep(0,N)
  samples <- matrix(0,N,d)
  #For allocation of memory
  p <- mvrnorm(N,rep(0,d),M);
  #Generate momentum for multivariate normal, this is an N by d quantity. I generate all at once here.
  for (k in 1:N){
    if (N==1){
      init.p <-p.epsilon <- p
    }
    else {
      init.p <-p.epsilon <- p[k,] 
    }
    #Take one momentum vector
    
    #Varies epsilon to avoid periodic behaviour and also 'encourage' ergodicity
    epsilon <- runif(1,init.epsilon*0.8,init.epsilon*1.2);
    #Varies epsilon to avoid periodic behaviour and also 'encourage' ergodicity
    
    
    # Leapfrog iterated 'L' times with stepsize epsilon
    for (t in 1:L){
      p.epsilon.half <- p.epsilon - (epsilon/2) * as.vector(diff.density(q.epsilon))
      q.epsilon <-  q.epsilon + ((epsilon / mass.vector) * p.epsilon.half)
      p.epsilon <- p.epsilon.half - (epsilon/2) * as.vector(diff.density(q.epsilon))
    }
    
    # Potential energy at old state
    old.u <- -log(q.density(init.q))
    # Kinetic energy at old state
    old.k <- sum((init.p^2/mass.vector))/2
    # Old Hamiltonian
    old.h[k] <- old.u+old.k
    # The energies at the proposed states
    proposal.u <- -log(q.density(q.epsilon))
    proposal.k <- sum((p.epsilon^2)/mass.vector)/2
    # The Hamiltonian at the proposed state
    proposal.h[k] <- proposal.u+proposal.k
    # The Metropolis step for accepting the proposed new state
    if (runif(1) < exp(old.h[k]-proposal.h[k])) {
      # Accept proposed state
      init.p=p.epsilon
      init.q=q.epsilon
    }
    else {
    # Stay at current state and increase rejections tally
    q.epsilon=init.q
    rejections= rejections+1
    }
    # Update the samples and advance the progress bar
    samples[k,] <- init.q;
  }
  if(output==1){
  setTxtProgressBar(pb, k)
  cat("Total samples after Burn in:", total.samples, 
      "\nStep Size used was (+-20%)", init.epsilon, 
      "\nLeapFrog iterations were", L, 
      "\nMean Of samples=", apply(samples,2,mean), 
      "\n(Co)variance=", cov(samples), 
      "\nRejectionRate=",rejections/N)
  #after=Sys.time()-now;
  #print(after)
  close(pb)
  }
  samples[((burnin+1):N),]
}


# Some examples
out.univariate <- HMC(total.samples = 10000, 
          q.density = function(x) dnorm(x,0,1), 
          M=1, 
          q = 0, 
          epsilon = 0.05, 
          L = 20, 
          diff.density = function(x) x, 
          burnin = 100)

out.bivariate <- HMC(total.samples = 10000, 
                     q.density = function(l) dmvnorm(l, c(0, 0),matrix(c(1, 0.95, 0.95, 1), 2, 2)), 
                     q = c(-2,-2), 
                     M=diag(2), 
                     epsilon = 0.18, 
                     L = 20, 
                     diff.density = function(x) {solve(matrix(c(1, 0.95, 0.95, 1), 2, 2))%*%as.matrix(x)}, 
                     burnin = 0)


out.multidimension <- 
  HMC(total.samples = 500, 
      q.density = function(l) dmvnorm(l,rep(0,150),diag(seq(from=0.02,to=1,length=150)^2)), 
      q = rep(0,150), 
      M=diag(150), 
      epsilon = 0.014, 
      L = 100, 
      diff.density = function(x) {solve(diag(seq(from=0.02,to=1,length=150)^2))%*%as.matrix(x)}, 
      burnin = 0)




#Bivariate Normal
JJ <- HMC(total.samples = 10000,q.density = function(l) dmvnorm(l,c(0,0),matrix(c(1,0.95,0.95,1),2,2)),q = c(-2,-2),M=diag(2),epsilon = 0.18,L = 20,diff.density = function(x) {solve(matrix(c(1,0.95,0.95,1),2,2))%*%as.matrix(x)}, burnin = 0)

JJ <- HMCVer1(total.samples = 1000,q.density = function(l) dmvnorm(l,c(0,0,0),matrix(c(1,0.9,0.8,0.9,1,0.7,0.8,0.7,1),3,3)),M=diag(3),q = c(10,10,-10),epsilon = 0.05,L = 60,diff.density = function(x) {solve(matrix(c(1,0.9,0.8,0.9,1,0.7,0.8,0.7,1),3,3))%*%(as.matrix(x)-as.matrix(c(20,4,-10)))}, burnin = 0)
#50 dimensions
JJ <- HMC(total.samples = 500,q.density = function(l) dmvnorm(l,rep(0,150),diag(seq(from=0.02,to=1,length=150)^2)),q = rep(0,150),M=diag(150),epsilon = 0.014,L = 100,diff.density = function(x) {solve(diag(seq(from=0.02,to=1,length=150)^2))%*%as.matrix(x)}, burnin = 0)
#Make matrix 50 by 50
MatrixMake=seq(from=0.02,to=1,length=150)^2
plot(MatrixMake[1],var(JJ[,1]),xlim=c(0,1),ylim=c(0,1),xlab="Real Variance of coordinate",ylab="Sample Variance",main="Hamiltonian Monte Carlo")
for (t in 2:150) {
  points(MatrixMake[t],var(JJ[,t]))
}
lines(MatrixMake,MatrixMake)
plot(MatrixMake,apply(JJ,2,mean),xlab="Real Variance",ylab="Mean",main="Hamiltonian Monte Carlo")
