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
