# Code for using letter plots
setwd("HMC/")

# Load the Gaussians specifying the letters
letter.models <- source("Output/output")[[1]]



MetropolisHastingC = function(TamplesN, densityInt, sigma, initial,thin) {
  # now=Sys.time()
  #Multivariate Normal with Samples, integration, (sigma in this case for 2-dim),initial conditions, thinning (1 if no thinning)
  SamplesN=thin*TamplesN;
  #Works out how much thinning
  OldPos=initial;
  OutputSamples= matrix(0,SamplesN,length(initial));
  # pb <- txtProgressBar(min = 0, max = SamplesN, style = 3);
  #Allocation
  Dimensions=length(initial);
  for (i in 1:SamplesN)
  {
    Proposal<- as.vector(OldPos)+as.vector(runif(Dimensions,-sigma,sigma))
    #Proposal Value using uniform proposal (2-dimesional))
    #rmvnorm(1,as.vector(OldPos),sigma);
    ProposedValue=densityInt(as.vector(Proposal));
    #Proposed Value Density
    OldValue=densityInt(as.vector(OldPos));
    #Old Value Density
    if (runif(1)<= (ProposedValue/OldValue))
      #MH rejection
    {
      OldPos=Proposal;
    }
    OutputSamples[i,]=OldPos
    #setTxtProgressBar(pb, i)
  }
  if (TamplesN==1)
  {
    FinalOutputSamples=OutputSamples[thin,]
  }
  else {
    #Size=dim(unique(OutputSamples))
    FinalOutputSamples=OutputSamples[seq(1,dim(OutputSamples)[1],thin),] }
  #Thinning
  #plot(FinalOutputSamples[,1],FinalOutputSamples[,2],type='l',col=2,xlab="X",ylab="Y",xlim=c(-2,2),ylim=c(-2,2))
  #points(FinalOutputSamples[,1],FinalOutputSamples[,2],pch=4)
  #cat("rejectionRate=",1-(Size[1]/SamplesN))
  #after=Sys.time()-now;
  #print(after)
  return(FinalOutputSamples)
  #close(pb)
}

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


# A function to map letters from a word according to the key
WordForm <- function(word) {
  # Split the word into individual characters
  chars <- strsplit(word, "")[[1]]
  out <- matrix(0, nrow=length(chars))
  # Map each letter to its key
  for(i in 1:length(chars)) {
    if(any(all.letters$key[all.letters$letter==chars[i]])) {
    out[i,1] <- all.letters$key[all.letters$letter==chars[i]]
    } else out[i,1] <- 55
  }
  out
}



# The function to generate a dataset of random points from a distribution
# whose probability density spells a given word
WordPrint <- function(word,samples) {
  # Arguments
  # word - the word you want plotted
  # samples - the number of data points you want simulated
  letters.key <- WordForm(word)
  # An array of the initial values for the Hamiltonian Monte Carlo
  last.h = array(0,c(nchar(word),7,2))
  # A matrix containing the simulated sample
  out <- matrix(0,nrow=samples,ncol=5)
  # The number of letters
  n <- length(letters.key)
  # Give each letter equal weight
  eq_probs <- seq(1/n,1,length=n)
  progress.bar <- txtProgressBar(min = 0, max = samples, style = 3)
  for (i in 1:samples) {
  # Randomly go to a single letter
  letter <- findInterval(runif(1),eq_probs)+1
  j <- letters.key[letter]
  # Displace letters so they can be seen side-by-side
  displacement <- (0:(n-1))*18-20
  # Extract the number of components the letter is made of 
  atts <- (length(letter.models[[j]])-1)/3
  # Find the mixing proportions
  lambdas <- cumsum(letter.models[[j]][2:(2+atts-1)])
  # Go to one component based on its lambda probability
  group <-  findInterval(runif(1),lambdas)+1
  # Extract the (displaced) mean and variances for this component
  mean.group <- letter.models[[j]][1+atts+group][[1]]+c(displacement[letter],0)
  var.group <- letter.models[[j]][1+2*atts+group][[1]]
  # If we have no initial value set them to the prior mean
  if (last.h[letter,group,1]==0 & last.h[letter,group,2]==0) {
    last.h[letter,group,]=as.vector(mean.group)
  }
  # Simulate a single datum from the selected Gaussian mixture
  out[i,1:2] <- HMC(total.samples = 1,
             q.density = function(l) dmvnorm(l,as.vector(mean.group),var.group), 
             q = as.vector(last.h[letter,group,]), M=diag(2),
             epsilon = 0.05, 
             L = 60,
             diff.density = function(x) {
               solve(var.group)%*%(as.matrix(x)-as.matrix(mean.group))}, 
             burnin = 0, output = 0)
  # Store the value to be used as the initial value when this mixture is used again
  last.h[letter,group,] <- out[i,1:2]
  # Store a unique colour for each mixture to see where the mixtures are
  out[i,3:5] <- c(letter/n,group/6,0.2)
  setTxtProgressBar(progress.bar, i)
  }
  plot(out[,1],out[,2],col=rgb(out[,3],out[,4],out[,5]))
  return(out)
  close(progress.bar)
}

WordPrintMH <- function(word, samples, thin = 1) {
  # Arguments
  # word - the word you want plotted
  # samples - the number of data points you want simulated
  letters.key <- WordForm(word)
  # An array of the initial values for the Hamiltonian Monte Carlo
  last.h = array(0,c(nchar(word),7,2))
  # A matrix containing the simulated sample
  out <- matrix(0,nrow=samples,ncol=5)
  # The number of letters
  n <- length(letters.key)
  # Give each letter equal weight
  eq_probs <- seq(1/n,1,length=n)
  progress.bar <- txtProgressBar(min = 0, max = samples, style = 3)
  for (i in 1:samples) {
    # Randomly go to a single letter
    letter <- findInterval(runif(1),eq_probs)+1
    j <- letters.key[letter]
    # Displace letters so they can be seen side-by-side
    displacement <- (0:(n-1))*18-20
    # Extract the number of components the letter is made of 
    atts <- (length(letter.models[[j]])-1)/3
    # Find the mixing proportions
    lambdas <- cumsum(letter.models[[j]][2:(2+atts-1)])
    # Go to one component based on its lambda probability
    group <-  findInterval(runif(1),lambdas)+1
    # Extract the (displaced) mean and variances for this component
    mean.group <- letter.models[[j]][1+atts+group][[1]]+c(displacement[letter],0)
    var.group <- letter.models[[j]][1+2*atts+group][[1]]
    # If we have no initial value set them to the prior mean
    if (last.h[letter,group,1]==0 & last.h[letter,group,2]==0) {
      last.h[letter,group,]=as.vector(mean.group)
    }
    # Simulate a single datum from the selected Gaussian mixture
    out[i,1:2] <- MetropolisHastingC(TamplesN = 1,
                                     densityInt = function(l) dmvnorm(l,as.vector(mean.group),var.group), 
                      initial = as.vector(last.h[letter,group,]), 
                      sigma = sqrt(sqrt(max(as.vector(var.group)))), 
                      thin = thin)
    # Store the value to be used as the initial value when this mixture is used again
    last.h[letter,group,] <- out[i,1:2]
    # Store a unique colour for each mixture to see where the mixtures are
    out[i,3:5] <- c(letter/n,group/6,0.2)
    setTxtProgressBar(progress.bar, i)
  }
  plot(out[,1],out[,2],col=rgb(out[,3],out[,4],out[,5]))
  return(out)
  close(progress.bar)
}



samples <- 2500
word <- "OxWaSP"
HMC.word <- WordPrint(word,samples)
MH.word <- WordPrintMH(word,samples, thin = 1)
x11()
plot(0,type="n",xlim = range(HMC.word[,1]),ylim = range(HMC.word[,2]*2), 
     xlab = "",ylab = "", xaxt = "n", yaxt = "n", bty = "n")
text(x = mean(range(HMC.word[,1])), y = 55, labels = "Hamiltonian Monte Carlo")
text(x = mean(range(HMC.word[,1])), y = 5, labels = "Metropolis-Hastings")
abline(a = max(HMC.word[,2]), b = 0, col = rgb (0.7, 0.7, 0.7, 0.7))
for (i in 1:nrow(HMC.word)) {
  points(HMC.word[i,1],HMC.word[i,2]+25,col=rgb(HMC.word[i,3],HMC.word[i,4],HMC.word[i,5]))
  points(MH.word[i,1],MH.word[i,2],col=rgb(MH.word[i,3],MH.word[i,4],MH.word[i,5]))
  Sys.sleep(0.005)
}

