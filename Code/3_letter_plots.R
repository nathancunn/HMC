# Code for using letter plots
setwd("HMC/")

# Load the Gaussians specifying the letters
letter.models <- source("Output/output")[[1]]

# Creating a key matching the letters to their location in letter.models
all.letters <- data.frame(key=1:54)
for (i in 1:54){
all.letters$letter[i] <- letter.models[[i]]$letter
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




WordPrint <- function(word,samples) {
  letters.key <- WordForm(word)
  last.h = array(0,c(nchar(word),7,2))
  out <- matrix(0,nrow=samples,ncol=2)
  n <- length(letters.key)
  # Give each letter equal weight
  eq_probs <- seq(1/n,1,length=n)
  # Randomly go to a single letter
  for (i in 1:samples) {
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
  out[i,] <- HMC(total.samples = 1,
             q.density = function(l) dmvnorm(l,as.vector(mean.group),var.group), 
             q = as.vector(last.h[letter,group,]), M=diag(2),
             epsilon = 0.05, 
             L = 60,
             diff.density = function(x) {
               solve(var.group)%*%(as.matrix(x)-as.matrix(mean.group))}, 
             burnin = 0, output = 0)
  last.h[letter,group,] <- out[i,]
  }
  out
}

plot(WordPrint("OxWaSP",1000))



