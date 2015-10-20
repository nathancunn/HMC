library(rstan)
m <- stan_model(model_code = 'data {
                int<lower=1> K;
                int<lower=1> d;
                int<lower=1> ncomponents;
                row_vector[ncomponents] lambda;
                matrix[ncomponents,d] mu;
                matrix[d*ncomponents,d] Sigmas;
                }
                parameters {
                matrix[K,d] x;
                real U;
                }
                transformed parameters {
                row_vector[d] z;
                for(k in 1:ncomponents) {
                if(U<=lambda[ncomponents-k+1])
                z <- x[ncomponents-k+1];
                }
                }
                model {
                matrix[d,d] sigma;
                for(i in 1:ncomponents) {
                for(r in 1:d) {
                for(c in 1:d) {
                sigma[r,c] <- Sigmas[(d*i-2+r),c];
                }
                }
                x[i] ~ multi_normal(mu[i], sigma);
                }
                
                U ~ uniform(0,1);
                }')

# Creating words
# Need wordform function
letters <- wordform("ab")
init <- 1
d <- 2 # number of dimensions
for(i in 1:nrow(letters)){
letter <- A[[letters[i,1]]]
n <- length(letter)
N <- nrow(letters)
K <- (n-1)/3
if(i==1) {mu <- letter[[2+K]]+c((i-1)*20,0);  # mean vectors row by row
          Sigmas <- letter[[2+2*K]];
          lambda <- unlist(letter[2:(1+K)])/N;
          init <- 2
} else {
  lambda <- append(lambda, unlist(letter[2:(1+K)])/N) # mixing proportions
  init <- 1
}
for(j in init:K) { 
  mu <- rbind(mu,letter[[1+K+j]]+c((i-1)*20,0))
  Sigmas <- rbind(Sigmas,letter[[1+2*K+j]])
  
}
}
lambda <- cumsum(lambda)
ncomponents <- length(lambda)

f <- sampling(m, iter = 1000,chains=2)

plot(0,type="n",xlim=c(0,40),ylim=c(0,40),asp=1)
for(i in 1:K) {
  points(f@sim$samples[[2]][i][[1]],f@sim$samples[[2]][i+K][[1]])
  
}


# 150 dimensions
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
m <- stan_model(model_code = 'data {
                int<lower=0> d;
                row_vector[d] mu;
                matrix[d,d] sigma;
                }
                parameters {
                row_vector[d] x;
                }
                model {
                x ~ multi_normal(mu, sigma);
                }')

mu <- rep(0,150)
#mu = c(0,0)
sigma <- diag(seq(from=0.02, to=1, length = 150)^2)
#sigma = matrix(c(1,0.95,0.95,1),2,2))
d <- length(mu)
f <- sampling(m, iter = 5000,chains=1)
var(f@sim$samples[[1]]$`x[150]`)
hist(f@sim$samples[[1]]$`x[2]`)
plot(f@sim$samples[[2]]$`z[1]`,f@sim$samples[[2]]$`z[2]`,asp = 1)

MatrixMake=seq(from=0.02,to=1,length=150)^2
plot(MatrixMake[1],var(f@sim$samples[[1]]$`x[1]`),xlim=c(0,1),ylim=c(0,1),xlab="Real Variance of coordinate",ylab="Sample Variance",main="NUTS Monte Carlo",pch=4)
Information=f@sim$samples[[1]]
for (t in 2:150) {
  points(MatrixMake[t],var(Information[[t]]),pch=4)
}
lines(MatrixMake,MatrixMake)
plot(MatrixMake[1],mean(f@sim$samples[[1]]$`x[1]`),xlim=c(0,1),ylim=c(-0.3,0.3),xlab="Real Variance of coordinate",ylab="Sample Mean",main="NUTS Monte Carlo",pch=4)
for (t in 2:150) {
points(MatrixMake[t],mean(Information[[t]]),pch=4)}
lines(MatrixMake,rep(0,150))
attributes(Information)
attributes(Information)$sampler_params$n_leapfrog_
ATTRIBUTESNUTSEPSILON=attributes(Information)$sampler_params$stepsize__
NUTSLEAPFROG=attributes(Information)$sampler_params$n_leapfrog_