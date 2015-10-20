m <- stan_model(model_code = 'data {
                int<lower=1> K;
                int<lower=1> d;
                row_vector[K] lambda;
                matrix[K,d] mu;
                matrix[d*K,d] Sigmas;
                }
                parameters {
                matrix[K,d] x;
                real U;
                }
                transformed parameters {
                row_vector[d] z;
                for(k in 1:K) {
                if(U<=lambda[K-k+1])
                z <- x[K-k+1];
                }
                }
                model {
                matrix[d,d] sigma;
                for(i in 1:K) {
                for(r in 1:d) {
                for(c in 1:d) {
                sigma[r,c] <- Sigmas[(d*(i-1)+r),c];
                }
                }
                x[i] ~ multi_normal(mu[i], sigma);
                }
                
                U ~ uniform(0,1);
                }')

letter <- A[[8]]
K <- (length(letter)-1)/3
lambda <- cumsum(letter[2:(1+K)]) # mixing proportions
d <- 2 # number of dimensions
mu <- letter[[2+K]] # mean vectors row by row
for(i in 2:K) { 
  mu <- rbind(mu,letter[[1+K+i]])
}
Sigmas <- letter[[2+2*K]]
for(i in 2:K) {
  Sigmas <- rbind(Sigmas,letter[[1+2*K+i]])
}


f <- sampling(m, iter = 5000,chains=2)
hist(f@sim$samples[[1]]$`z[3]`)
plot(f@sim$samples[[2]]$`z[1]`,f@sim$samples[[2]]$`z[2]`,asp = 1)


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
d <- 150
mu <- rep(0,150)
sigma <- diag(seq(from=0.02, to=1, length = 150)^2)

f <- sampling(m, iter = 5000,chains=2)
var(f@sim$samples[[1]]$`x[150]`)
hist(f@sim$samples[[1]]$`x[2]`)
plot(f@sim$samples[[2]]$`z[1]`,f@sim$samples[[2]]$`z[2]`,asp = 1)
