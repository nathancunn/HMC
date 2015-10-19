dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXXFLAGS=-O3", file = M, sep = "\n", append = TRUE)

install.packages("rstan", dependencies = TRUE)


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())





K <- 2 # number of mixing groups
lambda <- c(0.3,1) # mixing proportions
d <- 2 # number of dimensions
mu <- matrix(c(0,10,20,30), nrow = 2, byrow = T) # mean vectors row by row
Sigmas <- matrix(data = c(1,0.9,0.9,1,1,-0.9,-0.9,1), ncol = 2, byrow = T) # covariance matrices row by row
m <- stan_model(model_code = 'data {
                int<lower=1> K;
                int<lower=1> d;
                row_vector[K] lambda;
                matrix[K,d] mu;
                matrix[d*K,K] Sigmas;
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
f <- sampling(m, iter = 10000,chains=1)
hist(f@sim$samples[[1]]$`z[1]`)
plot(f@sim$samples[[1]]$`z[1]`,f@sim$samples[[1]]$`z[2]`)
