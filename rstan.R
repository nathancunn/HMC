dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXXFLAGS=-O3", file = M, sep = "\n", append = TRUE)

install.packages("rstan", dependencies = TRUE)


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


mu <- c(0,0)
mu <- as.vector(mu)
Sigma <- matrix(c(1,0.9,0.9,1),nrow=2)
m <- stan_model(model_code = 'data {row_vector[2] mu; matrix[2,2] Sigma;} parameters {row_vector[2] d;} model {d ~ multi_normal(mu,Sigma);}')

f <- sampling(m, iter = 1000)
pairs(f)
x <- f@sim$samples[[1]][[1]]
y <- f@sim$samples[[1]][[2]]

x2 <- f@sim$samples[[4]][[1]]
y2 <- f@sim$samples[[4]][[2]]
