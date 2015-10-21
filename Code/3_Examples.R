# Some examples
# Sampling from a univariate standard Gaussian
out.univariate <- HMC(total.samples = 1000, 
                      q.density = function(x) dnorm(x,0,1), 
                      M=1, 
                      q = 0, 
                      epsilon = 0.05, 
                      L = 20, 
                      diff.density = function(x) x, 
                      burnin = 100)


# Sampling from a bivariate Gaussian
bivariate.density <- function(l) {
  dmvnorm(l, c(0, 0),matrix(c(1, 0.95, 0.95, 1), 2, 2))
}
bivariate.diff <- function(x) {
  solve(matrix(c(1, 0.95, 0.95, 1), 2, 2))%*%as.matrix(x)
}
out.bivariate <- HMC(total.samples = 10000, 
                     q.density = bivariate.density, 
                     q = c(-2,-2), 
                     M = diag(2), 
                     epsilon = 0.18, 
                     L = 20, 
                     diff.density = bivariate.diff, 
                     burnin = 0)

# Sampling from a 150-dimension Gaussian
multi.density <- function(l) {
  dmvnorm(l,rep(0,150),diag(seq(from=0.02,to=1,length=150)^2))
}
multi.diff <- function(x) {
  solve(diag(seq(from=0.02,to=1,length=150)^2))%*%as.matrix(x)
}
out.multidimension <- HMC(total.samples = 500, 
                          q.density = multi.density, 
                          q = rep(0,150), 
                          M=diag(150), 
                          epsilon = 0.014, 
                          L = 100, 
                          diff.density = multi.diff, 
                          burnin = 0)