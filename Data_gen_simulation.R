rm(list = ls())
# Generate the data
library(extraDistr)
library(MASS)
L1.local.true = 6 # true number of local groups in population 1
L2.local.true = 7 # true number of local groups in population 2
L3.local.true = 5 # true number of local groups in population 3

L.global.true = 8 # true number of global groups in population
alpha0 = 25
m0 = 25

alpha.true = rgamma(n = 1, shape = alpha0, rate = 1) # True alpha for generating data
m.true = rgamma(n = 1, shape = m0, rate = 1) # True gamma for generating data

# True weights
beta.true = as.numeric(rdirichlet(n = 1, alpha = rep(m.true/L.global.true, L.global.true))) # True beta
pi1.true = as.numeric(rdirichlet(n = 1, alpha = rep(alpha.true/L1.local.true, L1.local.true))) # True pi1
pi2.true = as.numeric(rdirichlet(n = 1, alpha = rep(alpha.true/L2.local.true, L2.local.true))) # True pi2
pi3.true = as.numeric(rdirichlet(n = 1, alpha = rep(alpha.true/L3.local.true, L3.local.true))) # True pi3


# Sample sizes
n1 = 100 # First population
n2 = 110 # Second population
n3 = 115 # Third population

# True local-level indicators
t1.true = sample(1:L1.local.true, size = n1, prob = pi1.true, replace = TRUE)
t2.true = sample(1:L2.local.true, size = n2, prob = pi2.true, replace = TRUE)
t3.true = sample(1:L3.local.true, size = n3, prob = pi3.true, replace = TRUE)
# True global-level indicators
k1.true = sample(1:L.global.true, size = L.global.true, prob = beta.true, replace = TRUE)
k2.true = sample(1:L.global.true, size = L.global.true, prob = beta.true, replace = TRUE)
k3.true = sample(1:L.global.true, size = L.global.true, prob = beta.true, replace = TRUE)
# Draw data from Normal populations
p.global = 2 # Dimension of global variable
p.local1 = 1 # Dimension of local variable in population 1
p.local2 = 2 # Dimension of local variable in population 2
p.local3 = 3 # Dimension of local variable in population 3
# True variances for global variables
SigmaG0 = matrix(rinvgamma(p.global * L.global.true, alpha = 2, beta = 1), nrow = p.global)
# True variances for local variable in population 1
Sigma1L0 = matrix(rinvgamma(p.local1 * L1.local.true, alpha = 2, beta = 1), nrow = p.local1)
# True variances for local variable in population 2
Sigma2L0 = matrix(rinvgamma(p.local2 * L2.local.true, alpha = 2, beta = 1), nrow = p.local2)
# True variances for local variable in population 3
Sigma3L0 = matrix(rinvgamma(p.local3 * L3.local.true, alpha = 2, beta = 1), nrow = p.local3)
# prior mean to simulate the true mean of the global variables from its prior distribution
phi0 = rep(0, p.global)
# Parameter controlling the spread of the true mean of the global variables
lambda.global = 0.1
lambda.global.inv = 1/lambda.global

# Parameter controlling the spread of the true mean of the local variables
lambda.local = 0.1
lambda.local.inv = 1/lambda.local

# prior mean to simulate the true mean of the local variables from its prior distribution
mu = list(c(0),  c(0, 0), c(0, 0, 0))
library(MASS)
# Simulate the true mean for the local variables
mean1.local.true = sapply(1:L1.local.true, function(j){rnorm(n = 1, mean = mu[[1]], sd = sqrt(lambda.local.inv * Sigma1L0[,j]))}) # True mean for local variable in population 1

mean2.local.true =  t(sapply(1:L2.local.true, function(j){mvrnorm(n = 1, mu = mu[[2]], Sigma = diag(lambda.local.inv * Sigma2L0[,j], nrow = p.local2))})) # True mean for local variable in population 2
mean2.local.true = lapply(seq_len(nrow(mean2.local.true)), function(i) mean2.local.true[i, ])

mean3.local.true = t(sapply(1:L3.local.true, function(j){mvrnorm(n = 1, mu = mu[[3]], Sigma = diag(lambda.local.inv * Sigma3L0[,j], nrow = p.local3))})) # True mean for local variable in population 3
mean3.local.true = lapply(seq_len(nrow(mean3.local.true)), function(i) mean3.local.true[i, ])

# Simulate the true mean for the global variables
mean.global.true = t(sapply(1:L.global.true, function(j){mvrnorm(n = 1, mu = phi0, Sigma = diag(lambda.global.inv * SigmaG0[,j], nrow = p.global))}))# True mean local
mean.global.true = lapply(seq_len(nrow(mean.global.true)), function(i) mean.global.true[i, ])

## With the true mean and true variances, simulate the date from the different populations
X1 = sapply(1:n1, function(j){
  c(rnorm(n = 1, mean = mean1.local.true[t1.true[j]], sd = sqrt(Sigma1L0[, t1.true[j]])),
    mvrnorm(n = 1, mu = mean.global.true[[k1.true[t1.true][j]]], Sigma = diag(SigmaG0[, k1.true[t1.true][j]], p.global)))
}) # Data from population 1
X2 = sapply(1:n2, function(j){
  c(mvrnorm(n = 1, mu = mean2.local.true[[t2.true[j]]], Sigma = diag(Sigma2L0[ ,t2.true[j]], p.local2)),
    mvrnorm(n = 1, mu = mean.global.true[[k2.true[t2.true][j]]], Sigma = diag(SigmaG0[, k2.true[t2.true][j]], p.global)))
})    # Data from population 2 
X3 = sapply(1:n3, function(j){
  c(mvrnorm(n = 1, mu = mean3.local.true[[t3.true[j]]], Sigma = diag(Sigma3L0[, t3.true[j]], p.local3)),
    mvrnorm(n = 1, mu = mean.global.true[[k3.true[t3.true][j]]], Sigma = diag(SigmaG0[ , k3.true[t3.true][j]], p.global)))
})  # Data from population 3 


