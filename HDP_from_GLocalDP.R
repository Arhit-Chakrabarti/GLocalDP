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
# True variances for global variables
SigmaG0 = matrix(rinvgamma(p.global * L.global.true, alpha = 2, beta = 1), nrow = p.global)
# prior mean to simulate the true mean of the global variables from its prior distribution
phi0 = rep(0, p.global)
# Parameter controlling the spread of the true mean of the global variables
lambda.global = 0.1
lambda.global.inv = 1/lambda.global
library(MASS)

# Simulate the true mean for the global variables
mean.global.true = t(sapply(1:L.global.true, function(j){mvrnorm(n = 1, mu = phi0, Sigma = diag(lambda.global.inv * SigmaG0[,j], nrow = p.global))}))# True mean local
mean.global.true = lapply(seq_len(nrow(mean.global.true)), function(i) mean.global.true[i, ])

## With the true mean and true variances, simulate the date from the different populations
X1 = sapply(1:n1, function(j){
  c(mvrnorm(n = 1, mu = mean.global.true[[k1.true[t1.true][j]]], Sigma = diag(SigmaG0[, k1.true[t1.true][j]], p.global)))
}) # Data from population 1
X2 = sapply(1:n2, function(j){
  c(mvrnorm(n = 1, mu = mean.global.true[[k2.true[t2.true][j]]], Sigma = diag(SigmaG0[, k2.true[t2.true][j]], p.global)))
})    # Data from population 2 
X3 = sapply(1:n3, function(j){
  c(mvrnorm(n = 1, mu = mean.global.true[[k3.true[t3.true][j]]], Sigma = diag(SigmaG0[ , k3.true[t3.true][j]], p.global)))
})  # Data from population 3 


library(Rcpp)
sourceCpp("functions.cpp") # Load the Rcpp functions needed for the sampler
################################################################################
# Some auxillary functions
################################################################################
# Function to perform log-sum-exp trick while calculating probabilities to avoid numerical issues
log_sum <- function(x){
  exp(x - max(x))
}
# Indicator function
indic <- function(x, t){
  ifelse(x == t, 1, 0)
}
# Least-square method of clustering function
getDahl <- function(z.list, z.true){
  membershipMatrices <- lapply(z.list, function(x){
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
  })
  
  membershipAverage <- Reduce("+", membershipMatrices)/length(z.list)
  
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  
  DahlIndex <- which.min(SqError) # Which index post burn-in which gives minimun Dahl index
  
  library(fossil)
  library(aricode)
  if(is.null(z.true)){
    return(list(DahlIndex = DahlIndex))
  }else{
    return(list(rand = rand.index(z.list[[DahlIndex]], z.true), # Rand Index by Vanilla DP
                DahlIndex = DahlIndex,
                adj.rand = aricode::ARI(z.list[[DahlIndex]], z.true)))
  }
}
# Function to sample the concentration parameter alpha using MH
# @pars argument should be a named list with a0, b0, J, pi1, pi2, pi3
# @ycurrent corresponds to the previous value of alpha
sample_alpha <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  pi1 = as.numeric(pars$pi1); pi2 = as.numeric(pars$pi2); pi3 = as.numeric(pars$pi3)
  a0 = pars$a0; b0 = pars$b0
  J = pars$J
  L = length(pi1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = - b0 * ycand + (a0 - 1)*log(ycand) + sum(((ycand/L) - 1) * log(pi1)) + sum((ycand/L - 1) * log(pi2)) + sum((ycand/L - 1) * log(pi3)) + J * lgamma(ycand) - (J * L * lgamma(ycand/L)) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = - b0 * ycurrent + (a0 - 1)*log(ycurrent) + sum(((ycurrent/L) - 1) * log(pi1)) + sum((ycurrent/L - 1) * log(pi2)) + sum((ycurrent/L - 1) * log(pi3)) + J * lgamma(ycurrent) - (J * L * lgamma(ycurrent/L))  + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}
# Function to sample the concentration parameter gamma using MH
# @pars argument should be a named list with a1, b1 and beta
# @ycurrent corresponds to the previous value of m
sample_m <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  a1 = pars$a1; b1 = pars$b1
  beta = pars$beta
  L = length(beta) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a1, rate = b1) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = - b1 * ycand + (a1 - 1)*log(ycand) + sum(((ycand/L) - 1) * log(beta)) + lgamma(ycand) - (L * lgamma(ycand/L)) + dgamma(ycurrent, shape = a1, rate = b1, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = - b1 * ycurrent + (a1 - 1)*log(ycurrent) + sum(((ycurrent/L) - 1) * log(beta)) + lgamma(ycurrent) - (L * lgamma(ycurrent/L)) + dgamma(ycand, shape = a1, rate = b1, log = TRUE)
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}

###############################################################################
# Save global variables as a list across the populations
X.global = list(X1, X2, X3) 
p.global = nrow(X.global[[1]]) # Dimension of global variables
n1 = ncol(X.global[[1]]) # Sample size in population 1
n2 = ncol(X.global[[2]]) # Sample size in population 2
n3 = ncol(X.global[[3]]) # Sample size in population 3

###############################################################################
# Set the parameters of the GLocal DP
###############################################################################
L = 10 # Truncation level of GLocal DP 
J = 3 # Number of populations

# Specify the hyper-priors
alpha0 = 0.1; beta0 = 0.1; m0 = rep(0, p.global) # Prior hyper-parameters corresponding to global parameters
alpha.j0 = 0.1; beta.j0 = 0.1; m.j0 = list(0, c(0, 0), c(0, 0, 0)) # Prior hyper-parameters corresponding to local parameters

lambda.global = 0.01 # Global variable prior-precision
lambda.local = rep(0.01, J) # Local variable prior-precision acoss J populations
library(extraDistr)
# MCMC initial values of parameters
alpha.start = 2 # Starting value of concentration parameter alpha
m.start = 2 # Starting value of concentration parameter gamma

library(extraDistr)
beta.start = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting beta
pi.start = replicate(n = J, list())
pi.start[[1]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi1
pi.start[[2]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi2
pi.start[[3]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi3

# Initialization of local-level indicators for the different populations
t.start = replicate(n = J, list())
t.start[[1]] = sample(1:L, size = n1, prob = pi.start[[1]], replace = TRUE)
t.start[[2]] = sample(1:L, size = n2, prob = pi.start[[2]], replace = TRUE)
t.start[[3]] = sample(1:L, size = n3, prob = pi.start[[3]], replace = TRUE)

# Initialization of global-level indicators for the different populations
k.start = replicate(n = J, list())
k.start[[1]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)
k.start[[2]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)
k.start[[3]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)

num_iter = 50000 # Number of MCMC iterations

# Initialize a list to store the MCMC samples of t_ji
t.samples.HDP <- replicate(n = J, list(replicate(n = num_iter, list(0))))

t.samples.HDP[[1]][[1]] = t.start[[1]]  # Initialize the MCMC starting value for population 1
t.samples.HDP[[2]][[1]] = t.start[[2]]  # Initialize the MCMC starting value for population 2
t.samples.HDP[[3]][[1]] = t.start[[3]]  # Initialize the MCMC starting value for population 3

# Initialize a list to store the MCMC samples of k_jt
k.samples.HDP <- replicate(n = J, list(replicate(n = num_iter, list(0))))

k.samples.HDP[[1]][[1]] = k.start[[1]]  # Initialize the MCMC starting value for population 1
k.samples.HDP[[2]][[1]] = k.start[[2]]  # Initialize the MCMC starting value for population 2
k.samples.HDP[[3]][[1]] = k.start[[3]]  # Initialize the MCMC starting value for population 3

# Initialize a list to store the MCMC samples of pi_j
pi.samples.HDP <- replicate(n = J, list(replicate(n = num_iter, list())))
pi.samples.HDP[[1]][[1]] = pi.start[[1]]  # Initialize the MCMC starting value for population 1
pi.samples.HDP[[2]][[1]] = pi.start[[2]]  # Initialize the MCMC starting value for population 2
pi.samples.HDP[[3]][[1]] = pi.start[[3]]  # Initialize the MCMC starting value for population 3

# Initialize a list to store the MCMC samples of beta
beta.samples.HDP <- list(); beta.samples.HDP[[1]] <- beta.start
# Initialize List to store MCMC samples of global parameter (means and variances)
M.k.samples.HDP <- list(); Tau2.k.samples.HDP <- list() 
# Initialize List to store MCMC samples of alpha
alpha.samples.HDP <- list(); alpha.samples.HDP[[1]] = alpha.start
# Initialize List to store MCMC samples of gamma
m.samples.HDP <- list(); m.samples.HDP[[1]] = m.start

################################################################################
# SAMPLING 
################################################################################
time.start = Sys.time()
for(iter in 2:num_iter){
  # Printing the iterations
  if(iter == 2){
    cat(paste0("Iteration: ", iter-1, "\n"))
  }
  if(iter %% floor((10/100)*(num_iter + 1)) == 0) {
    cat(paste0("Iteration: ", iter, "\n"))
  }
  
  m.HDP = replicate(n = J, list(0))
  
  for(j in 1:J){
    for(l in 1:L){
      m.HDP[[j]][l]= sum(t.samples.HDP[[j]][[iter - 1]] == l)
    }
  }
  
  
  for(j in 1:J){
    pi.samples.HDP[[j]][[iter]] = as.numeric(rdirichlet(n = 1, alpha = m.HDP[[j]] + alpha.samples.HDP[[iter - 1]]/L))
  }
  
  for(j in 1:J){
    pi.samples.HDP[[j]][[iter]] = (pi.samples.HDP[[j]][[iter]] + 1e-8)/sum(pi.samples.HDP[[j]][[iter]] + 1e-8)
  }
  
  d = replicate(n = J, list(0))
  for(j in 1:J){
    for(k in 1:L){
      d[[j]][k] = sum(k.samples.HDP[[j]][[iter - 1]] == k)
      
    } 
  }
  
  d.all = Reduce(`+`, d)
  
  beta.samples.HDP[[iter]] = as.numeric(rdirichlet(n = 1, alpha = d.all + m.samples.HDP[[iter - 1]]/L))
  
  beta.samples.HDP[[iter]] = (beta.samples.HDP[[iter]] + 1e-8)/sum(beta.samples.HDP[[iter]] + 1e-8)
  
  n.k = replicate(n = J, list(0))
  for(j in 1:J){
    for(k in 1:L){
      n.k[[j]][k] = sum(k.samples.HDP[[j]][[iter - 1]][t.samples.HDP[[j]][[iter - 1]]] == k)
    }
  }
  
  n.k.all = Reduce(`+`, n.k)
  
  tau2.k.samples <- matrix(0, p.global, L)
  m.k.samples <- matrix(0, p.global, L)
  
  x.k.bar = list()
  
  for(k in 1:L){
    if(n.k.all[k] == 0){
      x.k.bar[[k]] = rep(0, nrow(X.global[[1]]))
    }else{
      x.k.bar[[k]] = (rowSums(X.global[[1]][, k.samples.HDP[[1]][[iter - 1]][t.samples.HDP[[1]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[2]][, k.samples.HDP[[2]][[iter - 1]][t.samples.HDP[[2]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[3]][, k.samples.HDP[[3]][[iter - 1]][t.samples.HDP[[3]][[iter - 1]]] == k, drop = FALSE]))/n.k.all[k]
    }
  }
  
  S.k = matrix(0, nrow = p.global, ncol = L)
  Z.k = matrix(0, nrow = p.global, ncol = L)
  alpha.hat = matrix(0, nrow = p.global, ncol = L); beta.hat = matrix(0, nrow = p.global, ncol = L); lambda.hat = matrix(0, nrow = p.global, ncol = L); m.k.hat = list()
  
  for(k in 1:L){
    # k = 1
    if(n.k.all[k] != 0){
      S.k.j = matrix(0, nrow = p.global, ncol = J)
      for(j in 1:J){
        S.k.j[,j] = rowSums((X.global[[j]][, k.samples.HDP[[j]][[iter - 1]][t.samples.HDP[[j]][[iter - 1]]] == k, drop = FALSE] - x.k.bar[[k]])^2)
        
      }
      S.k[,k] = rowSums(S.k.j)
      Z.k[,k] = ((lambda.global * n.k.all[k])/(lambda.global + n.k.all[k])) * rowSums(matrix((x.k.bar[[k]] - m0)^2, ncol = 1))
    }
    alpha.hat[, k] = alpha0 + 0.5 * n.k.all[k]
    beta.hat[, k] = beta0 + 0.5 * (S.k[ , k] + Z.k[ , k])
    lambda.hat[, k] = lambda.global + n.k.all[k]
    m.k.hat[[k]] = ((lambda.global * m0) + (n.k.all[k] * x.k.bar[[k]]))/lambda.hat[, k]
  }
  
  for(i in 1:p.global){
    for(k in 1:L){
      tau2.k.samples[i, k] = rinvgamma(n = 1, alpha = alpha.hat[i, k], beta = beta.hat[i, k])
      m.k.samples[i, k] = rnorm(n = 1, mean = m.k.hat[[k]][i], sd = sqrt(tau2.k.samples[i, k]/lambda.hat[i, k]))
    }
  }
  
  
  Tau2.k.samples.HDP[[iter]] = tau2.k.samples
  M.k.samples.HDP[[iter]] = m.k.samples
  
  # Calculate the probabilities of table indices and sample
  n = c(n1, n2, n3)
  
  # Define an array of global cluster variance samples 
  Tau2.k.samples.array = array(0, dim = c(p.global, p.global, L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:ncol(Tau2.k.samples.HDP[[iter]])){
    Tau2.k.samples.array[,,j] = diag(Tau2.k.samples.HDP[[iter]][,j], nrow = p.global)
  }
  
  t.prob = replicate(n = J, list())
  for(j in 1:J){
    t.prob[[j]] = matrix(0, nrow = n[j], ncol = L)
  }
  
  for(j in 1:J){
    t.prob[[j]] <- prob_exponent_no_local2(pi_samples = pi.samples.HDP[[j]][[iter]], phi_samples = matrix(m.k.samples, nrow = nrow(X.global[[j]]), L),
                                           X_global = t(X.global[[j]]), k = k.samples.HDP[[j]][[iter - 1]], Sigma_G = Tau2.k.samples.array)
  }
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){
    t.prob[[j]] = calc_probability_log_sum_exp_normalized(t.prob[[j]])
  }
  
  for(j in 1:J){
    t.samples.HDP[[j]][[iter]] = sample_my(t.prob[[j]])
  }
  
  
  # Calculate the probabilities of dish indices and sample
  prob = array(0, dim = c(L, L, J))
  
  for(j in 1:J){
    prob[, , j] = prob_exponent_mv_global(beta_samples = beta.samples.HDP[[iter]], 
                                          X_global = t(X.global[[j]]), 
                                          phi_samples = matrix(m.k.samples, nrow = nrow(X.global[[j]]), L), 
                                          Sigma_G = Tau2.k.samples.array,
                                          t_samples = t.samples.HDP[[j]][[iter]])}
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){prob[, ,j] = calc_probability_log_sum_exp_normalized(prob[, ,j])}
  
  for(j in 1:J){
    k.samples.HDP[[j]][[iter]] = sample_my(prob[, , j])}
  
  sampled.alpha = sample_alpha(ycurrent = alpha.samples.HDP[[iter - 1]], 
                               pars = list(a0 = 0.1, b0 = 0.1,
                                           pi1 = pi.samples.HDP[[1]][[iter]],
                                           pi2 = pi.samples.HDP[[2]][[iter]],
                                           pi3 = pi.samples.HDP[[3]][[iter]], J = 3))
  alpha.samples.HDP[[iter]] = sampled.alpha$out
  
  sampled.m = sample_m(ycurrent = m.samples.HDP[[iter - 1]], pars = list(a1 = 0.1, b1 = 0.1, beta = beta.samples.HDP[[iter]]))
  m.samples.HDP[[iter]] = sampled.m$out
  
}
time.end = Sys.time()
time.end - time.start
################################################################################
# CALCULATE LOG-LIKELIHOOD TO PLOT TRACEPLOT OF LOG-LIKELIHOOD
################################################################################
ll.HDP = 0

for(iter in 2:num_iter){
  # Printing the iterations
  if(iter == 2){
    cat(paste0("Iteration: ", iter-1, "\n"))
  }
  if(iter %% floor((10/100)*(num_iter + 1)) == 0) {
    cat(paste0("Iteration: ", iter, "\n"))
  }
  
  # Define an array of global cluster variance samples 
  Tau2.k.samples.array = array(0, dim = c(p.global, p.global, L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:ncol(Tau2.k.samples.HDP[[iter]])){
    Tau2.k.samples.array[,,j] = diag(Tau2.k.samples.HDP[[iter]][, j], nrow = p.global)
  }
  # To calculate the LL Psi-samples are replaced by the M.j.samples (Mean of Local NIG posterior)
  # Phi-samples are replaced by M.k.samples (Mean of Global NIG posterior)
  
  logsum = logll_no_local2(pi.samples.HDP[[1]][[iter]],
                           t.samples.HDP[[1]][[iter]], matrix(M.k.samples.HDP[[iter]], nrow = p.global, ncol = L), X.global[[1]], k.samples.HDP[[1]][[iter]], 
                           Sigma_G = Tau2.k.samples.array) + 
    logll_no_local2(pi.samples.HDP[[2]][[iter]],
                    t.samples.HDP[[2]][[iter]], matrix(M.k.samples.HDP[[iter]], nrow = p.global, ncol = L), X.global[[2]], k.samples.HDP[[2]][[iter]], 
                    Sigma_G = Tau2.k.samples.array) + 
    logll_no_local2(pi.samples.HDP[[3]][[iter]],
                    t.samples.HDP[[3]][[iter]], matrix(M.k.samples.HDP[[iter]], nrow = p.global, ncol = L), X.global[[3]], k.samples.HDP[[3]][[iter]], 
                    Sigma_G = Tau2.k.samples.array)
  
  m1.tilde = count_my(t.samples.HDP[[1]][[iter]], L)
  m2.tilde = count_my(t.samples.HDP[[2]][[iter]], L)
  m3.tilde = count_my(t.samples.HDP[[3]][[iter]], L)
  
  logsum = logsum + sum(m1.tilde * log(pi.samples.HDP[[1]][[iter]])) + sum(m2.tilde * log(pi.samples.HDP[[2]][[iter]])) + sum(m3.tilde * log(pi.samples.HDP[[3]][[iter]]))
  
  m.tilde = count_my(k.samples.HDP[[1]][[iter]], L) + count_my(k.samples.HDP[[2]][[iter]], L) + count_my(k.samples.HDP[[3]][[iter]], L) 
  
  logsum = logsum + sum(m.tilde * log(beta.samples.HDP[[iter]]))
  
  ll.HDP[iter] = logsum
}

################################################################################
# PLOT TRACEPLOT OF LOG-LIKELIHOOD AND ACF
################################################################################
burn = 25000 # Burn-in 
samples.thin = seq((burn + 1), (num_iter), by = 25)

log_like <- ll.HDP[samples.thin]

library(tidyverse)
# RUN TrueLL_MVNIG.R BEFORE THIS PLOT
ll_plot <- data.frame(x = 1:length(log_like), y = log_like) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(x = "Iteration post burn-in", y = "") +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )


library(forecast)
ACF_plot <-  ggAcf(x = log_like, lag.max = 40) + labs(title = "", y = "") +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )


gridExtra::grid.arrange(ll_plot, ACF_plot, ncol = 2)


## PLOT THE TRACEPLOT OF THE CONCENTRATION PARAMETERS ALPHA AND GAMMA
library(latex2exp)
alpha_plot <- data.frame(alpha = unlist(alpha.samples.HDP[samples.thin]),
                         Iteration = 1:length(samples.thin)) %>% ggplot(aes(x = Iteration, y = alpha)) + geom_line() +
  labs(
    x = "Iteration post burn-in", y = TeX("$\\alpha$")) +
  theme_classic() +
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),
    axis.title.y = element_text(size=16, face="bold", colour = "black"),
    axis.text.x = element_text(size=16, face="bold", colour = "black"),
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  )


gamma_plot <- data.frame(m = unlist(m.samples.HDP[samples.thin]),
                         Iteration = 1:length(samples.thin)) %>% ggplot(aes(x = Iteration, y = m)) + geom_line() +
  labs(
    x = "Iteration post burn-in", y = TeX("$\\gamma$")) +
  theme_classic() +
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),
    axis.title.y = element_text(size=16, face="bold", colour = "black"),
    axis.text.x = element_text(size=16, face="bold", colour = "black"),
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  )

gridExtra::grid.arrange(alpha_plot, gamma_plot, ncol = 2)
################################################################################
# CLUSTERING PERFORMANCE
################################################################################
x.limit.lower <- min(X1[1, ], X2[1, ], X3[1, ])
x.limit.upper <- max(X1[1, ], X2[1, ], X3[1, ])

y.limit.lower <- min(X1[2, ], X2[2, ], X3[2, ])
y.limit.upper <- max(X1[2, ], X2[2, ], X3[2, ])


myvalues = c("1" = "#F8766D",
             "2" = "#00BA38",
             "3" = "#619CFF", 
             "4" = "blueviolet",
             "5" = "cyan4",
             "6" = "#E6AB02",
             "7" = "#E36EF6",
             "8" = "bisque4",
             "9" = "coral4",
             "10" = "darkslateblue")
###########################################################
# POPULATION 1 GLOBAL LEVEL CLUSTERING OF GLOBAL VARIABLES
###########################################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = k.samples.HDP[[1]][[iter]][t.samples.HDP[[1]][[iter]]]
}

samples <- samples.thin
k1.rand = getDahl(index[samples],
                  k1.true[t1.true])

k1.rand
best.index = samples[k1.rand$DahlIndex]

library(tidyverse)

g1 = data.frame(x = X.global[[1]][1,], y = X.global[[1]][2,], cluster = factor(index[[best.index]])) %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point() + labs(x = "Gobal variable 1", y = "Global variable 2", title = "Population 1", subtitle = paste0("Adjusted Rand Index = ", round(k1.rand$adj.rand, 4))) + scale_color_manual(values = myvalues) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  )

g1
###########################################################
# POPULATION 2 GLOBAL LEVEL CLUSTERING OF GLOBAL VARIABLES
###########################################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = k.samples.HDP[[2]][[iter]][t.samples.HDP[[2]][[iter]]]
}

samples <- samples.thin 
k2.rand = getDahl(index[samples],
                  k2.true[t2.true])
k2.rand
best.index = samples[k2.rand$DahlIndex]

g2 = data.frame(x = X.global[[2]][1,], y = X.global[[2]][2,], cluster = factor(index[[best.index]])) %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point() + labs(x = "Gobal variable 1", y = "Global variable 2", title = "Population 2", subtitle = paste0("Adjusted Rand Index = ", round(k2.rand$adj.rand,4))) + scale_color_manual(values = myvalues) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  )

g2
###########################################################
# POPULATION 3 GLOBAL LEVEL CLUSTERING OF GLOBAL VARIABLES
###########################################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = k.samples.HDP[[3]][[iter]][t.samples.HDP[[3]][[iter]]]
}

samples <- samples.thin 
k3.rand = getDahl(index[samples],
                  k3.true[t3.true])

k3.rand
best.index = samples[k3.rand$DahlIndex]

g3 = data.frame(x = X.global[[3]][1,], y = X.global[[3]][2,], cluster = factor(index[[best.index]])) %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point() + labs(x = "Gobal variable 1", y = "Global variable 2", title = "Population 3", subtitle = paste0("Adjusted Rand Index = ", round(k3.rand$adj.rand, 4))) + scale_color_manual(values = myvalues) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  )

g3
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)

plot.HDP = gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = textGrob("HDP", gp = gpar(fontface = "bold", fontsize = 20)))
