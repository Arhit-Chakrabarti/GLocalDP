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
X.global = list(X1, X2[-c(1,2), ], X3[-c(1, 2, 3), ]) 
p.global = nrow(X.global[[1]]) # Dimension of global variables
n1 = ncol(X.global[[1]]) # Sample size in population 1
n2 = ncol(X.global[[2]]) # Sample size in population 2
n3 = ncol(X.global[[3]]) # Sample size in population 3

# Save local variables as a list across the populations
X.local = list(X2[c(1, 2), ], X3[c(1,2,3), ]) 
# Save local variables as a list of matrices across the populations
X.local.matrix = list(X2[c(1, 2), ], X3[c(1,2,3), ]) 
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
t.samples <- replicate(n = J, list(replicate(n = num_iter, list(0))))

t.samples[[1]][[1]] = t.start[[1]] # Initialize the MCMC starting value for population 1
t.samples[[2]][[1]] = t.start[[2]] # Initialize the MCMC starting value for population 2
t.samples[[3]][[1]] = t.start[[3]] # Initialize the MCMC starting value for population 3

# Initialize a list to store the MCMC samples of k_jt
k.samples <- replicate(n = J, list(replicate(n = num_iter, list(0))))

k.samples[[1]][[1]] = k.start[[1]] # Initialize the MCMC starting value for population 1
k.samples[[2]][[1]] = k.start[[2]] # Initialize the MCMC starting value for population 2
k.samples[[3]][[1]] = k.start[[3]] # Initialize the MCMC starting value for population 3

# Initialize a list to store the MCMC samples of pi_j
pi.samples <- replicate(n = J, list(replicate(n = num_iter, list())))
pi.samples[[1]][[1]] = pi.start[[1]] # Initialize the MCMC starting value for population 1
pi.samples[[2]][[1]] = pi.start[[2]] # Initialize the MCMC starting value for population 2
pi.samples[[3]][[1]] = pi.start[[3]] # Initialize the MCMC starting value for population 3

# Initialize a list to store the MCMC samples of beta
beta.samples <- list()
beta.samples[[1]] <- beta.start # Initialize the MCMC starting value for beta

# Initialize List to store MCMC samples of global parameter (means and variances)
M.k.samples <- list(); Tau2.k.samples <- list() 
# Initialize List to store MCMC samples of local parameter (means and variances)
M.j.samples <- replicate(n = 2, list(replicate(n = num_iter, list())))
Tau2.j.samples <- replicate(n = 2, list(replicate(n = num_iter, list())))
# Initialize List to store MCMC samples of alpha
alpha.samples <- list(); alpha.samples[[1]] = alpha.start # Starting value of alpha
# Initialize List to store MCMC samples of gamma
m.samples <- list(); m.samples[[1]] = m.start # Starting value of gamma

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
  
  m = replicate(n = J, list(0))
  for(j in 1:J){
    for(l in 1:L){
      m[[j]][l]= sum(t.samples[[j]][[iter - 1]] == l)
    }
  }
  
  
  for(j in 1:J){
    pi.samples[[j]][[iter]] = as.numeric(rdirichlet(n = 1, alpha = m[[j]] + alpha.samples[[iter - 1]]/L))
  }
  
  for(j in 1:J){
    pi.samples[[j]][[iter]] = (pi.samples[[j]][[iter]] + 1e-8)/sum(pi.samples[[j]][[iter]] + 1e-8)
  }
  
  d = replicate(n = J, list(0))
  for(j in 1:J){
    for(k in 1:L){
      d[[j]][k] = sum(k.samples[[j]][[iter - 1]] == k)
      
    } 
  }
  
  d.all = Reduce(`+`, d)
  
  beta.samples[[iter]] = as.numeric(rdirichlet(n = 1, alpha = d.all + m.samples[[iter - 1]]/L))
  
  beta.samples[[iter]] = (beta.samples[[iter]] + 1e-8)/sum(beta.samples[[iter]] + 1e-8)
  
  n.k = replicate(n = J, list(0))
  for(j in 1:J){
    for(k in 1:L){
      n.k[[j]][k] = sum(k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]] == k)
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
      x.k.bar[[k]] = (rowSums(X.global[[1]][, k.samples[[1]][[iter - 1]][t.samples[[1]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[2]][, k.samples[[2]][[iter - 1]][t.samples[[2]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[3]][, k.samples[[3]][[iter - 1]][t.samples[[3]][[iter - 1]]] == k, drop = FALSE]))/n.k.all[k]
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
        S.k.j[,j] = rowSums((X.global[[j]][, k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]] == k, drop = FALSE] - x.k.bar[[k]])^2)
        
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
  
  
  Tau2.k.samples[[iter]] = tau2.k.samples
  M.k.samples[[iter]] = m.k.samples
  
  n.j.t = replicate(n = 3, list(0));
  
  for(t in 1:L){
    for(j in 1:J){
      n.j.t[[j]][t] = sum(t.samples[[j]][[iter - 1]] == t)
    }
    
  }
  
  x.j.bar.local = replicate(n = 3, list())
  
  for(t in 1:L){
    x.j.bar.local[[1]][[t]] = 0
  }
  
  for(t in 1:L){
    if(n.j.t[[2]][t] == 0){
      x.j.bar.local[[2]][[t]] = rep(0, nrow(X.local[[1]]))
    }else{
      x.j.bar.local[[2]][[t]] = rowSums(X.local[[1]][, t.samples[[2]][[iter - 1]] == t, drop = FALSE])/n.j.t[[2]][t]
    }
  }
  
  for(t in 1:L){
    if(n.j.t[[3]][t] == 0){
      x.j.bar.local[[3]][[t]] = rep(0, nrow(X.local[[2]]))
    }else{
      x.j.bar.local[[3]][[t]] = rowSums(X.local[[2]][, t.samples[[3]][[iter - 1]] == t, drop = FALSE])/n.j.t[[3]][t]
    }
  }
  
  S.2.local = matrix(0, nrow = nrow(X.local[[1]]), L)
  Z.2.local = matrix(0, nrow = nrow(X.local[[1]]), L)
  alpha.2.hat = matrix(0, nrow = nrow(X.local[[1]]), L); beta.2.hat =matrix(0, nrow = nrow(X.local[[1]]), L)
  lambda.2.hat = matrix(0, nrow = nrow(X.local[[1]]), L); m.2.hat = matrix(0, nrow = nrow(X.local[[1]]), L)
  
  for(t in 1:L){
    if(n.j.t[[2]][t] != 0){
      S.2.local[, t] = rowSums((X.local.matrix[[1]][, t.samples[[2]][[iter - 1]] == t, drop = FALSE] - x.j.bar.local[[2]][[t]])^2)
      Z.2.local[, t] = ((lambda.local[2] * n.j.t[[2]][t])/(lambda.local[2] + n.j.t[[2]][t])) * rowSums(matrix((x.j.bar.local[[2]][[t]] - m.j0[[2]])^2, ncol = 1))
      alpha.2.hat[, t] = alpha.j0 + 0.5 * n.j.t[[2]][t] 
      beta.2.hat[, t] = beta.j0 + 0.5 * (S.2.local[, t] + Z.2.local[, t])
      lambda.2.hat[, t] = lambda.local[1] + n.j.t[[2]][t]
      m.2.hat[, t] = ((lambda.local[2] * m.j0[[2]]) + (n.j.t[[2]][t] * x.j.bar.local[[2]][[t]]))/lambda.2.hat[, t]
    }else{
      alpha.2.hat[, t] = alpha.j0 
      beta.2.hat[, t] = beta.j0 
      lambda.2.hat[, t] = lambda.local[1] 
      m.2.hat[, t] = ((lambda.local[2] * m.j0[[2]]) + (n.j.t[[2]][t] * x.j.bar.local[[2]][[t]]))/lambda.2.hat[, t]
    }
  }
  
  S.3.local = matrix(0, nrow = nrow(X.local[[2]]), L)
  Z.3.local = matrix(0, nrow = nrow(X.local[[2]]), L)
  alpha.3.hat = matrix(0, nrow = nrow(X.local[[2]]), L); beta.3.hat =matrix(0, nrow = nrow(X.local[[2]]), L)
  lambda.3.hat = matrix(0, nrow = nrow(X.local[[2]]), L); m.3.hat = matrix(0, nrow = nrow(X.local[[2]]), L)
  
  for(t in 1:L){
    if(n.j.t[[3]][t] != 0){
      S.3.local[, t] = rowSums((X.local.matrix[[2]][, t.samples[[3]][[iter - 1]] == t, drop = FALSE] - x.j.bar.local[[3]][[t]])^2)
      Z.3.local[, t] = ((lambda.local[3] * n.j.t[[3]][t])/(lambda.local[3] + n.j.t[[3]][t])) * rowSums(matrix((x.j.bar.local[[3]][[t]] - m.j0[[3]])^2, ncol = 1))
      alpha.3.hat[, t] = alpha.j0 + 0.5 * n.j.t[[3]][t] 
      beta.3.hat[, t] = beta.j0 + 0.5 * (S.3.local[, t] + Z.3.local[, t])
      lambda.3.hat[, t] = lambda.local[1] + n.j.t[[3]][t]
      m.3.hat[, t] = ((lambda.local[3] * m.j0[[3]]) + (n.j.t[[3]][t] * x.j.bar.local[[3]][[t]]))/lambda.3.hat[, t]
    }else{
      alpha.3.hat[, t] = alpha.j0
      beta.3.hat[, t] = beta.j0
      lambda.3.hat[, t] = lambda.local[1]
      m.3.hat[, t] = ((lambda.local[3] * m.j0[[3]]) + (n.j.t[[3]][t] * x.j.bar.local[[3]][[t]]))/lambda.3.hat[, t]
    }
  }
  
  m.j.samples <- list(matrix(0, nrow = nrow(X.local.matrix[[1]]), ncol = L),
                      matrix(0, nrow = nrow(X.local.matrix[[2]]), ncol = L))
  
  tau2.j.samples <- list(matrix(0, nrow = nrow(X.local.matrix[[1]]), ncol = L),
                         matrix(0, nrow = nrow(X.local.matrix[[2]]), ncol = L))
  
  for(i in 1:nrow(alpha.2.hat)){
    tau2.j.samples[[1]][i, ] = rinvgamma(n = L, alpha = alpha.2.hat[i, ], beta = beta.2.hat[i, ])
    m.j.samples[[1]][i, ] = rnorm(n = L, mean = m.2.hat[i, ], sd = sqrt(tau2.j.samples[[2]][i,] /lambda.2.hat[i,]))
  }
  for(i in 1:nrow(alpha.3.hat)){
    tau2.j.samples[[2]][i, ] = rinvgamma(n = L, alpha = alpha.3.hat[i, ], beta = beta.3.hat[i, ])
    m.j.samples[[2]][i, ] = rnorm(n = L, mean = m.3.hat[i, ], sd = sqrt(tau2.j.samples[[2]][i,] /lambda.3.hat[i,]))
  }
  
  for(j in 1:2){
    M.j.samples[[j]][[iter]] = m.j.samples[[j]]
    Tau2.j.samples[[j]][[iter]] = tau2.j.samples[[j]]
  }
  
  # Calculate the probabilities of table indices and sample
  n = c(n1, n2, n3)
  
  # Define an array of global cluster variance samples 
  Tau2.k.samples.array = array(0, dim = c(p.global, p.global, L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:ncol(Tau2.k.samples[[iter]])){
    Tau2.k.samples.array[,,j] = diag(Tau2.k.samples[[iter]][,j], nrow = p.global)
  }
  # Define an array of local cluster variance samples (for population 2) which is 2-variate Normal
  Tau2.j.samples.array2 = array(0, dim = c(nrow(X.local[[1]]), nrow(X.local[[1]]), L))
  for(j in 1:ncol(Tau2.j.samples[[1]][[iter]])){
    Tau2.j.samples.array2[,,j] = diag(Tau2.j.samples[[1]][[iter]][, j], nrow = nrow(X.local[[1]]))
  }
  # Define an array of local cluster variance samples (for population 3) which is 3-variate Normal
  Tau2.j.samples.array3 = array(0, dim = c(nrow(X.local[[2]]), nrow(X.local[[2]]), L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:ncol(Tau2.j.samples[[2]][[iter]])){
    Tau2.j.samples.array3[,,j] = diag(Tau2.j.samples[[2]][[iter]][, j], nrow = nrow(X.local[[2]]))
  }
  Tau2.j.samples.array.list <- list(Tau2.j.samples.array2, 
                                    Tau2.j.samples.array3)
  t.prob = replicate(n = J, list())
  for(j in 1:J){
    t.prob[[j]] = matrix(0, nrow = n[j], ncol = L)
  }
  
  t.prob[[1]] <- prob_exponent_no_local2(pi_samples = pi.samples[[1]][[iter]], 
                                         phi_samples = matrix(m.k.samples, nrow = nrow(X.global[[1]]), L),
                                         X_global = t(X.global[[1]]), 
                                         k = k.samples[[1]][[iter - 1]],
                                         Sigma_G = Tau2.k.samples.array)
  
  for(j in 2:J){
    t.prob[[j]] = prob_exponent_mv_local(pi_samples = pi.samples[[j]][[iter]], 
                                         X_local = t(X.local[[j-1]]), 
                                         psi_samples = matrix(M.j.samples[[j-1]][[iter]], nrow = nrow(X.local[[j-1]]), ncol = L),
                                         Sigma = Tau2.j.samples.array.list[[j-1]],
                                         phi_samples = matrix(m.k.samples, nrow = nrow(X.global[[j]]), L),
                                         X_global = t(X.global[[j]]), k = k.samples[[j]][[iter - 1]], 
                                         Sigma_G = Tau2.k.samples.array)
  }
  
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){
    t.prob[[j]] = calc_probability_log_sum_exp_normalized(t.prob[[j]])
  }
  
  for(j in 1:J){
    t.samples[[j]][[iter]] = sample_my(t.prob[[j]])
  }
  
  
  # Calculate the probabilities of dish indices and sample
  prob = array(0, dim = c(L, L, J))
  
  for(j in 1:J){
    prob[, , j] = prob_exponent_mv_global(beta_samples = beta.samples[[iter]], 
                                          X_global = t(X.global[[j]]), 
                                          phi_samples = matrix(m.k.samples, nrow = nrow(X.global[[j]]), L), 
                                          Sigma_G = Tau2.k.samples.array,
                                          t_samples = t.samples[[j]][[iter]])}
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){prob[, ,j] = calc_probability_log_sum_exp_normalized(prob[, ,j])}
  
  for(j in 1:J){
    k.samples[[j]][[iter]] = sample_my(prob[, , j])}
  
  sampled.alpha = sample_alpha(ycurrent = alpha.samples[[iter - 1]], pars = list(a0 = 0.1, b0 = 0.1,
                                                                                 pi1 = pi.samples[[1]][[iter]],
                                                                                 pi2 = pi.samples[[2]][[iter]],
                                                                                 pi3 = pi.samples[[3]][[iter]], J = 3))
  alpha.samples[[iter]] = sampled.alpha$out
  
  sampled.m = sample_m(ycurrent = m.samples[[iter - 1]], pars = list(a1 = 0.1, b1 = 0.1, beta = beta.samples[[iter]]))
  m.samples[[iter]] = sampled.m$out
  
}
time.end = Sys.time()
time.end - time.start

################################################################################
# CALCULATE LOG-POSTERIOR TO PLOT TRACEPLOT AND ACF OF LOG-POSTERIOR
################################################################################
library(mvtnorm)
ll = 0

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
  for(j in 1:ncol(Tau2.k.samples[[iter]])){
    Tau2.k.samples.array[,,j] = diag(Tau2.k.samples[[iter]][, j], nrow = p.global)
  }
  # Define an array of local cluster variance samples (for population 2) which is 2-variate Normal
  Tau2.j.samples.array2 = array(0, dim = c(nrow(X.local[[1]]), nrow(X.local[[1]]), L))
  for(j in 1:ncol(Tau2.j.samples[[1]][[iter]])){
    Tau2.j.samples.array2[,,j] = diag(Tau2.j.samples[[1]][[iter]][, j], nrow = nrow(X.local[[1]]))
  }
  # Define an array of local cluster variance samples (for population 3) which is 3-variate Normal
  Tau2.j.samples.array3 = array(0, dim = c(nrow(X.local[[2]]), nrow(X.local[[2]]), L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:ncol(Tau2.j.samples[[2]][[iter]])){
    Tau2.j.samples.array3[,,j] = diag(Tau2.j.samples[[2]][[iter]][, j], nrow = nrow(X.local[[2]]))
  }
  # To calculate the LL Psi-samples are replaced by the M.j.samples (Mean of Local NIG posterior)
  # Phi-samples are replaced by M.k.samples (Mean of Global NIG posterior)
  
  logsum = logll_no_local2(pi_samples = pi.samples[[1]][[iter]],
                           t_samples = t.samples[[1]][[iter]],
                           phi_samples = matrix(M.k.samples[[iter]], nrow = p.global, ncol = L),
                           X_global = X.global[[1]], 
                           k_samples = k.samples[[1]][[iter]], 
                           Sigma_G = Tau2.k.samples.array) + 
    logll_mv_local2(pi.samples[[2]][[iter]], X.local[[1]],
                    matrix(M.j.samples[[1]][[iter]], nrow = 2, ncol = L),  
                    Sigma = Tau2.j.samples.array2,
                    t.samples[[2]][[iter]], 
                    matrix(M.k.samples[[iter]], nrow = p.global, ncol = L), 
                    X.global[[2]], k.samples[[2]][[iter]], 
                    Sigma_G = Tau2.k.samples.array) + 
    logll_mv_local2(pi.samples[[3]][[iter]], X.local[[2]],
                    matrix(M.j.samples[[2]][[iter]], nrow = 3, ncol = L),  
                    Sigma = Tau2.j.samples.array3,
                    t.samples[[3]][[iter]], 
                    matrix(M.k.samples[[iter]], nrow = p.global, ncol = L), 
                    X.global[[3]], k.samples[[3]][[iter]], 
                    Sigma_G = Tau2.k.samples.array)
  
  m1.tilde = count_my(t.samples[[1]][[iter]], L)
  m2.tilde = count_my(t.samples[[2]][[iter]], L)
  m3.tilde = count_my(t.samples[[3]][[iter]], L)
  
  logsum = logsum + sum(m1.tilde * log(pi.samples[[1]][[iter]])) + sum(m2.tilde * log(pi.samples[[2]][[iter]])) + sum(m3.tilde * log(pi.samples[[3]][[iter]]))
  
  m.tilde = count_my(k.samples[[1]][[iter]], L) + count_my(k.samples[[2]][[iter]], L) + count_my(k.samples[[3]][[iter]], L) 
  
  logsum = logsum + sum(m.tilde * log(beta.samples[[iter]]))
  
  
  Local_density2 = 0
  for(l in 1:ncol(M.j.samples[[1]][[iter]])){
    Local_density2[l] <- dmvnorm(M.j.samples[[1]][[iter]][, l], mean = m.j0[[2]], sigma = diag(1/(lambda.local[2]) * Tau2.j.samples[[1]][[iter]][, l], length(m.j0[[2]])), log = TRUE)
  }
  
  Local_density2_sum = sum(Local_density2) + sum(dinvgamma(Tau2.j.samples[[2]][[iter]][1, ], alpha = alpha0, beta = beta0, log = TRUE)) + sum(dinvgamma(Tau2.j.samples[[2]][[iter]][2, ], alpha = alpha0, beta = beta0, log = TRUE))
  
  Local_density3 = 0
  for(l in 1:ncol(M.j.samples[[2]][[iter]])){
    Local_density3[l] <- dmvnorm(M.j.samples[[2]][[iter]][, l], mean = m.j0[[3]], sigma = diag(1/(lambda.local[3]) * Tau2.j.samples[[2]][[iter]][, l], length(m.j0[[3]])), log = TRUE)
  }
  
  Local_density3_sum = sum(Local_density3) + sum(dinvgamma(Tau2.j.samples[[2]][[iter]][1, ], alpha = alpha0, beta = beta0, log = TRUE)) + sum(dinvgamma(Tau2.j.samples[[2]][[iter]][2, ], alpha = alpha0, beta = beta0, log = TRUE)) + sum(dinvgamma(Tau2.j.samples[[2]][[iter]][3, ], alpha = alpha0, beta = beta0, log = TRUE))
  
  ## Component 1
  Local_density_sum = Local_density2_sum + Local_density3_sum
  
  Global_density = 0
  for(l in 1:ncol(M.k.samples[[iter]])){
    Global_density[l] <- dmvnorm(M.k.samples[[iter]][, l], mean = m0, sigma = diag(1/(lambda.global) * Tau2.k.samples[[iter]][, l], length(m0)), log = TRUE)
  }
  ## Component 2
  Global_density_sum = sum(Global_density) + sum(dinvgamma(Tau2.k.samples[[iter]][1, ], alpha = alpha0, beta = beta0, log = TRUE)) + sum(dinvgamma(Tau2.k.samples[[iter]][2, ], alpha = alpha0, beta = beta0, log = TRUE))
  ## Component 3
  pi_density_sum = ddirichlet(pi.samples[[1]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE) + ddirichlet(pi.samples[[2]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE) + ddirichlet(pi.samples[[3]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE)
  ## Component 4
  beta_density_sum = ddirichlet(beta.samples[[iter]], alpha = rep(m.samples[[iter]]/L, L), log = TRUE)
  
  ll[iter] = logsum + Local_density_sum + Global_density_sum + pi_density_sum + beta_density_sum + dgamma(alpha.samples[[iter]], shape = 0.1, rate = 0.1, log = TRUE) + dgamma(m.samples[[iter]], shape = 0.1, rate = 0.1, log = TRUE) 
}

################################################################################
# PLOT TRACEPLOT OF LOG-POSTERIOR AND ACF
################################################################################
burn = 25000 # Burn-in 
samples.thin = seq((burn + 1), (num_iter), by = 25)

log_like <- ll[samples.thin]

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
alpha_plot <- data.frame(alpha = unlist(alpha.samples[samples.thin]),
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


gamma_plot <- data.frame(m = unlist(m.samples[samples.thin]),
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
x.limit.lower <- min(X1[1, ], X2[3, ], X3[4, ])
x.limit.upper <- max(X1[1, ], X2[3, ], X3[4, ])

y.limit.lower <- min(X1[2, ], X2[4, ], X3[5, ])
y.limit.upper <- max(X1[2, ], X2[4, ], X3[5, ])


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
#######################################################
# GLOBAL LEVEL CLUSTERING
#######################################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = c(k.samples[[1]][[iter]][t.samples[[1]][[iter]]],
                    k.samples[[2]][[iter]][t.samples[[2]][[iter]]],
                    k.samples[[3]][[iter]][t.samples[[3]][[iter]]])
}


samples <- samples.thin
k.rand = getDahl(index[samples],
                 c(k1.true[t1.true],
                   k2.true[t2.true],
                   k3.true[t3.true]) 
)

best.index = samples[k.rand$DahlIndex]

singleton = as.numeric(names(which(table(index[[best.index]]) <= 0.0 * length(index[[best.index]]))))

singleton.index = which(index[[best.index]] %in% singleton)

z.estimated = index[[best.index]]

cluster.global <- data.frame(x = c(X.global[[1]][1, ],
                                   X.global[[2]][1, ],
                                   X.global[[3]][1, ]),
                             y = c(X.global[[1]][2, ],
                                   X.global[[2]][2, ],
                                   X.global[[3]][2, ]),
                             cluster = factor(z.estimated),
                             cancer = factor(c(rep("Population 1", length(X.global[[1]][1,])),
                                               rep("Population 2", length(X.global[[2]][1,])),
                                               rep("Population 3", length(X.global[[3]][1,])))))

if(length(singleton.index) > 0){
  z.estimated = z.estimated[-singleton.index_mmclust2]
  cluster.global = cluster.global[-singleton.index, ]
}else{
  z.estimated = z.estimated
  cluster.global = cluster.global
}

library(tidyverse)

cluster.global$population = factor(c(rep(paste0("Population 1\n", "ARI = ", round(aricode::ARI(cluster.global %>% filter(cancer == "Population 1") %>% pull(cluster),
                                                                                               k1.true[t1.true]), 3)), length(X.global[[1]][1,])),
                                     rep(paste0("Population 2\n", "ARI = ", round(aricode::ARI(cluster.global %>% filter(cancer == "Population 2") %>% pull(cluster),
                                                                                               k2.true[t2.true]), 3)), length(X.global[[2]][1,])),
                                     rep(paste0("Population 3\n", "ARI = ", round(aricode::ARI(cluster.global %>% filter(cancer == "Population 3") %>% pull(cluster),
                                                                                               k3.true[t3.true]), 3)), length(X.global[[3]][1,]))))

plot.globalLS <- cluster.global %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.8) + facet_wrap(~population, ncol = 1)+ labs(x = "Global variable 1", y = "Global variable 2", title = "Global level clustering")+ 
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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

#######################################################
# LOCAL LEVEL CLUSTERING
#######################################################
##########################################
# Population 2
##########################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = c(t.samples[[2]][[iter]])
}

samples <- samples.thin
k2.rand = getDahl(index[samples],
                  t2.true)
best.index = samples[k2.rand$DahlIndex]
singleton_2.2 = as.numeric(names(which(table(index[[best.index]]) <= 0.0 * length(index[[best.index]]))))
singleton.index_2.2 = which(index[[best.index]] %in% singleton_2.2)
z.estimated_2.2 = index[[best.index]]

cluster.local_2.2 <- data.frame(x = c(X.global[[2]][1, ]),
                                y = c(X.global[[2]][2, ]), 
                                cluster = factor(c(z.estimated_2.2)),
                                cancer = factor(c(rep("Population 2", length(X.global[[2]][1,])))))

if(length(singleton.index_2.2) > 0){
  z.estimated_2.2 = z.estimated_2.2[-singleton.index_2.2]
  cluster.local_2.2 = cluster.local_2.2[-singleton.index_2.2, ]
}else{
  z.estimated_2.2 = z.estimated_2.2
  cluster.local_2.2 = cluster.local_2.2
}

##########################################
# Population 3
##########################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = c(t.samples[[3]][[iter]])
}

samples <- samples.thin
k3.rand = getDahl(index[samples],
                  t3.true)

best.index = samples[k3.rand$DahlIndex]
singleton_2.3 = as.numeric(names(which(table(index[[best.index]]) <= 0.0 * length(index[[best.index]]))))
singleton.index_2.3 = which(index[[best.index]] %in% singleton_2.3)
z.estimated_2.3 = index[[best.index]]

cluster.local_2.3 <- data.frame(x = c(X.global[[3]][1, ]),
                                y = c(X.global[[3]][2, ]), 
                                cluster = factor(c(z.estimated_2.3)),
                                cancer = factor(c(rep("Population 3", length(X.global[[3]][1,])))))

if(length(singleton.index_2.3) > 0){
  z.estimated_2.3 = z.estimated_2.3[-singleton.index_2.3]
  cluster.local_2.3 = cluster.local_2.3[-singleton.index_2.3, ]
}else{
  z.estimated_2.3 = z.estimated_2.3
  cluster.local_2.3 = cluster.local_2.3
}

cluster.local = bind_rows(cluster.local_2.2, cluster.local_2.3)

cluster.local$population = factor(c(rep(paste0("Population 2\n", "ARI = ", round(aricode::ARI(cluster.local %>% filter(cancer == "Population 2") %>% pull(cluster),
                                                                                              t2.true), 3)), length(X.global[[2]][1,])),
                                    rep(paste0("Population 3\n", "ARI = ", round(aricode::ARI(cluster.local %>% filter(cancer == "Population 3") %>% pull(cluster),
                                                                                              t3.true), 3)), length(X.global[[3]][1,]))))

cluster.local %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3) + facet_wrap(~population)+ labs(x = "Global variable 1", y = "Global variable 2", title = "Local level clustering") + 
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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


#################################################
# RELABELLING LOCAL CLUSTERS 
#################################################
#################################################
# POPULATION 2
#################################################
cluster.global.LS.pop2 = cluster.global %>% filter(cancer == "Population 2") 
cluster.local.LS.pop2 = cluster.local %>% filter(cancer == "Population 2") 
cluster.local.LS.pop2$cluster_org = cluster.local.LS.pop2$cluster

charaters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")

for(l in 1:length(unique(cluster.global.LS.pop2$cluster))){
  if(length(singleton.index_2.2) > 0){
    cluster.global.LS.pop2.subset = cluster.global.LS.pop2$cluster[-singleton.index_2.2]
  }else{
    cluster.global.LS.pop2.subset = cluster.global.LS.pop2$cluster
  }
  
  cluster.index = as.numeric(as.character(unique(cluster.global.LS.pop2.subset)[l]))
  
  data.index = which(cluster.global.LS.pop2.subset == cluster.index)
  
  labels.local = cluster.local.LS.pop2$cluster[which(cluster.global.LS.pop2.subset == cluster.index)]
  labels.local.unique =  as.character(unique(cluster.local.LS.pop2$cluster[which(cluster.global.LS.pop2.subset == cluster.index)]))
  
  if(length(labels.local.unique) > 1){
    labels.local.update = NULL
    for(i in 1:length(labels.local.unique)){
      new.label = paste0(as.character(cluster.index),charaters[i])
      labels.local.update = c(labels.local.update, new.label)
    }
  }else{
    labels.local.update = as.character(cluster.index)
  }
  
  labels.local = factor(labels.local, labels = labels.local.update)
  cluster.local.LS.pop2$cluster <- as.character(cluster.local.LS.pop2$cluster)
  cluster.local.LS.pop2$cluster[data.index] <- as.character(labels.local)
}

cluster.local.LS.pop2$cluster <- factor(cluster.local.LS.pop2$cluster)

cluster.local.LS.pop2$cluster

myvalues_local.LS2 = c("1" = "#F8766D",
                       "2" = "#00BA38",
                       "3" = "#619CFF", 
                       "4" = "blueviolet",
                       "5" = "cyan4",
                       "6" = "#E6AB02",
                       "7" = "#E36EF6",
                       "8" = "bisque4",
                       "9" = "coral4",
                       "10" = "darkslateblue",
                       
                       "5a" = "lightseagreen",
                       "5b" = "yellow3",
                       "5c" = "bisque",
                       
                       "9a" = "hotpink3",
                       "9b" = "sienna2",
                       "9c" = "dodgerblue4",
                       "9d" = "dodgerblue")

#################################################
# POPULATION 3
#################################################
cluster.global.LS.pop3 = cluster.global %>% filter(cancer == "Population 3") 
cluster.local.LS.pop3 = cluster.local %>% filter(cancer == "Population 3") 
cluster.local.LS.pop3$cluster_org = cluster.local.LS.pop3$cluster

charaters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")

for(l in 1:length(unique(cluster.global.LS.pop3$cluster))){
  if(length(singleton.index_2.3) > 0){
    cluster.global.LS.pop3.subset = cluster.global.LS.pop3$cluster[-singleton.index_2.3]
  }else{
    cluster.global.LS.pop3.subset = cluster.global.LS.pop3$cluster
  }
  
  cluster.index = as.numeric(as.character(unique(cluster.global.LS.pop3.subset)[l]))
  
  data.index = which(cluster.global.LS.pop3.subset == cluster.index)
  
  labels.local = cluster.local.LS.pop3$cluster[which(cluster.global.LS.pop3.subset == cluster.index)]
  labels.local.unique =  as.character(unique(cluster.local.LS.pop3$cluster[which(cluster.global.LS.pop3.subset == cluster.index)]))
  
  if(length(labels.local.unique) > 1){
    labels.local.update = NULL
    for(i in 1:length(labels.local.unique)){
      new.label = paste0(as.character(cluster.index),charaters[i])
      labels.local.update = c(labels.local.update, new.label)
    }
  }else{
    labels.local.update = as.character(cluster.index)
  }
  
  labels.local = factor(labels.local, labels = labels.local.update)
  cluster.local.LS.pop3$cluster <- as.character(cluster.local.LS.pop3$cluster)
  cluster.local.LS.pop3$cluster[data.index] <- as.character(labels.local)
}

cluster.local.LS.pop3$cluster <- factor(cluster.local.LS.pop3$cluster)

cluster.local.LS.pop3$cluster

myvalues_local.LS3 = c("1" = "#F8766D",
                       "2" = "#00BA38",
                       "3" = "#619CFF", 
                       "4" = "blueviolet",
                       "5" = "cyan4",
                       "6" = "#E6AB02",
                       "7" = "#E36EF6",
                       "8" = "bisque4",
                       "9" = "coral4",
                       "10" = "darkslateblue",
                       
                       "5a" = "lightseagreen",
                       "5b" = "yellow3",
                       "5c" = "bisque")

cluster.local.LS.combined = bind_rows(cluster.local.LS.pop2,
                                      cluster.local.LS.pop3)

cluster.local.LS.combined$cluster
myvalues_local.LS.combined = c(myvalues_local.LS2, myvalues_local.LS3)

plot.localLS <- cluster.local.LS.combined %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.8) + facet_wrap(~population, ncol = 1)+ labs(x = "Global variable 1", y = "Global variable 2", title = "Local level clustering") +  scale_color_manual(values = myvalues_local.LS.combined) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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


if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
gridExtra::grid.arrange(plot.globalLS, plot.localLS, ncol = 2)

#########################################################
# LOCAL LEVEL CLUSTERING OF LOCAL VARIABLES
#########################################################
##########################################
# POPULATION 2
##########################################
ll2 = data.frame(x = X.local[[1]][1, ], y = X.local[[1]][2, ], cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Population 2") %>% pull(cluster))) %>% 
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.8) + labs(x = "Local variable 1", y = "Local variable 2", title = "Population 2", subtitle = paste0("ARI = ", round(k2.rand$adj.rand, 4)))  + scale_color_manual(values = myvalues_local.LS.combined)  +
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

##########################################
# POPULATION 3
##########################################
ll3.1 = data.frame(x = X.local[[2]][1, ], y = X.local[[2]][2, ], cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Population 3") %>% pull(cluster))) %>% 
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.8) + labs(x = "Local variable 1", y = "Local variable 2", title = "Population 3", subtitle = paste0("ARI = ", round(k3.rand$adj.rand, 4)))  + scale_color_manual(values = myvalues_local.LS.combined)  +
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

ll3.2 = data.frame(x = X.local[[2]][2, ], y = X.local[[2]][3, ], cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Population 3") %>% pull(cluster))) %>% 
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.8) + labs(x = "Local variable 2", y = "Local variable 3", title = "Population 3", subtitle = paste0("ARI = ", round(k3.rand$adj.rand, 4)))  + scale_color_manual(values = myvalues_local.LS.combined)  +
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
gridExtra::grid.arrange(ll2, ll3.1, ll3.2, ncol = 3)

