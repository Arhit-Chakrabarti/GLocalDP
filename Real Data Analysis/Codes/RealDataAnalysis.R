library(data.table)
# # READING THE PHENOTYPE DATA
colon_cov = fread("./Real Data Analysis/Data/Colon.gz")
rectal_cov = fread("./Real Data Analysis/Data/Rectal.gz")
stomach_cov = fread("./Real Data Analysis/Data/Stomach.gz")
eoso_cov = fread("./Real Data Analysis/Data/Eoso.gz")


# READING THE GENE-EXPRESSION DATA
colon = fread("./Real Data Analysis/Data/Colon_counts.tsv.gz")
rectal = fread("./Real Data Analysis/Data/Rectal_counts.tsv.gz")
stomach = fread("./Real Data Analysis/Data/Stomach_counts.tsv.gz")
eoso = fread("./Real Data Analysis/Data/Eoso_counts.tsv.gz")

library(tidyverse)
common.genes <- as.vector(intersect(intersect(intersect(colon[c(1:60483), 1],
                                                        rectal[c(1:60483), 1]),
                                              stomach[c(1:60483), 1]),
                                    eoso[c(1:60483), 1]))

attributes(common.genes) <- NULL
common.genes = unlist(common.genes)

colon_data <- t(colon %>% filter(Ensembl_ID %in% common.genes) %>%
                  select(-Ensembl_ID))

rectal_data <- t(rectal %>% filter(Ensembl_ID %in% common.genes) %>%
                   select(-Ensembl_ID))

stomach_data <- t(stomach %>% filter(Ensembl_ID %in% common.genes) %>%
                    select(-Ensembl_ID))
eoso_data <- t(eoso %>% filter(Ensembl_ID %in% common.genes) %>%
                 select(-Ensembl_ID))

all_data_merged <- rbind(colon_data,
                         rectal_data,
                         stomach_data,
                         eoso_data)

library(uwot)
library(Matrix)
# CHANGE N_COMPONENTS TO 3 AND 5 
umap.all = uwot::umap(all_data_merged, n_neighbors = 100,
                      n_components = 2,
                      metric = "euclidean")

library(tidyverse)
# Segreggate the UMAP embeddings according to the colon, rectal, stomach, and esophageal cancers
umap.colon = umap.all[which(rownames(umap.all) %in% colon_cov$submitter_id.samples),]      # Colon Cancer
umap.rectal = umap.all[which(rownames(umap.all) %in% rectal_cov$submitter_id.samples),]    # Rectal Cancer   
umap.stomach = umap.all[which(rownames(umap.all) %in% stomach_cov$submitter_id.samples),]  # Stomach Cancer
umap.eoso = umap.all[which(rownames(umap.all) %in% eoso_cov$submitter_id.samples),]        # Esophageal Cancer

# Plot the Data
x.limit.lower <- min(umap.all[, 1])
x.limit.upper <- max(umap.all[, 1])
y.limit.lower <- min(umap.all[, 2])
y.limit.upper <- max(umap.all[, 2])

library(latex2exp)
plot1 = data.frame(umap.colon) %>% ggplot(aes(x = X1, y = X2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + ggtitle("Colon cancer") + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + theme_classic() +  
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
plot2 = data.frame(umap.rectal) %>% ggplot(aes(x = X1, y = X2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + ggtitle("Rectal cancer") + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + theme_classic() +  
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
plot3 = data.frame(umap.stomach) %>% ggplot(aes(x = X1, y = X2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + ggtitle("Stomach cancer") + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + theme_classic() +  
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
plot4 = data.frame(umap.eoso) %>% ggplot(aes(x = X1, y = X2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + ggtitle("Esophageal cancer") + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + theme_classic() +  
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
gridExtra::grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

####################################################################################################
# FINAL DATA GENERATION BEFORE RUNNING OUR SAMPLER
####################################################################################################
stomach_UMAP <- data.frame(submitter_id.samples = rownames(umap.stomach), umap.stomach) %>% na.omit()
rownames(umap.stomach) <- NULL

colon_UMAP <- data.frame(submitter_id.samples = rownames(umap.colon), umap.colon) %>% na.omit() 
rownames(colon_UMAP) <- NULL

rectal_UMAP <- data.frame(submitter_id.samples = rownames(umap.rectal), umap.rectal) %>% na.omit() 
rownames(rectal_UMAP) <- NULL

eoso_UMAP <- data.frame(submitter_id.samples = rownames(umap.eoso), umap.eoso) %>% na.omit() 
rownames(eoso_UMAP) <- NULL

# MERGE THE UMAP EMBEDDINGS WITH THE COVARIATE DATA
colon_merged <- left_join(colon_UMAP, colon_cov, by = "submitter_id.samples")
rectal_merged <- left_join(rectal_UMAP, rectal_cov, by = "submitter_id.samples")
stomach_merged <- left_join(stomach_UMAP, stomach_cov, by = "submitter_id.samples")
eoso_merged <- left_join(eoso_UMAP, eoso_cov, by = "submitter_id.samples")

# STOMACH CANCER ONLY CONTAINS THE GLOBAL VARIABLES AND NO LOCAL VARIABLE
stomach_full <- stomach_merged %>% dplyr::select(X1,
                                                 X2) %>% na.omit()

# RECTAL CANCER CONTAINS THE GLOBAL VARIABLES AND CEA AS A LOCAL VARIABLE
rectal_full <- rectal_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                               X1,
                                               X2) %>% na.omit()
# COLON CANCER CONTAINS THE GLOBAL VARIABLES, BMI AND CEA AS TWO LOCAL VARIABLES
colon_full <- colon_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                             bmi.exposures,
                                             X1,
                                             X2) %>% na.omit()
# ESOPHAGEAL CANCER CONTAINS THE GLOBAL VARIABLES AND CIGAREETES PER DAY AS A LOCAL VARIABLE
eoso_full <- eoso_merged %>% dplyr::select(cigarettes_per_day.exposures,
                                           X1,
                                           X2) %>% na.omit()
colnames(stomach_full) <- NULL
colnames(rectal_full) <- NULL
colnames(colon_full) <- NULL
colnames(eoso_full) <- NULL

# TAKE AS X1 THE DATA FOR STOMACH CANCER
X1 <- t(as.matrix(stomach_full))
# TAKE AS X2 THE DATA FOR RECTAL CANCER
X2 <- t(as.matrix(rectal_full))
# TAKE AS X3 THE DATA FOR COLON CANCER
X3 <- t(as.matrix(colon_full))
# TAKE AS X4 THE DATA FOR ESOPHAGEAL CANCER
X4 <- t(as.matrix(eoso_full))
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
# @pars argument should be a named list with a0, b0, J, pi1, pi2, pi3, pi4
# @ycurrent corresponds to the previous value of alpha
sample_alpha <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  pi1 = as.numeric(pars$pi1); pi2 = as.numeric(pars$pi2); pi3 = as.numeric(pars$pi3); pi4 = as.numeric(pars$pi4)
  a0 = pars$a0; b0 = pars$b0
  J = pars$J
  L = length(pi1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = - b0 * ycand + (a0 - 1)*log(ycand) + sum(((ycand/L) - 1) * log(pi1)) + sum((ycand/L - 1) * log(pi2)) + sum((ycand/L - 1) * log(pi3)) + sum((ycand/L - 1) * log(pi4)) + J * lgamma(ycand) - (J * L * lgamma(ycand/L)) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = - b0 * ycurrent + (a0 - 1)*log(ycurrent) + sum(((ycurrent/L) - 1) * log(pi1)) + sum((ycurrent/L - 1) * log(pi2)) + sum((ycurrent/L - 1) * log(pi3)) + sum((ycurrent/L - 1) * log(pi4)) + J * lgamma(ycurrent) - (J * L * lgamma(ycurrent/L))  + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
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
# Set the parameters of the GLocal DP
###############################################################################
L = 20 # Truncation level of GLocal DP 
J = 4 # Number of populations
p_G = nrow(X1) # Number of global variables

## Specify the hyper-priors
alpha0 = 0.1; beta0 = 0.1; m0 = rep(0, p_G) # Prior hyper-parameters corresponding to global parameters
alpha.j0 = 0.1; beta.j0 = 0.1; m.j0 = list(0, c(0, 0), 0) # Prior hyper-parameters corresponding to local parameters

# MCMC initial values of parameters
alpha.start = 1 # Starting value of concentration parameter alpha
m.start = 1 # Starting value of concentration parameter gamma

lambda.global = 1 # Global variable prior-precision
lambda.local = rep(1, J) # Local variable prior-precision acoss J populations

library(extraDistr)
library(MASS)

beta.start = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting beta
pi.start = replicate(n = J, list())
pi.start[[1]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi1
pi.start[[2]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi2
pi.start[[3]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi3
pi.start[[4]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi4

n1 = ncol(X1); n2 = ncol(X2); n3 = ncol(X3); n4 = ncol(X4) # Define the sample sizes for the four populations

# Initialization of local-level indicators for the different populations
t.start = replicate(n = J, list())
t.start[[1]] = sample(1:L, size = n1, prob = pi.start[[1]], replace = TRUE)
t.start[[2]] = sample(1:L, size = n2, prob = pi.start[[2]], replace = TRUE)
t.start[[3]] = sample(1:L, size = n3, prob = pi.start[[3]], replace = TRUE)
t.start[[4]] = sample(1:L, size = n4, prob = pi.start[[4]], replace = TRUE)
# Initialization of global-level indicators for the different populations
k.start = replicate(n = J, list())
k.start[[1]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)
k.start[[2]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)
k.start[[3]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)
k.start[[4]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)

num_iter = 100000 # Number of MCMC iterations

# Initialize a list to store the MCMC samples of t_ji
t.samples <- replicate(n = J, list(replicate(n = num_iter, list(0))))
t.samples[[1]][[1]] = t.start[[1]] # Initialize the MCMC starting value for population 1
t.samples[[2]][[1]] = t.start[[2]] # Initialize the MCMC starting value for population 2
t.samples[[3]][[1]] = t.start[[3]] # Initialize the MCMC starting value for population 3
t.samples[[4]][[1]] = t.start[[4]] # Initialize the MCMC starting value for population 4

# Initialize a list to store the MCMC samples of k_jt
k.samples <- replicate(n = J, list(replicate(n = num_iter, list(0))))
k.samples[[1]][[1]] = k.start[[1]] # Initialize the MCMC starting value for population 1
k.samples[[2]][[1]] = k.start[[2]] # Initialize the MCMC starting value for population 2
k.samples[[3]][[1]] = k.start[[3]] # Initialize the MCMC starting value for population 3
k.samples[[4]][[1]] = k.start[[4]] # Initialize the MCMC starting value for population 4

# Initialize a list to store the MCMC samples of pi_j
pi.samples <- replicate(n = J, list(replicate(n = num_iter, list())))
pi.samples[[1]][[1]] = pi.start[[1]] # Initialize the MCMC starting value for population 1
pi.samples[[2]][[1]] = pi.start[[2]] # Initialize the MCMC starting value for population 2
pi.samples[[3]][[1]] = pi.start[[3]] # Initialize the MCMC starting value for population 3
pi.samples[[4]][[1]] = pi.start[[4]] # Initialize the MCMC starting value for population 4

# Initialize a list to store the MCMC samples of beta
beta.samples <- list()
beta.samples[[1]] <- beta.start # Initialize the MCMC starting value for beta

# Initialize List to store MCMC samples of global parameter (means and variances)
M.k.samples <- list(); Tau2.k.samples <- list() 
# Initialize List to store MCMC samples of local parameter (means and variances)
M.j.samples <- replicate(n = J, list(replicate(n = num_iter, list())))
Tau2.j.samples <- replicate(n = J, list(replicate(n = num_iter, list())))
# Initialize List to store MCMC samples of alpha
alpha.samples <- list(); alpha.samples[[1]] = alpha.start # Starting value of alpha
# Initialize List to store MCMC samples of gamma
m.samples <- list(); m.samples[[1]] = m.start  # Starting value of gamma

library(Rcpp)
sourceCpp("functions.cpp")

library(utilities)
###############################################################################
# Save global variables as a list across the populations
###############################################################################
X.global = list(X1, X2[-1, ], X3[-c(1,2), ], X4[-1, ]) 
###############################################################################
# Save local variables as a list across the populations
###############################################################################
X.local = list(X2[1, ],
               X3[c(1, 2), ],
               X4[1, ])
###############################################################################
# Save local variables as a list of matrices across the populations
###############################################################################
X.local.matrix = list(X2[1, , drop = FALSE],
                      X3[c(1, 2), ],
                      X4[1, , drop = FALSE])
###############################################################################
# WE CENTER THE LOCAL VARIABLES 
X.local = list(utilities::rm.attr(scale(X2[1, ], scale = FALSE)),
               t(apply(X3[c(1, 2), ], 1, scale, scale = FALSE)),
               utilities::rm.attr(scale(X4[1, ], scale = FALSE)))

X.local.matrix = list(t(apply(X2[1, , drop = FALSE], 1, scale, scale = FALSE)),
                      t(apply(X3[c(1, 2), ], 1, scale, scale = FALSE)),
                      t(apply(X4[1, , drop = FALSE], 1, scale, scale = FALSE)))

###############################################################################
# NOTE: POPULATION 1 HAS NO LOCAL COVARIATE. ONLY GENE EXPRESSION DATA.
# NOTE: POPULATION 2 and 4 HAVE ONE LOCAL COVARIATE
# NOTE: POPULATION 3 HAS TWO LOCAL COVARIATES
################################################################################
# BEGIN MCMC SAMPLING USING PROPOSED GLOCAL DP
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
    m[[j]] = count_my(t.samples[[j]][[iter - 1]], L)
  }
  
  
  
  for(j in 1:J){
    pi.samples[[j]][[iter]] = as.numeric(rdirichlet(n = 1, alpha = m[[j]] + alpha.samples[[iter - 1]]/L))
    pi.samples[[j]][[iter]] = (pi.samples[[j]][[iter]] + 1e-8)
    pi.samples[[j]][[iter]] = (pi.samples[[j]][[iter]])/sum(pi.samples[[j]][[iter]])
  }
  
  d = replicate(n = J, list(0))
  for(j in 1:J){
    d[[j]] = count_my(k.samples[[j]][[iter - 1]], L)}
  
  d.all = Reduce(`+`, d)
  
  beta.samples[[iter]] = as.numeric(rdirichlet(n = 1, alpha = d.all + m.samples[[iter - 1]]/L))
  
  beta.samples[[iter]] = (beta.samples[[iter]] + 1e-8)/sum(beta.samples[[iter]] + 1e-8)
  
  n.k = replicate(n = J, list(0))
  for(j in 1:J){
    n.k[[j]] = count_my(k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]], L)}
  
  n.k.all = Reduce(`+`, n.k)
  
  tau2.k.samples <- list()
  m.k.samples <- list()
  
  x.k.bar = list()
  
  for(k in 1:L){
    if(n.k.all[k] == 0){
      x.k.bar[[k]] = rep(0, nrow(X.global[[1]]))
    }else{
      x.k.bar[[k]] = (rowSums(X.global[[1]][, k.samples[[1]][[iter - 1]][t.samples[[1]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[2]][, k.samples[[2]][[iter - 1]][t.samples[[2]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[3]][, k.samples[[3]][[iter - 1]][t.samples[[3]][[iter - 1]]] == k, drop = FALSE]) + 
                        rowSums(X.global[[4]][, k.samples[[4]][[iter - 1]][t.samples[[4]][[iter - 1]]] == k, drop = FALSE]))/n.k.all[k]
    }
  }
  
  x.k.bar.matrix <- matrix(unlist(x.k.bar), nrow = nrow(X.global[[1]]), ncol = L)
  
  S.k = 0
  Z.k = 0
  alpha.hat = 0; beta.hat = 0; lambda.hat = 0; m.k.hat = list()
  
  for(k in 1:L){
    if(n.k.all[k] == 0){
      S.k[k] = 0
      Z.k[k] = 0
    }else{
      S.k.j = 0
      for(j in 1:J){
        S.k.j[j] = sum((X.global[[j]][, k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]] == k] - x.k.bar.matrix[, k])^2)
        
      }
      S.k[k] = sum(S.k.j)
      Z.k[k] = ((lambda.global * n.k.all[k])/(lambda.global + n.k.all[k])) * sum((x.k.bar.matrix[, k] - m0)^2)
    }
    alpha.hat[k] = alpha0 + 0.5 * n.k.all[k] * nrow(X.global[[1]])
    beta.hat[k] = beta0 + 0.5 * (S.k[k] + Z.k[k])
    lambda.hat[k] = lambda.global + n.k.all[k]
    m.k.hat[[k]] = ((lambda.global * m0) + (n.k.all[k] * x.k.bar.matrix[, k]))/lambda.hat[k]
  }
  
  for(k in 1:L){
    tau2.k.samples[[k]] = rinvgamma(n = 1, alpha = alpha.hat[k], beta = beta.hat[k])
    m.k.samples[[k]] = mvrnorm(n = 1, mu = m.k.hat[[k]], Sigma = diag(tau2.k.samples[[k]]/lambda.hat[k], nrow = nrow(X.global[[1]])))
  }
  
  Tau2.k.samples[[iter]] = unlist(tau2.k.samples)
  M.k.samples[[iter]] = m.k.samples
  
  n.j.t = replicate(n = J, list(0));
  
  for(t in 1:L){
    for(j in 1:J){
      n.j.t[[j]][t] = sum(t.samples[[j]][[iter - 1]] == t)
    }
  }
  
  m.j.samples <- replicate(n = 3, list())
  tau2.j.samples <- replicate(n = 3, list(0))
  x.j.bar.local = replicate(n = 3, list())
  
  for(t in 1:L){
    if(n.j.t[[2]][t] == 0){
      x.j.bar.local[[1]][[t]] = 0
    }else{
      x.j.bar.local[[1]][[t]] = sum(X.local[[1]][t.samples[[2]][[iter - 1]] == t])/n.j.t[[2]][t]
    }
  }
  
  for(t in 1:L){
    if(n.j.t[[3]][t] == 0){
      x.j.bar.local[[2]][[t]] = rep(0, nrow(X.local[[2]]))
    }else{
      x.j.bar.local[[2]][[t]] = rowSums(X.local[[2]][, t.samples[[3]][[iter - 1]] == t, drop = FALSE])/n.j.t[[3]][t]
    }
  }
  
  for(t in 1:L){
    if(n.j.t[[4]][t] == 0){
      x.j.bar.local[[3]][[t]] = 0
    }else{
      x.j.bar.local[[3]][[t]] = sum(X.local[[3]][t.samples[[4]][[iter - 1]] == t])/n.j.t[[4]][t]
    }
  }
  
  S.j.local = replicate(n = 3, list(0))
  Z.j.local = replicate(n = 3, list(0))
  alpha.j.hat = replicate(n = 3, list(0)); beta.j.hat = replicate(n = 3, list(0));
  lambda.j.hat = replicate(n = 3, list(0)); m.j.hat = replicate(n = 3, list())
  
  for(j in 1:3){
    p = nrow(X.local.matrix[[j]])
    for(t in 1:L){
      if(n.j.t[[(j + 1)]][t] == 0){
        S.j.local[[j]][t] = 0
        Z.j.local[[j]][t] = 0
        alpha.j.hat[[j]][t] = alpha.j0 + 0.5 * n.j.t[[(j+1)]][t] * p
        beta.j.hat[[j]][t] = beta.j0 + 0.5 * (S.j.local[[j]][t] + Z.j.local[[j]][t])
        lambda.j.hat[[j]][t] = lambda.local[j] + n.j.t[[(j+1)]][t]
        m.j.hat[[j]][[t]] = ((lambda.local[j] * m.j0[[j]]) + (n.j.t[[(j+1)]][t] * x.j.bar.local[[j]][[t]]))/lambda.j.hat[[j]][t]
      }else{
        S.j.local[[j]][t] = sum((X.local.matrix[[j]][, t.samples[[(j+1)]][[iter - 1]] == t] - x.j.bar.local[[j]][[t]])^2)
        Z.j.local[[j]][t] = ((lambda.local[j] * n.j.t[[(j+1)]][t])/(lambda.local[j] + n.j.t[[(j+1)]][t])) * sum((x.j.bar.local[[j]][[t]] - m.j0[[j]])^2)
        alpha.j.hat[[j]][t] = alpha.j0 + 0.5 * n.j.t[[(j+1)]][t] * p
        beta.j.hat[[j]][t] = beta.j0 + 0.5 * (S.j.local[[j]][t] + Z.j.local[[j]][t])
        lambda.j.hat[[j]][t] = lambda.local[j] + n.j.t[[(j+1)]][t]
        m.j.hat[[j]][[t]] = ((lambda.local[j] * m.j0[[j]]) + (n.j.t[[(j+1)]][t] * x.j.bar.local[[j]][[t]]))/lambda.j.hat[[j]][t]
      }
    }
  }
  
  for(j in 1:3){
    if(j == 2){
      p = nrow(X.local.matrix[[j]])
      for(t in 1:L){
        tau2.j.samples[[j]][t] = rinvgamma(n = 1, alpha = alpha.j.hat[[j]][t], beta = beta.j.hat[[j]][t])
        m.j.samples[[j]][[t]] = mvrnorm(n = 1, mu = m.j.hat[[j]][[t]], Sigma = diag(tau2.j.samples[[j]][t] /lambda.j.hat[[j]][t], nrow = p))
      }
    }else{
      for(t in 1:L){
        tau2.j.samples[[j]][t] = rinvgamma(n = 1, alpha = alpha.j.hat[[j]][t], beta = beta.j.hat[[j]][t])
        m.j.samples[[j]][[t]] = rnorm(n = 1, mean = m.j.hat[[j]][[t]], sd = sqrt(tau2.j.samples[[j]][t] /lambda.j.hat[[j]][t]))
      }
      
    }
    
    M.j.samples[[j]][[iter]] = m.j.samples[[j]]
    Tau2.j.samples[[j]][[iter]] = tau2.j.samples[[j]]
  }
  
  
  # Calculate the probabilities of table indices and sample
  n = c(n1, n2, n3, n4)
  
  # Define an array of global cluster variance samples 
  Tau2.k.samples.array = array(0, dim = c(p_G, p_G, L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:length(Tau2.k.samples[[iter]])){
    Tau2.k.samples.array[,,j] = diag(Tau2.k.samples[[iter]][j], nrow = p_G)
  }
  # Define an array of local cluster variance samples (for population 2) which is 2-variate Normal
  Tau2.j.samples.array2 = array(0, dim = c(nrow(X.local[[2]]), nrow(X.local[[2]]), L))
  for(j in 1:length(Tau2.j.samples[[2]][[iter]])){
    Tau2.j.samples.array2[,,j] = diag(Tau2.j.samples[[2]][[iter]][j], nrow = nrow(X.local[[2]]))
  }
  
  t.prob = replicate(n = J, list())
  for(j in 1:J){
    t.prob[[j]] = matrix(0, nrow = n[j], ncol = L)
  }
  
  t.prob[[1]] <- prob_exponent_no_local2(pi_samples = pi.samples[[1]][[iter]],  
                                         phi_samples = matrix(unlist(m.k.samples),
                                                              nrow = nrow(X.global[[1]]), L),
                                         X_global = t(X.global[[1]]), 
                                         k = k.samples[[1]][[iter - 1]], 
                                         Sigma_G =  Tau2.k.samples.array)
  
  t.prob[[2]] <- prob_exponent_univ_local(pi_samples = pi.samples[[2]][[iter]], 
                                          X_local = as.numeric(X.local[[1]]), 
                                          psi_samples = as.numeric(m.j.samples[[1]]),
                                          Sigma = as.numeric(sqrt(tau2.j.samples[[1]])), 
                                          phi_samples = matrix(unlist(m.k.samples),
                                                               nrow = nrow(X.global[[1]]), L),
                                          X_global = t(X.global[[2]]), 
                                          k = k.samples[[2]][[iter - 1]], 
                                          Sigma_G =  Tau2.k.samples.array)
  
  t.prob[[3]] = prob_exponent_mv_local(pi_samples = pi.samples[[3]][[iter]], 
                                       X_local = t(X.local[[2]]), 
                                       psi_samples = matrix(unlist(M.j.samples[[2]][[iter]]), nrow = nrow(X.local[[2]]), ncol = L),
                                       Sigma = Tau2.j.samples.array2,
                                       phi_samples = matrix(unlist(m.k.samples), nrow = nrow(X.global[[2]]), L),
                                       X_global = t(X.global[[3]]),
                                       k = k.samples[[3]][[iter - 1]], 
                                       Sigma_G = Tau2.k.samples.array)
  
  t.prob[[4]] <- prob_exponent_univ_local(pi_samples = pi.samples[[4]][[iter]], 
                                          X_local = as.numeric(X.local[[3]]), 
                                          psi_samples = as.numeric(m.j.samples[[3]]),
                                          Sigma = as.numeric(sqrt(tau2.j.samples[[3]])), 
                                          phi_samples = matrix(unlist(m.k.samples),
                                                               nrow = nrow(X.global[[4]]), L),
                                          X_global = t(X.global[[4]]), 
                                          k = k.samples[[4]][[iter - 1]], 
                                          Sigma_G =  Tau2.k.samples.array)
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){
    t.prob[[j]] = calc_probability_log_sum_exp_normalized(t.prob[[j]])
  }
  
  for(j in 1:J){
    t.samples[[j]][[iter]] = sample_my(t.prob[[j]])
  }
  
  
  prob = array(0, dim = c(L, L, J))
  
  for(j in 1:J){
    prob[, , j] = prob_exponent_mv_global(beta_samples = beta.samples[[iter]], 
                                          X_global = t(X.global[[j]]), 
                                          phi_samples = matrix(unlist(m.k.samples), nrow = nrow(X.global[[j]]), L), 
                                          Sigma_G = Tau2.k.samples.array,
                                          t_samples = t.samples[[j]][[iter]])}
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){prob[, ,j] = calc_probability_log_sum_exp_normalized(prob[, ,j])}
  
  for(j in 1:J){
    k.samples[[j]][[iter]] = sample_my(prob[, , j])}
  
  sampled.alpha = sample_alpha(ycurrent = alpha.samples[[iter - 1]], 
                               pars = list(a0 = 0.1, b0 = 0.1, 
                                           pi1 = pi.samples[[1]][[iter]],
                                           pi2 = pi.samples[[2]][[iter]], 
                                           pi3 = pi.samples[[3]][[iter]],
                                           pi4 = pi.samples[[4]][[iter]],
                                           J = 4))
  alpha.samples[[iter]] = sampled.alpha$out
  
  sampled.m = sample_m(ycurrent = m.samples[[iter - 1]], pars = list(a1 = 0.1, b1 = 0.1, beta = beta.samples[[iter]]))
  m.samples[[iter]] = sampled.m$out
  
}
time.end = Sys.time()
time.end - time.start

################################################################################
# CALCULATE LOG-POSTERIOR TO PLOT TRACEPLOT AND ACF OF LOG-POSTERIOR
################################################################################
ll = 0
library(mvtnorm)
for(iter in 2:num_iter){
  # Printing the iterations
  if(iter == 2){
    cat(paste0("Iteration: ", iter-1, "\n"))
  }
  if(iter %% floor((10/100)*(num_iter + 1)) == 0) {
    cat(paste0("Iteration: ", iter, "\n"))
  }
  
  # Define an array of global cluster variance samples
  Tau2.k.samples.array = array(0, dim = c(p_G, p_G, L))
  # This is needed as the Rcpp function to calculate LL needs a cube for variance samples
  for(j in 1:length(Tau2.k.samples[[iter]])){
    Tau2.k.samples.array[,,j] = diag(Tau2.k.samples[[iter]][j], nrow = p_G)
  }
  # Define an array of local cluster variance samples (for population 2) which is 2-variate Normal
  Tau2.j.samples.array2 = array(0, dim = c(nrow(X.local[[2]]), nrow(X.local[[2]]), L))
  for(j in 1:length(Tau2.j.samples[[2]][[iter]])){
    Tau2.j.samples.array2[,,j] = diag(Tau2.j.samples[[2]][[iter]][j], nrow = nrow(X.local[[2]]))
  }
  
  # To calculate the LL Psi-samples are replaced by the M.j.samples (Mean of Local NIG posterior)
  # Phi-samples are replaced by M.k.samples (Mean of Global NIG posterior)
  logsum = logll_no_local2(pi_samples = pi.samples[[1]][[iter]],
                           t_samples = t.samples[[1]][[iter]],
                           phi_samples = matrix(unlist(M.k.samples[[iter]]), nrow = p_G, ncol = L),
                           X_global = X.global[[1]],
                           k_samples = k.samples[[1]][[iter]],
                           Sigma_G = Tau2.k.samples.array) +
    logll_univ_local2(pi.samples[[2]][[iter]], X.local[[1]], unlist(M.j.samples[[1]][[iter]]),
                      Sigma = sqrt(Tau2.j.samples[[1]][[iter]]),
                      t.samples[[2]][[iter]], matrix(unlist(M.k.samples[[iter]]), nrow = p_G, ncol = L),
                      X.global[[2]], k.samples[[2]][[iter]],
                      Sigma_G = Tau2.k.samples.array) +
    logll_mv_local2(pi.samples[[3]][[iter]], X.local[[2]],
                    matrix(unlist(M.j.samples[[2]][[iter]]), nrow = nrow(X.local[[2]]), ncol = L),
                    Sigma = Tau2.j.samples.array2,
                    t.samples[[3]][[iter]],
                    matrix(unlist(M.k.samples[[iter]]), nrow = p_G, ncol = L),
                    X.global[[3]], k.samples[[3]][[iter]],
                    Sigma_G = Tau2.k.samples.array) +
    logll_univ_local2(pi.samples[[4]][[iter]], X.local[[3]], unlist(M.j.samples[[3]][[iter]]),
                      Sigma = sqrt(Tau2.j.samples[[3]][[iter]]),
                      t.samples[[4]][[iter]], matrix(unlist(M.k.samples[[iter]]), nrow = p_G, ncol = L),
                      X.global[[4]], k.samples[[4]][[iter]],
                      Sigma_G = Tau2.k.samples.array)
  
  m1.tilde = count_my(t.samples[[1]][[iter]], L)
  m2.tilde = count_my(t.samples[[2]][[iter]], L)
  m3.tilde = count_my(t.samples[[3]][[iter]], L)
  m4.tilde = count_my(t.samples[[4]][[iter]], L)
  
  logsum = logsum + sum(m1.tilde * log(pi.samples[[1]][[iter]])) + sum(m2.tilde * log(pi.samples[[2]][[iter]])) + sum(m3.tilde * log(pi.samples[[3]][[iter]])) + sum(m4.tilde * log(pi.samples[[4]][[iter]]))
  
  m.tilde = count_my(k.samples[[1]][[iter]], L) + count_my(k.samples[[2]][[iter]], L) + count_my(k.samples[[3]][[iter]], L) + count_my(k.samples[[4]][[iter]], L)
  
  logsum = logsum + sum(m.tilde * log(beta.samples[[iter]]))
  
  Local_density2 = 0
  for(l in 1:length(M.j.samples[[1]][[iter]])){
    Local_density2[l] <- dnorm(M.j.samples[[1]][[iter]][[l]], mean = m.j0[[1]], sd = sqrt(1/(lambda.local[1]) * Tau2.j.samples[[1]][[iter]])[l], log = TRUE)
  }
  Local_density2_sum = sum(Local_density2) + sum(dinvgamma(Tau2.j.samples[[1]][[iter]], alpha = alpha0, beta = beta0, log = TRUE))
  
  Local_density3 = 0
  for(l in 1:length(M.j.samples[[2]][[iter]])){
    Local_density3[l] <- dmvnorm(M.j.samples[[2]][[iter]][[l]], mean = m.j0[[2]], sigma = diag(1/(lambda.local[2]) * Tau2.j.samples[[2]][[iter]][l], length(m.j0[[2]])), log = TRUE)
  }
  
  Local_density3_sum = sum(Local_density3) + sum(dinvgamma(Tau2.j.samples[[2]][[iter]], alpha = alpha0, beta = beta0, log = TRUE))
  
  Local_density4 = 0
  for(l in 1:length(M.j.samples[[3]][[iter]])){
    Local_density4[l] <- dnorm(M.j.samples[[3]][[iter]][[l]], mean = m.j0[[3]], sd = sqrt(1/(lambda.local[3]) * Tau2.j.samples[[3]][[iter]])[l], log = TRUE)
  }
  
  Local_density4_sum = sum(Local_density4) + sum(dinvgamma(Tau2.j.samples[[3]][[iter]], alpha = alpha0, beta = beta0, log = TRUE))
  
  ## Component 1
  Local_density_sum = Local_density2_sum + Local_density3_sum + Local_density4_sum
  
  Global_density = 0
  for(l in 1:length(M.k.samples[[iter]])){
    Global_density[l] <- dmvnorm(M.k.samples[[iter]][[l]], mean = m0, sigma = diag(1/(lambda.global) * Tau2.k.samples[[iter]][l], length(m0)), log = TRUE)
  }
  ## Component 2
  Global_density_sum = sum(Global_density) + sum(dinvgamma(Tau2.k.samples[[iter]], alpha = alpha0, beta = beta0, log = TRUE))
  ## Component 3
  pi_density_sum = ddirichlet(pi.samples[[1]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE) + ddirichlet(pi.samples[[2]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE) + ddirichlet(pi.samples[[3]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE) + ddirichlet(pi.samples[[4]][[iter]], alpha = rep(alpha.samples[[iter]]/L, L), log = TRUE)
  ## Component 4
  beta_density_sum = ddirichlet(beta.samples[[iter]], alpha = rep(m.samples[[iter]]/L, L), log = TRUE)
  
  ll[iter] = logsum + Local_density_sum + Global_density_sum + pi_density_sum + beta_density_sum + dgamma(alpha.samples[[iter]], shape = 0.1, rate = 0.1, log = TRUE) + dgamma(m.samples[[iter]], shape = 0.1, rate = 0.1, log = TRUE) 
}

################################################################################
# PLOT TRACEPLOT OF LOG-LIKELIHOOD AND ACF
################################################################################
burn = 25000 # Burn-in 
samples.thin = seq((burn + 1), (num_iter), by = 75)

log_like <- ll[samples.thin]

library(tidyverse)
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
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

library(forecast)
ACF_plot <-  ggAcf(x = log_like, lag.max = 40) +
  ggtitle("") +
  labs(y = "") +
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

gridExtra::grid.arrange(ll_plot, ACF_plot, ncol = 2)

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
x.limit.lower <- min(X1[1, ], X2[2, ], X3[3, ], X4[2, ])
x.limit.upper <- max(X1[1, ], X2[2, ], X3[3, ], X4[2, ])

y.limit.lower <- min(X1[2, ], X2[3, ], X3[4, ], X4[3, ])
y.limit.upper <- max(X1[2, ], X2[3, ], X3[4, ], X4[3, ])

###########################################################
# GLOBAL LEVEL CLUSTERING OF GLOBAL VARIABLES
###########################################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = c(k.samples[[1]][[iter]][t.samples[[1]][[iter]]],
                    k.samples[[2]][[iter]][t.samples[[2]][[iter]]],
                    k.samples[[3]][[iter]][t.samples[[3]][[iter]]],
                    k.samples[[4]][[iter]][t.samples[[4]][[iter]]])
}
samples = samples.thin
k.rand = getDahl(index[samples],
                 NULL 
)

best.index = samples[k.rand$DahlIndex]
posterior_samples <- matrix(0, nrow = length(samples), ncol = length(index[samples][[1]]))

for(i in 1:length(samples)){
  posterior_samples[i, ] = index[samples][[i]]
}

library(mcclust.ext)
sim_mat = comp.psm(posterior_samples)
par(mfrow = c(1, 1))
plotpsm(sim_mat)

reorder_dismat <-  function(dismat, groups, order.groups=NULL){
  # Use correlation between variables as distance
  order.dis   = integer(0)
  J           = length(unique(groups))
  if(is.null(order.groups)){
    order.j   = 1:J
  } else {
    order.j   = order.groups
  }
  for (j in order.j){
    groups.j  = which(groups==j)
    dd        = as.dist((1-dismat[groups.j, groups.j])/2)
    hc        = hclust(dd)
    order.dis = c(order.dis, hc$order+length(order.dis))
  }
  dismat      = dismat[order.dis, order.dis]
  dismat      = dismat[nrow(dismat):1,]
}

dissimlar_stable = sim_mat; N_S_map = ncol(sim_mat)
library(scales)
library(tidyverse)
library(reshape2)
dismat      = round(dissimlar_stable, 2)
dismat      = reorder_dismat(dismat,groups=rep(1, N_S_map))
plot_dismat = reshape2::melt(dismat)

ggplot(data=plot_dismat, aes(x=factor(Var1), y=factor(Var2), fill=value))+ 
  geom_tile() + theme_bw()+ 
  scale_y_discrete(limits = seq(N_S_map, 1),
                   breaks = floor(seq(1,N_S_map,length.out = 9)), 
                   labels = floor(seq(1,N_S_map,length.out = 9))) +
  scale_x_discrete(breaks = floor(seq(1,N_S_map,length.out = 9)), 
                   labels = floor(seq(1,N_S_map,length.out = 9))) +
  xlab("observation")+ylab("observation")+
  scale_fill_gradientn(colours = c("white", "yellow", "red"), 
                       values = rescale(c(0,0.25,1)), space = "Lab", name="")+
  theme(legend.position = "right", text = element_text(size=20))


###############################################################
# Some Pre-processing before plotting the estimated clustering
###############################################################
singleton = as.numeric(names(which(table(index[[best.index]]) <= 0.0 * length(index[[best.index]]))))

singleton.index = which(index[[best.index]] %in% singleton)

z.estimated = index[[best.index]]

cluster.global <- data.frame(x = c(X.global[[1]][1, ],
                                   X.global[[2]][1, ],
                                   X.global[[3]][1, ],
                                   X.global[[4]][1, ]),
                             y = c(X.global[[1]][2, ],
                                   X.global[[2]][2, ],
                                   X.global[[3]][2, ],
                                   X.global[[4]][2, ]),
                             cluster = factor(z.estimated),
                             cancer = factor(c(rep("Stomach Cancer", length(X.global[[1]][1,])),
                                               rep("Rectal Cancer", length(X.global[[2]][1,])),
                                               rep("Colon Cancer", length(X.global[[3]][1,])),
                                               rep("Esophageal Cancer", length(X.global[[4]][1,])))))

if(length(singleton.index) > 0){
  z.estimated = z.estimated[-singleton.index_mmclust2]
  cluster.global = cluster.global[-singleton.index, ]
}else{
  z.estimated = z.estimated
  cluster.global = cluster.global
}


myvalues =  c("1" = "#F8766D",
              "2" = "#00BA38",
              "3" = "#619CFF", 
              "4" = "blueviolet",
              "5" = "cyan4",
              "6" = "#E6AB02",
              "7" = "#E36EF6",
              "8" = "bisque4",
              "9" = "coral4",
              "10" = "darkslateblue",
              
              "11" = "lightseagreen",
              "12" = "#E69F00", 
              "13" = "#AA3377",
              "14" = "sienna3",
              "15" = "hotpink",
              "16" = "sienna4",
              "17" = "hotpink3",
              "18" = "sienna1",
              "19" = "dodgerblue4",
              "20" = "bisque2",
              
              "21" = "darkgreen",
              "22" = "orange", 
              "23" = "maroon2",
              "24" = "sienna2",
              "25" = "hotpink4",
              "26" = "sienna3",
              "27" = "brown4",
              "28" = "sienna1",
              "29" = "dodgerblue3",
              "30" = "bisque3")

library(patchwork)
library(tidyverse)
library(latex2exp)

plot.global.LS <- cluster.global %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3) + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + facet_wrap(~cancer, ncol = 2) + 
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

plot.global.LS

#######################################################
# LOCAL LEVEL CLUSTERING
#######################################################
##########################################
# RECTAL CANCER
##########################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = c(t.samples[[2]][[iter]])
}

samples <- samples.thin
posterior_samples <- matrix(0, nrow = length(samples), ncol = length(index[samples][[1]]))
for(i in 1:length(samples)){
  posterior_samples[i, ] = index[samples][[i]]
}
library(mcclust.ext)
sim_mat = comp.psm(posterior_samples)
par(mfrow = c(1, 1))
plotpsm(sim_mat)

dissimlar_stable = sim_mat; N_S_map = ncol(sim_mat)
library(scales)
library(tidyverse)
dismat      = round(dissimlar_stable, 2)
dismat      = reorder_dismat(dismat,groups=rep(1, N_S_map))
plot_dismat = reshape2::melt(dismat)

ggplot(data=plot_dismat, aes(x=factor(Var1), y=factor(Var2), fill=value))+ 
  geom_tile() + theme_bw()+ 
  scale_y_discrete(limits = seq(N_S_map, 1),
                   breaks = floor(seq(1,N_S_map,length.out = 9)), 
                   labels = floor(seq(1,N_S_map,length.out = 9))) +
  scale_x_discrete(breaks = floor(seq(1,N_S_map,length.out = 9)), 
                   labels = floor(seq(1,N_S_map,length.out = 9))) +
  xlab("observation")+ylab("observation")+
  scale_fill_gradientn(colours = c("white", "yellow", "red"), 
                       values = rescale(c(0,0.25,1)), space = "Lab", name="")+
  theme(legend.position = "right", text = element_text(size=20))

k2.rand = getDahl(index[samples],
                  NULL)
best.index = samples[k2.rand$DahlIndex]
singleton_2.1 = as.numeric(names(which(table(index[[best.index]]) <= 0.02 * length(index[[best.index]]))))
singleton.index_2.1 = which(index[[best.index]] %in% singleton_2.1)
z.estimated_2.1 = index[[best.index]]

cluster.local_2.1 <- data.frame(x = c(X.global[[2]][1, ]),
                                y = c(X.global[[2]][2, ]), 
                                cluster = factor(c(z.estimated_2.1)),
                                cancer = factor(c(rep("Rectal Cancer", length(X.global[[2]][1,])))))

if(length(singleton.index_2.1) > 0){
  z.estimated_2.1 = z.estimated_2.1[-singleton.index_2.1]
  cluster.local_2.1 = cluster.local_2.1[-singleton.index_2.1, ]
}else{
  z.estimated_2.1 = z.estimated_2.1
  cluster.local_2.1 = cluster.local_2.1
}
##########################################
# COLON CANCER
##########################################
index <- list()
for(iter in 2:num_iter){
  index[[iter]] = c(t.samples[[3]][[iter]])
}

samples <- samples.thin
posterior_samples <- matrix(0, nrow = length(samples), ncol = length(index[samples][[1]]))
for(i in 1:length(samples)){
  posterior_samples[i, ] = index[samples][[i]]
}
library(mcclust.ext)
sim_mat = comp.psm(posterior_samples)
par(mfrow = c(1, 1))
plotpsm(sim_mat)

dissimlar_stable = sim_mat; N_S_map = ncol(sim_mat)
library(scales)
library(tidyverse)
dismat      = round(dissimlar_stable, 2)
dismat      = reorder_dismat(dismat,groups=rep(1, N_S_map))
plot_dismat = reshape2::melt(dismat)

ggplot(data=plot_dismat, aes(x=factor(Var1), y=factor(Var2), fill=value))+ 
  geom_tile() + theme_bw()+ 
  scale_y_discrete(limits = seq(N_S_map, 1),
                   breaks = floor(seq(1,N_S_map,length.out = 9)), 
                   labels = floor(seq(1,N_S_map,length.out = 9))) +
  scale_x_discrete(breaks = floor(seq(1,N_S_map,length.out = 9)), 
                   labels = floor(seq(1,N_S_map,length.out = 9))) +
  xlab("observation")+ylab("observation")+
  scale_fill_gradientn(colours = c("white", "yellow", "red"), 
                       values = rescale(c(0,0.25,1)), space = "Lab", name="")+
  theme(legend.position = "right", text = element_text(size=20))

k3.rand = getDahl(index[samples],
                  NULL)
best.index = samples[k3.rand$DahlIndex]
singleton_2.2 = as.numeric(names(which(table(index[[best.index]]) <= 0.02 * length(index[[best.index]]))))
singleton.index_2.2 = which(index[[best.index]] %in% singleton_2.2)
z.estimated_2.2 = index[[best.index]]


cluster.local_2.2 <- data.frame(x = c(X.global[[3]][1, ]),
                                y = c(X.global[[3]][2, ]), 
                                cluster = factor(c(z.estimated_2.2)),
                                cancer = factor(c(rep("Colon Cancer", length(X.global[[3]][1,])))))

if(length(singleton.index_2.2) > 0){
  z.estimated_2.2 = z.estimated_2.2[-singleton.index_2.2]
  cluster.local_2.2 = cluster.local_2.2[-singleton.index_2.2, ]
}else{
  z.estimated_2.2 = z.estimated_2.2
  cluster.local_2.2 = cluster.local_2.2
}

cluster.local = bind_rows(cluster.local_2.1, cluster.local_2.2)

cluster.local %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3) + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + facet_wrap(~cancer, ncol = 2) + scale_color_manual(values = myvalues) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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
# RELABELLING LOCAL CLUSTERS FOR RECTAL CANCER
#################################################
cluster.global.LS.rectal = cluster.global %>% filter(cancer == "Rectal Cancer") 
cluster.local.LS.rectal = cluster.local %>% filter(cancer == "Rectal Cancer") 
cluster.local.LS.rectal$cluster_org = cluster.local.LS.rectal$cluster

charaters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")

for(l in 1:length(unique(cluster.global.LS.rectal$cluster))){
  if(length(singleton.index_2.1) > 0){
    cluster.global.LS.rectal.subset = cluster.global.LS.rectal$cluster[-singleton.index_2.1]
  }else{
    cluster.global.LS.rectal.subset = cluster.global.LS.rectal$cluster
  }
  
  cluster.index = as.numeric(as.character(unique(cluster.global.LS.rectal.subset)[l]))
  
  data.index = which(cluster.global.LS.rectal.subset == cluster.index)
  
  labels.local = cluster.local.LS.rectal$cluster[which(cluster.global.LS.rectal.subset == cluster.index)]
  labels.local.unique =  as.character(unique(cluster.local.LS.rectal$cluster[which(cluster.global.LS.rectal.subset == cluster.index)]))
  
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
  cluster.local.LS.rectal$cluster <- as.character(cluster.local.LS.rectal$cluster)
  cluster.local.LS.rectal$cluster[data.index] <- as.character(labels.local)
}

cluster.local.LS.rectal$cluster <- factor(cluster.local.LS.rectal$cluster)

cluster.local.LS.rectal$cluster

myvalues_local.LS = c("1" = "#F8766D",
                      "2" = "#00BA38",
                      "3" = "#619CFF", 
                      "4" = "blueviolet",
                      "5" = "cyan4",
                      "6" = "#E6AB02",
                      "7" = "#E36EF6",
                      "8" = "bisque4",
                      "9" = "coral4",
                      "7a" = "yellow3",
                      "7b" = "sienna1",
                      
                      "11" = "lightseagreen",
                      "12" = "#E69F00", 
                      "13" = "#AA3377",
                      "14" = "sienna3",
                      "15" = "hotpink",
                      "16" = "sienna4",
                      "17" = "hotpink3",
                      "18" = "sienna1",
                      "19" = "dodgerblue4",
                      "20" = "bisque2",
                      
                      "21" = "darkgreen",
                      "22" = "orange", 
                      "23" = "maroon2",
                      "24" = "sienna2",
                      "25" = "hotpink4",
                      "26" = "sienna3",
                      "27" = "brown4",
                      "28" = "sienna1",
                      "29" = "dodgerblue3",
                      "30" = "bisque3")

cluster.local.LS.rectal %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3) + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + scale_color_manual(values = myvalues_local.LS) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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

library(aricode)
cluster.global.LS.colon = cluster.global %>% filter(cancer == "Colon Cancer") 
cluster.local.LS.colon = cluster.local %>% filter(cancer == "Colon Cancer") 
cluster.local.LS.colon$cluster_org = cluster.local.LS.colon$cluster

for(l in 1:length(unique(cluster.global.LS.colon$cluster))){
  if(length(singleton.index_2.2) > 0){
    cluster.global.LS.colon.subset = cluster.global.LS.colon$cluster[-singleton.index_2.2]
  }else{
    cluster.global.LS.colon.subset = cluster.global.LS.colon$cluster
  }
  
  cluster.index = as.numeric(as.character(unique(cluster.global.LS.colon.subset)[l]))
  
  data.index = which(cluster.global.LS.colon.subset == cluster.index)
  
  labels.local = cluster.local.LS.colon$cluster[which(cluster.global.LS.colon.subset == cluster.index)]
  labels.local.unique =  as.character(unique(cluster.local.LS.colon$cluster[which(cluster.global.LS.colon.subset == cluster.index)]))
  
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
  cluster.local.LS.colon$cluster <- as.character(cluster.local.LS.colon$cluster)
  cluster.local.LS.colon$cluster[data.index] <- as.character(labels.local)
}

cluster.local.LS.colon$cluster <- factor(cluster.local.LS.colon$cluster)
cluster.local.LS.colon$cluster 

myvalues_local.LS = c("1" = "#F8766D",
                      "2" = "#00BA38",
                      "3" = "#619CFF", 
                      "4" = "blueviolet",
                      "5" = "cyan4",
                      "6" = "#E6AB02",
                      "7" = "#E36EF6",
                      "8" = "bisque4",
                      "9" = "coral4",
                      "7a" = "yellow3",
                      "7b" = "sienna1",
                      "7c" = "#AA3377",
                      
                      "11" = "lightseagreen",
                      "12" = "#E69F00", 
                      "13" = "#AA3377",
                      "14" = "sienna3",
                      "15" = "hotpink",
                      "16" = "sienna4",
                      "17" = "hotpink3",
                      "18" = "sienna1",
                      "19" = "dodgerblue4",
                      "20" = "bisque2",
                      
                      "21" = "darkgreen",
                      "22" = "orange", 
                      "23" = "maroon2",
                      "24" = "sienna2",
                      "25" = "hotpink4",
                      "26" = "sienna3",
                      "27" = "brown4",
                      "28" = "sienna1",
                      "29" = "dodgerblue3",
                      "30" = "bisque3"
)

cluster.local.LS.colon %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.7) + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + scale_color_manual(values = myvalues_local.LS) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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

cluster.local.LS.combined = bind_rows(cluster.local.LS.colon,
                                      cluster.local.LS.rectal)

cluster.local.LS.combined$cluster
N_cluster_local = length(unique(cluster.local.LS.combined $cluster))

myvalues_local.LS.combined = c("1" = "#F8766D",
                               "2" = "#00BA38",
                               "3" = "#619CFF", 
                               "4" = "blueviolet",
                               "5" = "cyan4",
                               "6" = "#E6AB02",
                               "7" = "#E36EF6",
                               "8" = "bisque4",
                               "9" = "coral4",
                               "7a" = "yellow3",
                               "7b" = "sienna1",
                               "7c" = "#AA3377",
                               
                               "11" = "lightseagreen",
                               "12" = "#E69F00", 
                               "13" = "#AA3377",
                               "14" = "sienna3",
                               "15" = "hotpink",
                               "16" = "sienna4",
                               "17" = "hotpink3",
                               "18" = "sienna1",
                               "19" = "dodgerblue4",
                               "20" = "bisque2",
                               
                               "21" = "darkgreen",
                               "22" = "orange", 
                               "23" = "maroon2",
                               "24" = "sienna2",
                               "25" = "hotpink4",
                               "26" = "sienna3",
                               "27" = "brown4",
                               "28" = "sienna1",
                               "29" = "dodgerblue3",
                               "30" = "bisque3"
)

plot.local.LS <- cluster.local.LS.combined %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3, alpha = 0.9) + labs(x = TeX("$UMAP_1$"), y = TeX("$UMAP_2$")) + facet_wrap(~cancer, ncol = 2) + scale_color_manual(values = myvalues_local.LS.combined) +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
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

plot.global.LS # GLOBAL-LEVEL CLUSTERING PLOT OF GLOBAL VARIABLES
plot.local.LS  # LOCAL-LEVEL CLUSTERING PLOT OF GLOBAL VARIABLES

################################################################################
## LOCAL VARIABLE PLOTS
################################################################################
cluster <- data.frame(y = as.numeric(rectal_full[, 1][-singleton.index_2.1]),
                      cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Rectal Cancer") %>% pull(cluster)))
cluster %>% group_by(cluster) %>% summarise(median(y))

localVarplot1 <- cluster %>% group_by(cluster) %>% ggplot(aes(x = y, col = cluster)) + geom_density(size = 1) + labs(x = "CEA", y = "Density", title = "Rectal Cancer")  + scale_color_manual(values = myvalues_local.LS.combined)  + theme_classic() +  
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
  ) + xlim(0, 100)



cluster <- data.frame(x = colon_full[-singleton.index_2.2, c(1)],
                      y = colon_full[-singleton.index_2.2, c(2)],
                      cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Colon Cancer") %>% pull(cluster)))

## CEA for Colon Cancer Patients
data.frame(y = as.numeric(colon_full[, 1][-singleton.index_2.2]),
           cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Colon Cancer") %>% pull(cluster))) %>% group_by(cluster) %>% summarise(mean(y))

## BMI for Colon Cancer Patients
data.frame(y = as.numeric(colon_full[, 2][-singleton.index_2.2]),
           cluster = factor(cluster.local.LS.combined %>% filter(cancer == "Colon Cancer") %>% pull(cluster))) %>% group_by(cluster) %>% summarise(mean(y))

localVarplot2 <- cluster  %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3) + labs(x = "CEA", y = "BMI", title = "Colon Cancer") + scale_color_manual(values = myvalues_local.LS.combined)  + theme_classic() + 
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
  ) + ylim(15, 70) + xlim(0, 300)

library(gridExtra)
grid.arrange(localVarplot1, localVarplot2, nrow = 2)

################################################################################
## DE ANALYSIS GLOBAL
################################################################################
library(tidyverse)
colon_full_with_id <- colon_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                                     bmi.exposures,
                                                     X1,
                                                     X2,
                                                     submitter_id.samples) %>% na.omit()


colon_full_with_id$cluster <- cluster.global %>% filter(cancer == "Colon Cancer") %>% dplyr::select(cluster) %>% pull()

colon_cluster_subset <- as.data.frame(colon)[ ,which(colnames(colon) %in% colon_full_with_id$submitter_id.samples)]
data = as.matrix(colon_cluster_subset[c(1:60483), ])
rownames(data) <- common.genes
library(scran)
clusters <- cluster.global %>% filter(cancer == "Colon Cancer") %>% dplyr::select(cluster) %>% pull()

N_DE=6
cluster.markers = scran::findMarkers(x=data, pval.type = "any", 
                                     groups=as.numeric(as.character(clusters)), "binom")

clust_VI_stable = as.numeric(as.character(clusters))
markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = cbind(t(data[markers,]), clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(mean))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Cluster", sort(unique(clust_VI_stable))), "FDR")
library(kableExtra)
Table_3[-1, ]

Boxplot_DE_colon <- function(markers_cell = markers_cell){
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  BoxPlot_Data = rbind(cbind(ord_core_stable_mat[,c(1,7)], markers[1]),
                       cbind(ord_core_stable_mat[,c(2,7)], markers[2]),
                       cbind(ord_core_stable_mat[,c(3,7)], markers[3]),
                       cbind(ord_core_stable_mat[,c(4,7)], markers[4]),
                       cbind(ord_core_stable_mat[,c(5,7)], markers[5]),
                       cbind(ord_core_stable_mat[,c(6,7)], markers[6]))
  
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  BoxPlot_Data           = data.frame(BoxPlot_Data)
  BoxPlot_Data$value     = as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster    = factor(BoxPlot_Data$cluster)
  
  ggplot(BoxPlot_Data, aes(y=value, x=gene, col=cluster)) + geom_boxplot() +
    # labs(title = "Colon Cancer") +
    theme_bw() + 
    scale_color_manual(values = myvalues) + theme(axis.title=element_blank(), 
                                                  axis.text =element_text(size=12, angle = 45, vjust = 0.5, 
                                                                          hjust=1), strip.text = element_text(size=30),
                                                  plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ))
}

# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_Boxplot_colon = Boxplot_DE_colon(markers_cell = markers_cell) 
Plot_Boxplot_colon

rectal_full_with_id <- rectal_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                                       X1,
                                                       X2,
                                                       submitter_id.samples) %>% na.omit()

rectal_full_with_id$cluster <- cluster.global_mcclust2 %>% filter(cancer == "Rectal Cancer") %>% dplyr::select(cluster) %>% pull()

rectal_cluster_subset <- as.data.frame(rectal)[, which(colnames(rectal) %in% rectal_full_with_id$submitter_id.samples)]

data = as.matrix(rectal_cluster_subset[c(1:60483), ])
rownames(data) <- common.genes
library(scran)
clusters <- cluster.global %>% filter(cancer == "Rectal Cancer") %>% dplyr::select(cluster) %>% pull()

N_DE=6
cluster.markers = scran::findMarkers(x=data, pval.type = "any", 
                                     groups=as.numeric(as.character(clusters)), "binom")

clust_VI_stable = as.numeric(as.character(clusters))
markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = cbind(t(data[markers,]), clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(mean))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Cluster", sort(unique(clust_VI_stable))), "FDR")
Table_3[-1,]

Boxplot_DE_rectal <- function(markers_cell = markers_cell){
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  BoxPlot_Data = rbind(cbind(ord_core_stable_mat[,c(1,7)], markers[1]),
                       cbind(ord_core_stable_mat[,c(2,7)], markers[2]),
                       cbind(ord_core_stable_mat[,c(3,7)], markers[3]),
                       cbind(ord_core_stable_mat[,c(4,7)], markers[4]),
                       cbind(ord_core_stable_mat[,c(5,7)], markers[5]),
                       cbind(ord_core_stable_mat[,c(6,7)], markers[6]))
  
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  BoxPlot_Data           = data.frame(BoxPlot_Data)
  BoxPlot_Data$value     = as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster    = factor(BoxPlot_Data$cluster)
  
  ggplot(BoxPlot_Data, aes(y=value, x=gene, col=cluster)) + geom_boxplot() + 
    # labs(title = "Rectal Cancer") +
    theme_bw() + 
    scale_color_manual(values = myvalues) + theme(axis.title=element_blank(), 
                                                  axis.text =element_text(size=12, angle = 45, vjust = 0.5, 
                                                                          hjust=1), strip.text = element_text(size=30),
                                                  plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ))
}
# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_Boxplot_rectal = Boxplot_DE_rectal(markers_cell = markers_cell) 
Plot_Boxplot_rectal

stomach_full_with_id <- stomach_merged %>% dplyr::select(X1,
                                                         X2,
                                                         submitter_id.samples) %>% na.omit()

stomach_full_with_id$cluster <- cluster.global %>% filter(cancer == "Stomach Cancer") %>% dplyr::select(cluster) %>% pull()

stomach_cluster_subset <- as.data.frame(stomach)[, which(colnames(stomach) %in% stomach_full_with_id$submitter_id.samples)]

data = as.matrix(stomach_cluster_subset[c(1:60483), ])
rownames(data) <- common.genes
library(scran)

clusters <- cluster.global %>% filter(cancer == "Stomach Cancer") %>% dplyr::select(cluster) %>% pull()

N_DE=6
cluster.markers = scran::findMarkers(x=data, pval.type = "any", 
                                     groups=as.numeric(as.character(clusters)), "binom")

clust_VI_stable = as.numeric(as.character(clusters))
markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = cbind(t(data[markers,]), clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(mean))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Cluster", as.character(unique(clusters))), "FDR")
Table_3[-1,]

Boxplot_DE_stomach <- function(markers_cell = markers_cell){
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  BoxPlot_Data = rbind(cbind(ord_core_stable_mat[,c(1,7)], markers[1]),
                       cbind(ord_core_stable_mat[,c(2,7)], markers[2]),
                       cbind(ord_core_stable_mat[,c(3,7)], markers[3]),
                       cbind(ord_core_stable_mat[,c(4,7)], markers[4]),
                       cbind(ord_core_stable_mat[,c(5,7)], markers[5]),
                       cbind(ord_core_stable_mat[,c(6,7)], markers[6]))
  
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  BoxPlot_Data           = data.frame(BoxPlot_Data)
  BoxPlot_Data$value     = as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster    = factor(BoxPlot_Data$cluster)
  
  ggplot(BoxPlot_Data, aes(y=value, x=gene, col=cluster)) + geom_boxplot() + 
    # labs(title = "Stomach Cancer") +
    theme_bw() + 
    scale_color_manual(values = myvalues) + theme(axis.title=element_blank(), 
                                                  axis.text =element_text(size=12, angle = 45, vjust = 0.5, 
                                                                          hjust=1), strip.text = element_text(size=30),
                                                  plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ))
}

# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_Boxplot_stomach = Boxplot_DE_stomach(markers_cell = markers_cell)
Plot_Boxplot_stomach

eoso_full_with_id <- eoso_merged %>% dplyr::select(cigarettes_per_day.exposures,
                                                   X1,
                                                   X2,
                                                   submitter_id.samples) %>% na.omit()

eoso_full_with_id$cluster <- cluster.global %>% filter(cancer == "Esophageal Cancer") %>% dplyr::select(cluster) %>% pull()

eoso_cluster_subset <- as.data.frame(eoso)[, which(colnames(eoso) %in% eoso_full_with_id$submitter_id.samples)]

data = as.matrix(eoso_cluster_subset[c(1:60483), ])
rownames(data) <- common.genes
library(scran)
clusters <- cluster.global %>% filter(cancer == "Esophageal Cancer") %>% dplyr::select(cluster) %>% pull()
N_DE=6
cluster.markers = scran::findMarkers(x=data, pval.type = "any", 
                                     groups=as.numeric(as.character(clusters)), "binom")

clust_VI_stable = as.numeric(as.character(clusters))
markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = cbind(t(data[markers,]), clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(mean))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Cluster", as.character(unique(clusters))), "FDR")
Table_3[-1,]

Boxplot_DE_eoso <- function(markers_cell = markers_cell){
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  BoxPlot_Data = rbind(cbind(ord_core_stable_mat[,c(1,7)], markers[1]),
                       cbind(ord_core_stable_mat[,c(2,7)], markers[2]),
                       cbind(ord_core_stable_mat[,c(3,7)], markers[3]),
                       cbind(ord_core_stable_mat[,c(4,7)], markers[4]),
                       cbind(ord_core_stable_mat[,c(5,7)], markers[5]),
                       cbind(ord_core_stable_mat[,c(6,7)], markers[6]))
  
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  BoxPlot_Data           = data.frame(BoxPlot_Data)
  BoxPlot_Data$value     = as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster    = factor(BoxPlot_Data$cluster)
  
  ggplot(BoxPlot_Data, aes(y=value, x=gene, col=cluster)) + geom_boxplot() + 
    # labs(title = "Esophageal Cancer") +
    theme_bw() + 
    scale_color_manual(values = myvalues) + theme(axis.title=element_blank(), 
                                                  axis.text =element_text(size=12, angle = 45, vjust = 0.5, 
                                                                          hjust=1), strip.text = element_text(size=30),
                                                  plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ))
}

# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_Boxplot_eoso = Boxplot_DE_eoso(markers_cell = markers_cell)
Plot_Boxplot_eoso

################################################################################
## DE ANALYSIS LOCAL
################################################################################
colon_full_with_id <- colon_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                                     bmi.exposures,
                                                     X1,
                                                     X2,
                                                     submitter_id.samples) %>% na.omit()

colon_full_with_id = colon_full_with_id[-singleton.index_2.2, ]
colon_full_with_id$cluster <- cluster.local %>% filter(cancer == "Colon Cancer") %>% dplyr::select(cluster) %>% pull()

colon_cluster_subset <- as.data.frame(colon)[, which(colnames(colon) %in% colon_full_with_id$submitter_id.samples)]

data = as.matrix(colon_cluster_subset[c(1:60483), ])
rownames(data) <- common.genes
library(scran)
clusters <- cluster.local.LS.colon$cluster_org
clusters_relabelled = cluster.local.LS.colon$cluster 
N_DE=6

cluster.markers = scran::findMarkers(x=data, pval.type = "any", 
                                     groups=as.numeric(as.character(clusters)), "binom")

clust_VI_stable = as.character(clusters_relabelled)
markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = data.frame(t(data[markers,]), clust_VI_stable = clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(median))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Cluster", sort(unique(clust_VI_stable))), "FDR")
Table_3[-1,]

Boxplot_DE <- function(markers_cell = markers_cell){
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  names(ord_core_stable_mat) <- NULL
  BoxPlot_Data = rbind(data.frame(cbind(ord_core_stable_mat[,c(1,7)], X3 = markers[1])),
                       data.frame(cbind(ord_core_stable_mat[,c(2,7)], X3 = markers[2])),
                       data.frame(cbind(ord_core_stable_mat[,c(3,7)], X3 = markers[3])),
                       data.frame(cbind(ord_core_stable_mat[,c(4,7)], X3 = markers[4])),
                       data.frame(cbind(ord_core_stable_mat[,c(5,7)], X3 = markers[5])),
                       data.frame(cbind(ord_core_stable_mat[,c(6,7)], X3 = markers[6])))
  
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  BoxPlot_Data           = data.frame(BoxPlot_Data)
  BoxPlot_Data$value     = as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster    = factor(BoxPlot_Data$cluster)
  
  ggplot(BoxPlot_Data, aes(y=value, x=gene, col=cluster)) + geom_boxplot() + 
    # labs(title = "Colon Cancer") +
    theme_bw() + scale_color_manual(values = myvalues_local.LS.combined) + theme(axis.title=element_blank(), 
                                                                                 axis.text =element_text(size=12, angle = 45, vjust = 0.5, 
                                                                                                         hjust=1), strip.text = element_text(size=30),
                                                                                 plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ))
}
# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_Boxplot_colon_local = Boxplot_DE(markers_cell = markers_cell) 
Plot_Boxplot_colon_local

rectal_full_with_id <- rectal_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                                       X1,
                                                       X2,
                                                       submitter_id.samples) %>% na.omit()

rectal_full_with_id = rectal_full_with_id[-singleton.index_2.1, ]
rectal_full_with_id$cluster <- cluster.local %>% filter(cancer == "Rectal Cancer") %>% dplyr::select(cluster) %>% pull()

rectal_cluster_subset <- as.data.frame(rectal)[, which(colnames(rectal) %in% rectal_full_with_id$submitter_id.samples)]

data = as.matrix(rectal_cluster_subset[c(1:60483), ])
rownames(data) <- common.genes
library(scran)
clusters <- cluster.local.LS.rectal$cluster_org
clusters_relabelled = cluster.local.LS.rectal$cluster 
N_DE=6
cluster.markers = scran::findMarkers(x=data, pval.type = "any", 
                                     groups=as.numeric(as.character(clusters)), "binom")

clust_VI_stable = as.character(clusters_relabelled)
markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = data.frame(t(data[markers,]), clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(mean))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Cluster", sort(unique(clust_VI_stable))), "FDR")
Table_3[-1,]

Boxplot_DE <- function(markers_cell = markers_cell){
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  names(ord_core_stable_mat) <- NULL
  BoxPlot_Data = rbind(data.frame(cbind(ord_core_stable_mat[,c(1,7)], X3 = markers[1])),
                       data.frame(cbind(ord_core_stable_mat[,c(2,7)], X3 = markers[2])),
                       data.frame(cbind(ord_core_stable_mat[,c(3,7)], X3 = markers[3])),
                       data.frame(cbind(ord_core_stable_mat[,c(4,7)], X3 = markers[4])),
                       data.frame(cbind(ord_core_stable_mat[,c(5,7)], X3 = markers[5])),
                       data.frame(cbind(ord_core_stable_mat[,c(6,7)], X3 = markers[6])))
  
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  BoxPlot_Data           = data.frame(BoxPlot_Data)
  BoxPlot_Data$value     = as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster    = factor(BoxPlot_Data$cluster)
  
  ggplot(BoxPlot_Data, aes(y=value, x=gene, col=cluster)) + geom_boxplot() + 
    # labs(title = "Rectal Cancer") +
    theme_bw() + scale_color_manual(values = myvalues_local.LS.combined) + theme(axis.title=element_blank(), 
                                                                                 axis.text =element_text(size=12, angle = 45, vjust = 0.5, 
                                                                                                         hjust=1), strip.text = element_text(size=30),
                                                                                 plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ))
}
# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_Boxplot_rectal_local = Boxplot_DE(markers_cell = markers_cell) 
Plot_Boxplot_rectal_local
################################################################################
## KM ANALYSIS (GLOBAL-LEVEL)
################################################################################
library(data.table)
# READING THE PHENOTYPE DATA
colon_cov = fread("./Real Data Analysis/Data/Colon.gz")
rectal_cov = fread("./Real Data Analysis/Data/Rectal.gz")
stomach_cov = fread("./Real Data Analysis/Data/Stomach.gz")
eoso_cov = fread("./Real Data Analysis/Data/Eoso.gz")

colon_surv = fread("./Real Data Analysis/Data/Colon_Survival.txt")
rectal_surv = fread("./Real Data Analysis/Data/Rectal_Survival.txt")
stomach_surv = fread("./Real Data Analysis/Data/Stomach_Survival.txt")
eoso_surv = fread("./Real Data Analysis/Data/Eoso_Survival.txt")
library(tidyverse)
# RENAMING THE VARIABLES
stomach_merged <- dplyr::rename(stomach_merged, x = X1)
stomach_merged <- dplyr::rename(stomach_merged, y = X2) 

cluster <-  cluster.global %>% filter(cancer == "Stomach Cancer")
stomach_merged2 <- left_join(cluster, stomach_merged, by = "x")
stomach_merged2 <- stomach_merged2 %>% dplyr::select(x, y.x, cluster, submitter_id.samples)

rectal_merged <- dplyr::rename(rectal_merged, x = X1)
rectal_merged <- dplyr::rename(rectal_merged, y = X2) 

cluster <-  cluster.global %>% filter(cancer == "Rectal Cancer")
rectal_merged2 <- left_join(cluster, rectal_merged, by = "x")
rectal_merged2 <- rectal_merged2 %>% dplyr::select(x, y.x, cluster, submitter_id.samples)

colon_merged <- dplyr::rename(colon_merged, x = X1)
colon_merged <- dplyr::rename(colon_merged, y = X2) 

cluster <-  cluster.global %>% filter(cancer == "Colon Cancer")
colon_merged2 <- left_join(cluster, colon_merged, by = "x")
colon_merged2 <- colon_merged2 %>% dplyr::select(x, y.x, cluster, submitter_id.samples)

eoso_merged <- dplyr::rename(eoso_merged, x = X1)
eoso_merged <- dplyr::rename(eoso_merged, y = X2) 

cluster <-  cluster.global %>% filter(cancer == "Esophageal Cancer")
eoso_merged2 <- left_join(cluster, eoso_merged, by = "x")
eoso_merged2 <- eoso_merged2 %>% dplyr::select(x, y.x, cluster, submitter_id.samples)

library(stringr)
colon_merged2$submitter_id.samples <- str_sub(colon_merged2$submitter_id.samples, end = -2)
rectal_merged2$submitter_id.samples <- str_sub(rectal_merged2$submitter_id.samples, end = -2)
stomach_merged2$submitter_id.samples <- str_sub(stomach_merged2$submitter_id.samples, end = -2)
eoso_merged2$submitter_id.samples <- str_sub(eoso_merged2$submitter_id.samples, end = -2)

library(tidyverse)
colon_surv <- dplyr::rename(colon_surv, submitter_id.samples = sample)
rectal_surv <- dplyr::rename(rectal_surv, submitter_id.samples = sample)
stomach_surv <- dplyr::rename(stomach_surv, submitter_id.samples = sample)
eoso_surv <- dplyr::rename(eoso_surv, submitter_id.samples = sample)

stomach_surv$submitter_id.samples <- str_sub(stomach_surv$submitter_id.samples, end = -2)
## MERGING THE SURVIVAL DATA FOR COLON CANCER
colon_merged_surv <- left_join(colon_merged2, colon_surv, by = "submitter_id.samples")
# SELECT ONLY THE RELEVANT VARIABLES FOR KM ANALYSIS
colon_merged_surv2 <- colon_merged_surv %>% dplyr::select(OS, OS.time, cluster)

require("survival")
# FIT THE KM CURVE
fit.colon <- survfit(Surv(OS.time, OS) ~ cluster, data = colon_merged_surv2)

mytheme <- theme_bw() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=24, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=20, face="bold", colour = "black"),    
    axis.title.y = element_text(size=20, face="bold", colour = "black"),    
    axis.text.x = element_text(size=20, face="bold", colour = "black"), 
    axis.text.y = element_text(size=20, face="bold", colour = "black"),
    strip.text.x = element_text(size = 16, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 16, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=18),
    legend.text=element_text(size=17)
  )

library("survminer")
os1.global = ggsurvplot(fit.colon, data = colon_merged_surv2, conf.int = F,
                        legend.title = "Cluster",
                        legend.labs = gtools::mixedsort(as.character(unique(colon_merged_surv2$cluster))),
                        palette = myvalues,
                        ggtheme = mytheme) + ggtitle("Overall Survival: Colon Cancer") 

## MERGING THE SURVIVAL DATA FOR RECTAL CANCER
rectal_merged_surv <- left_join(rectal_merged2, rectal_surv, by = "submitter_id.samples")
# SELECT ONLY THE RELEVANT VARIABLES FOR KM ANALYSIS
rectal_merged_surv2 <- rectal_merged_surv %>% dplyr::select(OS, OS.time, cluster)

require("survival")
# FIT THE KM CURVE
fit.rectal <- survfit(Surv(OS.time, OS) ~ cluster, data = rectal_merged_surv2)

library("survminer")
os2.global = ggsurvplot(fit.rectal, data = rectal_merged_surv2, conf.int = F,
                        legend.title = "Cluster",
                        legend.labs = gtools::mixedsort(as.character(unique(rectal_merged_surv2$cluster))),
                        palette = myvalues,
                        ggtheme = mytheme) + ggtitle("Overall Survival: Rectal Cancer")

## MERGING THE SURVIVAL DATA FOR STOMACH CANCER
stomach_merged_surv <- left_join(stomach_merged2, stomach_surv, by = "submitter_id.samples")
# SELECT ONLY THE RELEVANT VARIABLES FOR KM ANALYSIS
stomach_merged_surv2 <- stomach_merged_surv %>% dplyr::select(OS, OS.time, cluster)

require("survival")
# FIT THE KM CURVE
fit.stomach <- survfit(Surv(OS.time, OS) ~ factor(cluster), data = stomach_merged_surv2)


library("survminer")
os3.global = ggsurvplot(fit.stomach, data = stomach_merged_surv2, conf.int = F,
                        legend.title = "Cluster",
                        legend.labs = gtools::mixedsort(as.character(unique(stomach_merged_surv2$cluster))),
                        palette = myvalues,
                        ggtheme = mytheme
) + ggtitle("Overall Survival: Stomach Cancer")

# SOME ADDITONAL PRE-PROCESSING BEFORE KM ANALYSIS
eoso_surv$submitter_id.samples <- str_sub(eoso_surv$submitter_id.samples, end = -2)
## MERGING THE SURVIVAL DATA FOR ESOPHAGEAL CANCER
eoso_merged_surv <- left_join(eoso_merged2, eoso_surv, by = "submitter_id.samples")
# SELECT ONLY THE RELEVANT VARIABLES FOR KM ANALYSIS
eoso_merged_surv2 <- eoso_merged_surv %>% dplyr::select(OS, OS.time, cluster)

require("survival")
# FIT THE KM CURVE
fit.eoso <- survfit(Surv(OS.time, OS) ~ factor(cluster), data = eoso_merged_surv2)

library("survminer")
os4.global = ggsurvplot(fit.eoso, data = eoso_merged_surv2, conf.int = F,
                        legend.title = "Cluster",
                        legend.labs = gtools::mixedsort(as.character(unique(eoso_merged_surv2$cluster))),
                        palette = myvalues,
                        ggtheme = mytheme
) + ggtitle("Overall Survival: Esophageal Cancer")

if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
splots <- list()
splots[[1]] <- os4.global
splots[[2]] <- os1.global
splots[[3]] <- os3.global
splots[[4]] <- os2.global

res <- arrange_ggsurvplots(splots, print = TRUE,
                           ncol = 2, nrow = 2)
################################################################################
## KM ANALYSIS (LOCAL-LEVEL)
################################################################################
cluster <- cluster.local.LS.combined %>% filter(cancer == "Rectal Cancer")
rectal_merged3 <- left_join(cluster, rectal_merged, by = "x")
rectal_merged3 <- rectal_merged3 %>% dplyr::select(x, y.x, cluster, submitter_id.samples)

cluster <- cluster.local.LS.combined %>% filter(cancer == "Colon Cancer")
colon_merged3 <- left_join(cluster, colon_merged, by = "x")
colon_merged3 <- colon_merged3 %>% dplyr::select(x, y.x, cluster, submitter_id.samples)

library(stringr)
colon_merged3$submitter_id.samples <- str_sub(colon_merged3$submitter_id.samples, end = -2)
rectal_merged3$submitter_id.samples <- str_sub(rectal_merged3$submitter_id.samples, end = -2)

## MERGING THE SURVIVAL DATA FOR COLON CANCER
colon_merged_surv3 <- left_join(colon_merged3, colon_surv, by = "submitter_id.samples")
# SELECT ONLY THE RELEVANT VARIABLES
colon_merged_surv4 <- colon_merged_surv3 %>% dplyr::select(OS, OS.time, cluster)

require("survival")
# FIT THE KM CURVE
fit.colon <- survfit(Surv(OS.time, OS) ~ cluster, data = colon_merged_surv4)

library("survminer")
os1.local = ggsurvplot(fit.colon, data = colon_merged_surv4, conf.int = F,
                       legend.title = "Cluster",
                       legend.labs = as.character(c("12", "20", "3", "7a", "7b", "7c")),
                       palette = myvalues_local.LS.combined,
                       ggtheme = mytheme) + ggtitle("Overall Survival: Colon Cancer")

## MERGING THE SURVIVAL DATA FOR RECTAL CANCER
rectal_merged_surv3 <- left_join(rectal_merged3, rectal_surv, by = "submitter_id.samples")
# SELECT ONLY THE RELEVANT VARIABLES
rectal_merged_surv4 <- rectal_merged_surv3 %>% dplyr::select(OS, OS.time, cluster)

require("survival")
fit.rectal <- survfit(Surv(OS.time, OS) ~ cluster, data = rectal_merged_surv4)
library("survminer")
os2.local = ggsurvplot(fit.rectal, data = rectal_merged_surv4, conf.int = F,
                       legend.title = "Cluster",
                       legend.labs = as.character(c("12", "3", "7a", "7b")),
                       palette = myvalues_local.LS.combined,
                       ggtheme = mytheme
) + ggtitle("Overall Survival: Rectal Cancer")

if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
splots <- list()
splots[[1]] <- os1.local
splots[[2]] <- os2.local

res <- arrange_ggsurvplots(splots, print = TRUE,
                           ncol = 2, nrow = 1)
