###############################################################################
###############################################################################
# Simulate item response and response time data with item preknowledge effect
###############################################################################
###############################################################################

require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)

###############################################################################

set.seed(06152021)

# This is a code being used for simulating a data with item preknowledge. 
# N   : sample size
# n   : number of items
# pe  : proportion of examinees with item preknowledge
# pi  : proportion of compromised items
# d1  : item preknowledge effect on response accuracy
# d2  : item preknowledge effect on response time

N <- 200
n <- 30
pe <- 0.20
pi <- 0.50
d1 <- 1
d2 <- 1.5


  # 100*(1 - 1/exp(d1)), average % reduction in response time due to preknowledge
  # 100*(exp(d2)-1), average % increase in odds of getting the item correct

################################################################################
#            MODEL PARAMETER GENERATION
#
# These parameters are reported in one of the tables in the paper
################################################################################

  # MODEL PARAMETERS
  
  # Response time model item parameters

  beta  <- rnorm(n,4.19,0.38)
  alpha <- rnorm(n,1.42,0.37)
  
  # Response model item parameters
  
  a     <- rlnorm(n,0,.5)
  b     <- rnorm(n,0,1)
  
  # Sigma, correlation matrix among tau and theta
  
  Sigma <-  matrix(c(1,.3,.3,1),2,2)
  
  # Standard deviations for tau and theta

  sd_tau   <- 0.2
  sd_theta <- 1
  
  # Means for tau and theta
  
  mu_tau   <- 0
  mu_theta <- 0
  
  # Generate tau and theta for honest responses
  
  S <- cor2cov(rho=Sigma,sigma = c(sd_tau,sd_theta))
  
  
  tau_theta <- mvrnorm(N,c(mu_tau,mu_theta),S)
  
  tau   <- tau_theta[,1]
  theta <- tau_theta[,2]
  
  # Generate tau and theta for fraudulent responses
  
  tau_c   <- tau   + rnorm(N,d1,.1)
  theta_c <- theta + rnorm(N,d2,.1)
  
  # A vector for item status (0: not disclosed, 1:disclosed)
  
  C1    <- rep(c(0,1),n/2)
  C2    <- rep(0,n)
  
  # RESPONSE TIME GENERATION
  
  # Note that the response time data generated is
  # already on the log scale
  
  # Examinees with Item Preknowledge
  
  rt0 <- matrix(nrow = round(N*pe), ncol = n)
  
  for (i in 1:nrow(rt0)) {
    for (j in 1:n) {
      p_t = beta[j] - tau[i]
      p_c = beta[j] - tau_c[i]
      p   = p_t * (1 - C1[j]) + p_c * C1[j]
      rt0[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Examinees with no Item Preknowledge
  
  rt1 <- matrix(nrow = N-round(N*pe), ncol = n)
  
  for (i in 1:nrow(rt1)) {
    for (j in 1:n) {
      p_t = beta[j] - tau[i]
      p_c = beta[j] - tau_c[i]
      p   = p_t * (1 - C2[j]) + p_c * C2[j]
      rt1[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }

  # Combine the groups
  
  rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
              cbind(data.frame(exp(rt1)), gr = 0))
  
  colnames(rt)[1:n] <- paste0('RT',1:n)
  
  
  rt$ID <- 1:nrow(rt)
  
  # RESPONSE GENERATION
  
  # Examinees with Item Preknowledge
  
  r0 <- matrix(nrow = round(N*pe), ncol = n)
  
  for (i in 1:nrow(r0)) {
    for (j in 1:n) {
      p_t = theta[i] - b[j]
      p_c = theta_c[i] - b[j]
      p   = p_t * (1 - C1[j]) + p_c * C1[j]
      
      prob = exp(a[j]*p)/(1+exp(a[j]*p))
      
      r0[i, j] = (prob>runif(1,0,1))*1
    }
  }
  
  # Examinees with no Item Preknowledge

  r1 <- matrix(nrow = N-round(N*pe), ncol = n)
  
  for (i in 1:nrow(r1)) {
    for (j in 1:n) {
      p_t = theta[i] - b[j]
      p_c = theta_c[i] - b[j]
      p   = p_t * (1 - C2[j]) + p_c * C2[j]
      
      prob = exp(a[j]*p)/(1+exp(a[j]*p))
      
      r1[i, j] = (prob>runif(1,0,1))*1
    }
  }
  
  
  # Combine the groups
  
  r <- rbind(cbind(data.frame(r0), gr = 1),
              cbind(data.frame(r1), gr = 0))
  
  
  colnames(r)[1:n] <- paste0('R',1:n)
  
  r$ID <- 1:nrow(r)
  
  
# Combine Response Time and Item Response Data into one
  
  resp <- merge(rt,r,by=c('ID','gr'))
  
  
##############################################################################
  
  vnames <- NULL
  for(i in 3:(n+2)){
    vnames <- c(vnames,c(i,i+n))
  }
  
# Prepare response time data in the long format
  
  Y.long <- reshape(
    data        = resp[,-2],
    idvar       = "ID",
    varying     = colnames(resp)[vnames],
    timevar     = "Item",
    times       = 1:n,
    v.names      = c("R",'RT'),
    direction   = "long"
  )
  
  Y.long <- Y.long[order(Y.long$ID),]
  
  colnames(Y.long) <- c('ID','Item','RT','R') 
  
  Y.long$RT <- log(Y.long$RT)
  
  
# Remove missing data (relevant for real data analysis)
  
  Y.long <- na.omit(Y.long)
  
  beta0 = mean(Y.long$RT)
  
  alpha0 = sqrt(1/colMeans((log(rt[,1:n]) - matrix(beta0,N,n,byrow=T))^2,na.rm=TRUE))
  
# Input data list for Stan
  
  data_rt <- list(
    J               = n,
    I               = N,
    n_obs           = nrow(Y.long),
    ind_person_obs  = Y.long$ID,
    ind_item_obs    = Y.long$Item,
    RT              = Y.long$RT,
    R               = Y.long$R,
    beta0           = mean(Y.long$RT),
    v1              = N/2,
    v2              = N/2*mean(alpha0)
  )
  
##############################################################################


