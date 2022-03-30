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
d1 <- 0.5
H <- rbinom(N,1,pe)
C <- rbinom(n,1,pi)


  # 100*(1 - 1/exp(d1)), average % reduction in response time due to preknowledge
  # 100*(exp(d2)-1), average % increase in odds of getting the item correct

################################################################################
#            MODEL PARAMETER GENERATION
#
# These parameters are reported in one of the tables in the paper
################################################################################

################## RESPONSE TIME #############################################

# Time intensity parameters

beta  <- rnorm(n,4,.5)

# Time discrimination parameters

alpha <- rnorm(n,2,0.5) 

# Tau_t and tau_c

tau <- mvrnorm(N,
               mu = c(0,d1),
               Sigma = matrix(c(0.01,0.0105,0.0105,0.0225),2,2))

tau_t <- tau[,1]
tau_c <- tau[,2]

# Generate observed response times according to the model

rt <- matrix(nrow=N,ncol=n)

for(i in 1:N){
  for(j in 1:n){
    
    p_t <- beta[j] - tau_t[i]
    p_c <- beta[j] - tau_c[i]
    
    if(H[i] == 1 & C[j] == 1){
      rt[i,j] = exp(rnorm(1,p_c,1/alpha[j]))
    } else {
      rt[i,j] = exp(rnorm(1,p_t,1/alpha[j]))
    }
    
  }
}

# Convert it to data frame and add group membership and a unique ID

rt       <- as.data.frame(rt)
rt$group <- H
rt$id    <- 1:nrow(rt)

# Check the data

head(rt)

################################################################################
  
# RESPONSE GENERATION
  

# Time intensity parameters

b <- rnorm(n,0,1)

# Time discrimination parameters

a <- rlnorm(n,0,0.5) 

# Theta

theta <- rnorm(N,0,1)

# Generate responses

r <- matrix(nrow = N, ncol = n)
  
for (i in 1:N) {
  for (j in 1:n) {
      
    p    = theta[i] - b[j]
      
    prob = exp(a[j]*p)/(1+exp(a[j]*p))
      
    if(H[i] == 1 & C[j] == 1){
      r[i,j]  = rbinom(1,1,.9)
      } else {
        r[i,j] = (prob>runif(1,0,1))*1
      }
      
  }
}
  
# Convert it to data frame and add group membership and a unique ID

r       <- as.data.frame(r)
r$group <- H
r$id    <- 1:nrow(r)

# Check the data

head(r)
  
# Combine Response Time and Item Response Data into one
  
  resp <- merge(rt,r,by=c('id','group'))
  
  
##############################################################################
  
describeBy(rt[,1],rt$group)
describeBy(rt[,2],rt$group)
describeBy(rt[,3],rt$group)
describeBy(rt[,4],rt$group)
describeBy(rt[,5],rt$group)
describeBy(rt[,6],rt$group)


describeBy(r[,1],r$group)
describeBy(r[,2],r$group)
describeBy(r[,3],r$group)
describeBy(r[,4],r$group)
describeBy(r[,5],r$group)
describeBy(r[,6],r$group)

  
  
