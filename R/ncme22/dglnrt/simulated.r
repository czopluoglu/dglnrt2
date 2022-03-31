require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)

################################################################################
require(MASS)

set.seed(06202021)

N = 200    # number of examinees
n = 20     # number of items

# Time intensity parameters

beta  <- rnorm(n,4,.5)

# Time discrimination parameters

alpha <- rnorm(n,2,0.5) 

# Tau_t and tau_c

tau <- mvrnorm(N,
               mu = c(0,0.4),
               Sigma = matrix(c(0.01,0.0105,0.0105,0.0225),2,2))

cov2cor(matrix(c(0.01,0.0105,0.0105,0.0225),2,2))

tau_t <- tau[,1]
tau_c <- tau[,2]

# Randomly select (approximately) 20% of examinees as having item prekowledge

H <- rbinom(N,1,.2)

# Randomly select (approximately) 50% of items as compromised

C <- rbinom(n,1,.5)

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

# Reshape it to long format (for plotting purposes)

rt.long <- reshape(data        = rt,
                   idvar       = 'id',
                   varying     = list(colnames(rt)[1:n]),
                   timevar     = "Item",
                   times       = 1:n,
                   v.names      = "RT",
                   direction   = "long")

# Add item status

rt.long$compromised <- NA

for(j in 1:n){
  
  rt.long[rt.long$Item==j,]$compromised = C[j]
  
}

head(rt.long)


d <- rt.long
################################################################################



data_rt <- list(
  J              = length(unique(d$Item)),
  I              = length(unique(d$id)),
  n_obs          = nrow(d),
  p_loc          = d$id,
  i_loc          = d$Item,
  Y              = log(d$RT)
)


# Compile the model syntax

mod <- cmdstan_model(here('R/ncme22/dglnrt/dglnrt.stan'))

# Fit the model

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup   = 250,
  iter_sampling = 750,
  refresh = 10,
  adapt_delta = 0.99)

# Compile the output files into an rstan object


fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

# Estimation Time 

get_elapsed_time(stanfit)

(sum(get_elapsed_time(stanfit))/4)/3600

# Analyze the parameter estimates


T <- as.numeric(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary[,1])

table(H,T>.5)


C_ <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
table(C,C_>0.5)


tau_ <- summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary

tau_ <- matrix(tau_[,1],ncol=2,byrow=TRUE)

cor(tau[,1],tau_[,1])
plot(tau[,1],tau_[,1])


cor(tau[,2],tau_[,2])
plot(tau[,2],tau_[,2])

describe(tau_)

summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary



ipar <- summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary

ipar <- matrix(ipar[,1],ncol=2,byrow=TRUE)
ipar[,1] <- exp(ipar[,1])
  
cor(beta,ipar[,2])
plot(beta,ipar[,2])

cor(alpha,ipar[,1])
plot(alpha,ipar[,1])

describe(ipar)









