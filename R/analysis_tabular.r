require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)


#######################################################################

# rt is a N  x n tabular response time matrix
# C is a vector of length n for item compromise status

data_rt <- list(
  J               = n,
  I               = N,
  RT              = log(rt[,1:20]),
  C               = rep(c(0,1),10)
  
)

#######################################################################
# Read the Stan model syntax, this is same across all replications

mod <- cmdstan_model(here('R/dglnrt_tabular.stan'))
  
  fit <- mod$sample(
    data = data_rt,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 30,
    iter_sampling = 120,
    refresh = 10,
    adapt_delta = 0.99)
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())

  
  T <- as.numeric(summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary[,1])
  B <- as.numeric(summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary[,1])
  A <- as.numeric(summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary[,1])
  
  tau_t <- as.numeric(summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary[,1])
  tau_c <- as.numeric(summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary[,1])
  
  mean(T[1:30])
  mean(T[31:200])
  
  th <- quantile(T[31:200],.5)
  th = 0.5
  
  table(T[1:30]>th)
  table(T[31:200]>th)
  
  
  
  
  
  
  
  
