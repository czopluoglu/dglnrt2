require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)

#######################################################################
# Read the Stan model syntax, this is same across all replications

mod <- cmdstan_model(here('R/dglnrt2.stan'))
  
  fit <- mod$sample(
    data = data_rt,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 500,
    iter_sampling = 1000,
    refresh = 10,
    adapt_delta = 0.99)
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())

  
  T <- as.numeric(summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary[,1])
  C <- as.numeric(summary(stanfit, pars = c("C"), probs = c(0.025, 0.975))$summary[,1])
  as.numeric(summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary[,1])
  
  
