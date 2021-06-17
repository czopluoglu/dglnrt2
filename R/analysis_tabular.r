require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
require(mltools)

#######################################################################

# rt is a N  x n tabular response time matrix
# C is a vector of length n for item compromise status

data_rt <- list(
  J               = n,
  I               = N,
  RT              = log(rt[,1:n])
)

#######################################################################
# Read the Stan model syntax, this is same across all replications

mod <- cmdstan_model(here('R/dglnrt_tabular.stan'))
  
  fit <- mod$sample(
    data = data_rt,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 200,
    iter_sampling = 1000,
    refresh = 10,
    adapt_delta = 0.99)
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())

  
  T <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
  C <- summary(stanfit, pars = c("C"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=T[,1],actual=c(rep(1,75),rep(0,425)))
  
  mean(T[1:75,1])
  mean(T[76:500,1])
  
  table(T[1:75,1]>.7)
  table(T[76:500,1]>.7)
  
  
  C <- summary(stanfit, pars = c("C"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=C[,1],actual=rep(c(0,1),25))
  
  
  
  B <- as.numeric(summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary[,1])
  A <- as.numeric(summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary[,1])
  
  tau_t <- as.numeric(summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary[,1])
  tau_c <- as.numeric(summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary[,1])
  
  
  
  
  
  
  
