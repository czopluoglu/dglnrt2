require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)

################################################################################

# Import data (long format)
  
  d <- read.csv(here('data/data_ncme2022/formB/formB_sub.csv'))
  
  unique(d$pid)
  unique(d$iid)
  
  data_rt <- list(
    J              = 15,
    I              = 150,
    n_obs          = nrow(d),
    p_loc          = d$pid,
    i_loc          = d$iid,
    Y              = log(d$time)
  )
  

# Compile the model syntax

  mod <- cmdstan_model(here('R/ncme22/dglnrt.stan'))

# Fit the model
  
  fit <- mod$sample(
    data = data_rt,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 50,
    iter_sampling = 100,
    refresh = 10,
    adapt_delta = 0.99)

# Compile the output files into an rstan object
  
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())

# Estimation Time 
  
  (sum(get_elapsed_time(stanfit))/4)/3600
  
# Analyze the parameter estimates
  
  
  T <- as.numeric(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary[,1])
  T
  
  C <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
  C
  
  summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary
  
  
