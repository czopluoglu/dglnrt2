require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)

################################################################################

# Import data (long format)
  
  formA <- read.csv(here('data/data_ncme2022/formA/formA.csv'))
  
  # Select the first 1000 examinees
  
  d <- formA[formA$pid %in% 1:1000,]
  
  length(unique(d$pid))
  length(unique(d$iid))
  
################################################################################
  
  
  data_rt <- list(
    J              = length(unique(d$iid)),
    I              = length(unique(d$pid)),
    n_obs          = nrow(d),
    p_loc          = d$pid,
    i_loc          = d$iid,
    Y              = log(d$time)
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
  T
  
  C <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
  C
  
  summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary
  
  
