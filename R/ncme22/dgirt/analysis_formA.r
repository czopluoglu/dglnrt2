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

data_resp <- list(
  J              = length(unique(d$iid)),
  I              = length(unique(d$pid)),
  n_obs          = nrow(d),
  p_loc          = d$pid,
  i_loc          = d$iid,
  Y              = d$score
)


# Compile the model syntax

mod <- cmdstan_model(here('R/ncme22/dgirt/dgirt.stan'))

# Fit the model

fit <- mod$sample(
  data = data_resp,
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


  C_ <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])

  est.b <- summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary

  summary(stanfit, pars = c("mu_b"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_b"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary


  est.th <- matrix(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary[,1],
                   nrow=200,ncol=2,byrow=T)

################################################################################

par_name1 <- 'person[5,1]'
par_name2 <- 'person[5,2]'

mcmc_hist_by_chain(x = stanfit,pars=par_name1)
mcmc_hist_by_chain(x = stanfit,pars=par_name2)

par_name <- 'b[8]'
mcmc_hist_by_chain(x = stanfit,pars=par_name)

