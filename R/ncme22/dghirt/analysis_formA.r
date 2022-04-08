require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
require(bayesplot)

################################################################################

# Import data (long format)

formA <- read.csv(here('data/data_ncme2022/formA/formA.csv'))

  # Select the first 100 examinees and first 10 items

  # d <- formA[formA$pid %in% 1:100 & formA$iid %in% 1:10,]


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
  RT             = log(d$time),
  Y              = d$score
)


# Compile the model syntax

mod <- cmdstan_model(here('R/ncme22/dghirt/dghirt.stan'))

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
  View(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary)
  hist(T)
  
  T[which(T>0.8)]
  
  
  View(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary)
  
  C <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
  C
  
  C[which(C>0.8)]
  
  gr <- c(rep('operational',50),rep('pilot',121))
  
  describeBy(C,gr)
  
  plot(density(C[1:50]),type='l',xlim=c(0,1),ylim=c(0,3))
  points(density(C[51:171]),type='l')
  
  
  summary(stanfit, pars = c("mu_beta"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_beta"), probs = c(0.025, 0.975))$summary
  
  summary(stanfit, pars = c("mu_alpha"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_alpha"), probs = c(0.025, 0.975))$summary
  
  summary(stanfit, pars = c("sigma_taut"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_tauc"), probs = c(0.025, 0.975))$summary

  summary(stanfit, pars = c("omega_tau"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("omega_item_rt"), probs = c(0.025, 0.975))$summary
  
  
  View(summary(stanfit, pars = c("tau"), probs = c(0.025, 0.975))$summary)
  
  tau <- matrix(summary(stanfit, pars = c("tau"), probs = c(0.025, 0.975))$summary[,1],
                1000,2,byrow=T)
  describe(tau)
  
  item_rt <- matrix(summary(stanfit, pars = c("item_rt"), probs = c(0.025, 0.975))$summary[,1],
                171,2,byrow=T)
  
  View(summary(stanfit, pars = c("item_rt"), probs = c(0.025, 0.975))$summary)
  
  describe(item_rt)
  

  summary(stanfit, pars = c("mu_b"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_b"), probs = c(0.025, 0.975))$summary
  
  summary(stanfit, pars = c("mu_thetat"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_thetat"), probs = c(0.025, 0.975))$summary
  
  summary(stanfit, pars = c("mu_thetac"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_thetac"), probs = c(0.025, 0.975))$summary

  summary(stanfit, pars = c("omega_theta"), probs = c(0.025, 0.975))$summary
  
  View(summary(stanfit, pars = c("theta"), probs = c(0.025, 0.975))$summary)
  
  theta <- matrix(summary(stanfit, pars = c("theta"), probs = c(0.025, 0.975))$summary[,1],
                1000,2,byrow=T)
  
  describe(theta)
  
  View(summary(stanfit, pars = c("b2"), probs = c(0.025, 0.975))$summary)
  
  
  
  
  
  examinee_info <- read.csv("data/data_ncme2022/formB/examinee_info_FormB.csv")
  
  table(examinee_info$country[1:1000])
  table(examinee_info[which(T>0.95),]$country)
  
  table(examinee_info$test_center[1:1000])
  table(examinee_info[which(T>0.95),]$test_center)
  
  table(examinee_info$modality[1:1000])
  table(examinee_info[which(T>0.95),]$modality)
  
  table(examinee_info$voucher[1:1000])
  table(examinee_info[which(T>0.95),]$voucher)
  
  table(examinee_info$Flag.Condition[1:1000])
  table(examinee_info[which(T>0.95),]$Flag.Condition)
  
  
  par_name <- 'pC[141]'
  mcmc_hist_by_chain(x = stanfit,pars=par_name)
  

