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
  
  # d <- formA
  
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
    iter_warmup   = 200,
    iter_sampling = 500,
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
  T
  
  T[which(T>0.95)]
  
  
  View(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary)
  
  C <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
  C
  
  C[which(C>0.95)]
  
  gr <- c(rep('operational',50),rep('pilot',121))
  
  describeBy(C,gr)
  
  plot(density(C[1:50]),type='l',xlim=c(0,1))
  points(density(C[51:171]),type='l')
  
  
  summary(stanfit, pars = c("mu_beta"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_beta"), probs = c(0.025, 0.975))$summary
  
  summary(stanfit, pars = c("sigma_taut"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_tauc"), probs = c(0.025, 0.975))$summary
  
  summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("omega_I"), probs = c(0.025, 0.975))$summary
  
  
  View(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary)
  
  tau <- matrix(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary[,1],
         1000,2,byrow=T)
  
  describe(tau[which(T>0.95),])
  
  examinee_info <- read.csv("data/data_ncme2022/formA/examinee_info_FormA.csv")
  
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
  
  
  
  
  View(summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary)
  
  
