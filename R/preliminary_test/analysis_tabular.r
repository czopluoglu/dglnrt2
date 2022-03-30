require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
require(mltools)

#######################################################################
#######################################################################
#######################################################################
#
#      RESPONSE TIME ONLY
#
#######################################################################
#######################################################################
#######################################################################

# rt is an N  x n tabular response time matrix

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

  
  pH <- summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary
  pC <- summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=pH[,1],actual=rt$group)
  
  plot(density(pH[rt$group==0,1],adjust=2),xlim=c(0,1))
  points(density(pH[rt$group==1,1],adjust = 2),type='l',lty=2)

  
  pC
  auc_roc(preds=pC[,1],actual=C)

  plot(density(pC[C==1,1],adjust=2),xlim=c(0,1))
  points(density(pC[C==0,1],adjust = 2),type='l',lty=2)
  
  mean(pC[C==0,1])
  mean(pC[C==1,1])
  
  
  
  
  
  B <- as.numeric(summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary[,1])
  
  cor(B,b)
  plot(B,b)
  describe(b)
  describe(B)
  
  A <- as.numeric(summary(stanfit, pars = c("a"), probs = c(0.025, 0.975))$summary[,1])
  cor(A,a)
  plot(A,a)
  describe(a)
  describe(A)
  
  
  theta_t <- as.numeric(summary(stanfit, pars = c("theta_t"), probs = c(0.025, 0.975))$summary[,1])
  theta_c <- as.numeric(summary(stanfit, pars = c("theta_c"), probs = c(0.025, 0.975))$summary[,1])

  plot(theta,theta_t)
  cor(theta,theta_t)
  
#######################################################################
#######################################################################
#######################################################################
#
#      ITEM RESPONSES ONLY
#
#######################################################################
#######################################################################
#######################################################################
  
# r is an N  x n tabular item response matrix

  data_r <- list(
    J               = n,
    I               = N,
    R               = r[,1:n])
  
#######################################################################
# Read the Stan model syntax, this is same across all replications
  
  mod <- cmdstan_model(here('R/dgIRT_tabular.stan'))
  
  fit <- mod$sample(
    data = data_r,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 200,
    iter_sampling = 1000,
    refresh = 10,
    adapt_delta = 0.99)
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())
  
  
  pH <- summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=pH[,1],actual=H)
  
  hist(pH[,7])
  
  mean(pH[H==0,1])
  mean(pH[H==1,1])
  
  
  plot(density(pH[H==1,1],adjust=2),xlim=c(0,1))
  points(density(pH[H==0,1],adjust = 2),type='l',lty=2)
  
  
  pC <- summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary
  hist(pC[,7])
  
  auc_roc(preds=pC[,1],actual=C)
  
  
  require(bayesplot)
  
  mcmc_hist_by_chain(x = stanfit,pars='pC[29]')
  mcmc_hist(x = stanfit,pars='pC[29]')
  
  thetas <- summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary
  hist(thetas[,7])
  
  thetas <- matrix(thetas[,1],nrow=200,2,byrow=TRUE)
  
  plot(thetas[H==1,],xlim=c(-1,1),ylim=c(-1,1))
  cor(thetas[H==1,])
  colMeans(thetas[H==1,])
  
  
  plot(thetas[H==0,],xlim=c(-1,1),ylim=c(-1,1))
  cor(thetas[H==0,])
  colMeans(thetas[H==0,])
  
  
  
  B <- as.numeric(summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary[seq(2,60,2),1])
  cor(b,B)
  plot(b,B)
  
  
  
  A <- exp(as.numeric(summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary[seq(1,60,2),1]))
  cor(a,A)
  plot(a,A)
  


  
#######################################################################
#######################################################################
#######################################################################
#
#      RESPONSE TIME and ITEM RESPONSES TOGETHER 
#
#######################################################################
#######################################################################
#######################################################################
  
# r  is an N  x n tabular item response matrix
# rt is an N  x n tabular log response time matrix
  
data_rt <- list(
    J               = n,
    I               = N,
    R               = r[,1:n],
    RT              = log(rt[,1:n])
  )

#######################################################################

mod <- cmdstan_model(here('R/dghirt_tabular.stan'))

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


pH <- summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary

hist(pH[,7])
mean(pH[H==0,1])
mean(pH[H==1,1])
auc_roc(preds=pH[,1],actual=H)

plot(density(pH[H==1,1],adjust=2),xlim=c(0,1))
points(density(pH[H==0,1],adjust = 2),type='l',lty=2)

table(H,pH[,1]>.7)



pC <- summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary

hist(pC[,7])
auc_roc(preds=pC[,1],actual=C)

plot(density(pC[C==1,1],adjust=2),xlim=c(0,1))
points(density(pC[C==0,1],adjust = 2),type='l',lty=2)

mean(pC[C==0,1])
mean(pC[C==1,1])




summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("omega_I"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("omega_P2"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("omega_I2"), probs = c(0.025, 0.975))$summary



Al <- exp(summary(stanfit, pars = c('item'), probs = c(0.025, 0.975))$summary[seq(1,60,2),1])
cor(alpha,Al)
plot(alpha,Al)

Beta <- summary(stanfit, pars = c('item'), probs = c(0.025, 0.975))$summary[seq(2,60,2),1]
cor(beta,Beta)
plot(beta,Beta)



A <- exp(summary(stanfit, pars = c('item2'), probs = c(0.025, 0.975))$summary[seq(1,60,2),1])
cor(a,A)
plot(a,A)


B <- summary(stanfit, pars = c('item2'), probs = c(0.025, 0.975))$summary[seq(2,60,2),1]
cor(b,B)
plot(b,B)











  
