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

  
  T <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
  C <- summary(stanfit, pars = c("C"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=T[,1],actual=c(rep(1,75),rep(0,425)))
  
  plot(density(T[76:200,1],adjust=2),xlim=c(0,1))
  points(density(T[1:75,1],adjust = 2),type='l',lty=2)
  
  mean(T[1:75,1])
  mean(T[76:500,1])
  
  table(T[1:75,1]>.7)
  table(T[76:500,1]>.7)
  
  
  C <- summary(stanfit, pars = c("C"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=C[,1],actual=rep(c(0,1),25))
  
  
  
  
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

  data_rt <- list(
    J               = n,
    I               = N,
    R               = r[,1:n])
  
#######################################################################
# Read the Stan model syntax, this is same across all replications
  
  mod <- cmdstan_model(here('R/dgIRT_tabular.stan'))
  
  fit <- mod$sample(
    data = data_rt,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 100,
    iter_sampling = 500,
    refresh = 10,
    adapt_delta = 0.99)
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())
  
  
  T <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
  
  auc_roc(preds=T[,1],actual=c(rep(1,N*pe),rep(0,N-N*pe)))
  
  hist(T[,7])
  
  mean(T[1:(N*pe),1])
  mean(T[(N*pe+1):N,1])
  
  
  plot(density(T[(N*pe+1):N,1],adjust=2),xlim=c(0,1))
  points(density(T[1:(N*pe),1],adjust = 2),type='l',lty=2)
  
  
  table(T[1:(N*pe),1]>.6)
  table(T[(N*pe+1):N,1]>.6)
  
  
  C <- summary(stanfit, pars = c("C"), probs = c(0.025, 0.975))$summary
  hist(C[,7])
  
  auc_roc(preds=C[,1],actual=rep(c(0,1),n/2))
  
  
  
  B <- as.numeric(summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary[,1])
  cor(b,B)
  plot(b,B)
  
  
  
  A <- as.numeric(summary(stanfit, pars = c("a"), probs = c(0.025, 0.975))$summary[,1])
  cor(a,A)
  plot(a,A)
  
  theta_t <- as.numeric(summary(stanfit, pars = c("theta_t"), probs = c(0.025, 0.975))$summary[,1])
  cor(theta,theta_t)
  plot(theta,theta_t)
  
  
  theta_cc <- as.numeric(summary(stanfit, pars = c("theta_c"), probs = c(0.025, 0.975))$summary[,1])
  cor(theta_c,theta_cc)
  plot(theta_c,theta_cc)

  
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
  iter_warmup   = 500,
  iter_sampling = 3000,
  refresh = 10,
  adapt_delta = 0.99)

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())


T <- summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary

hist(T[,7])
mean(T[1:40,1])
mean(T[41:200,1])
auc_roc(preds=T[,1],actual=c(rep(1,40),rep(0,160)))


plot(density(T[41:200,1],adjust=2),xlim=c(0,1))
points(density(T[1:40,1],adjust = 2),type='l',lty=2)

table(T[1:40,1]>.65)
table(T[41:200,1]>.65)




C <- summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary

hist(C[,7])
auc_roc(preds=C[,1],actual=rep(c(0,1),15))

C


summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("omega_I"), probs = c(0.025, 0.975))$summary

A <- exp(summary(stanfit, pars = c('item'), probs = c(0.025, 0.975))$summary[seq(1,120,4),1])
cor(a,A)
plot(a,A)


B <- summary(stanfit, pars = c('item'), probs = c(0.025, 0.975))$summary[seq(2,120,4),1]
cor(b,B)
plot(b,B)

Al <- exp(summary(stanfit, pars = c('item'), probs = c(0.025, 0.975))$summary[seq(3,120,4),1])
cor(alpha,Al)
plot(alpha,Al)

Beta <- summary(stanfit, pars = c('item'), probs = c(0.025, 0.975))$summary[seq(4,120,4),1]
cor(beta,Beta)
plot(beta,Beta)














  
