# // Stan model syntax with modified deterministic gated item response (Rasch version)
# 
# data{
#   int <lower=1> I;                       // number of examinees          
#   int <lower=1> J;                       // number of items
#   int <lower=1> n_obs;                   // number of observations (I xJ - missing responses)
#   int <lower=1> p_loc[n_obs];            // person indicator   
#   int <lower=1> i_loc[n_obs];            // item indicator
#   int <lower=0,upper=1> Y[n_obs];       // vector of item responses
# }
# 
# parameters {
#   real mu_b;                 // mean for item difficulty parameters
#   real<lower=0> sigma_b;     // sd for item difficulty parameters
#   
#   real<lower=0> sigma_thetat;// sd for theta_t
#   
#   real mu_thetac;            // mean for theta_c
#   real<lower=0> sigma_thetac;// sd for theta_c
#   
#   corr_matrix[2] omega_P;    // 2 x 2 correlation matrix for person parameters
#   
#   vector<lower=0,upper=1>[J] pC; // vector of length J for the probability of item compromise status
#   
#   vector<lower=0,upper=1>[I] pH; // vector of length I for the probability of examinee item peknowledge 
#   
#   ordered[2] person[I];          // an array with length I for person specific latent parameters
#   // Each array has two elements
#   // first element is theta_t
#   // second element is theta_c
#   // ordered vector assures that tau_c > tau_t for every person
#   // to make sure chains are exploring the same mode and 
#   // multiple chains do not go east and west leading multi-modal posteriors
#   
#   
#   vector[J] b;    // vector of item difficulty parameters
#   
# }
# 
# 
# transformed parameters{
#   
#   vector[2] mu_P;                        // vector for mean vector of person parameters 
#   vector[2] scale_P;                     // vector of standard deviations for person parameters
#   cov_matrix[2] Sigma_P;                 // covariance matrix for person parameters
#   
#   mu_P[1] = 0;
#   mu_P[2] = mu_thetac;
#   
#   scale_P[1] = sigma_thetat;               
#   scale_P[2] = sigma_thetac;
#   
#   Sigma_P = quad_form_diag(omega_P, scale_P); 
# }
# 
# 
# 
# model{
#   
#   sigma_thetat  ~ exponential(1);
#   sigma_thetac  ~ exponential(1);
#   
#   omega_P       ~ lkj_corr(1);
#   mu_thetac     ~ normal(0,1);
#   person        ~ multi_normal(mu_P,Sigma_P);
#   
#   mu_b      ~ normal(0,1);
#   sigma_b   ~ exponential(1);
#   b         ~ normal(mu_b,sigma_b);
#   
#   pC ~ beta(1,1);
#   pH ~ beta(1,1);
#   
#   for (i in 1:n_obs) {
#     
#     // b[i_loc[i]] represents b-parameter of the (i_loc[i])th item
#     
#     //person[p_loc[i],1] represents parameter theta_t of the (p_loc[i])th person
#     //person[p_loc[i],2] represents parameter theta_c of the (p_loc[i])th person
#     
#     
#     real p_t = person[p_loc[i],1] - b[i_loc[i]];  // non-cheating response
#     real p_c = person[p_loc[i],2] - b[i_loc[i]];  // cheating response
#     
#     // log of probability densities for each combination of two discrete parameters
#     // (C,T) = {(0,0),(0,1),(1,0),(1,1)}
#     
#     real lprt1 = log1m(pC[i_loc[i]]) + log1m(pH[p_loc[i]]) + bernoulli_logit_lpmf(Y[i] | p_t);  // T = 0, C=0
#     real lprt2 = log1m(pC[i_loc[i]]) + log(pH[p_loc[i]])   + bernoulli_logit_lpmf(Y[i] | p_t);  // T = 1, C=0
#     real lprt3 = log(pC[i_loc[i]])   + log1m(pH[p_loc[i]]) + bernoulli_logit_lpmf(Y[i] | p_t);  // T = 0, C=1
#     real lprt4 = log(pC[i_loc[i]])   + log(pH[p_loc[i]])   + bernoulli_logit_lpmf(Y[i] | p_c);  // T = 1, C=1 
#     
#     target += log_sum_exp([lprt1, lprt2, lprt3, lprt4]);
#   }
#   
# }

################################################################################
require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
require(irtoys)
require(bayesplot)

set.seed(03302022)
################################################################################

N = 200    # number of examinees
n = 50     # number of items

# Item difficulty parameters

b <- rnorm(n,0,1)

# Theta_t and theta_c

th <- mvrnorm(N,
               mu = c(0,1),
               Sigma = matrix(c(1,0.5,0.5,1),2,2))

th_t <- th[,1]
th_c <- th[,2]

# Randomly select (approximately) 30% of examinees as having item prekowledge

H <- rbinom(N,1,.5)

# Randomly select (approximately) 50% of items as compromised

C <- rbinom(n,1,.5)

# Generate observed responses according to the Rasch model

r <- matrix(nrow=N,ncol=n)

for(i in 1:N){
  for(j in 1:n){
    
    p_t <- exp(th_t[i] - b[j])/(1+exp(th_t[i] - b[j]))
    p_c <- exp(th_c[i] - b[j])/(1+exp(th_c[i] - b[j]))
    
    if(H[i] == 1 & C[j] == 1){
      r[i,j] = rbinom(1,1,p_c)
    } else {
      r[i,j] = rbinom(1,1,p_t)
    }
    
  }
}

# Convert it to data frame and add group membership and a unique ID

r       <- as.data.frame(r)
r$group <- H
r$id    <- 1:nrow(r)

# Check the data

head(r)

# Reshape it to long format (for plotting purposes)

r.long <- reshape(data        = r,
                   idvar       = 'id',
                   varying     = list(colnames(r)[1:n]),
                   timevar     = "Item",
                   times       = 1:n,
                   v.names      = "R",
                   direction   = "long")

# Add item status

r.long$compromised <- NA

for(j in 1:n){
  
  r.long[r.long$Item==j,]$compromised = C[j]
  
}


d <- r.long

################################################################################

  
  data_resp <- list(
    J              = length(unique(d$Item)),
    I              = length(unique(d$id)),
    n_obs          = nrow(d),
    p_loc          = d$id,
    i_loc          = d$Item,
    Y              = d$R
  )
  

# Compile the model syntax

  mod <- cmdstan_model(here('R/ncme22/dgirt/dgirt.stan'))

# Fit the model

# Compile the output files into an rstan object
  
  
  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())

# Estimation Time 
  
  get_elapsed_time(stanfit)
  
  (sum(get_elapsed_time(stanfit))/4)/3600
  
# Analyze the parameter estimates
  
  
  T <- as.numeric(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary[,1])
  
  table(H,T>.6)
  
  
  C_ <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
  table(C,C_>0.5)
  
  est.b <- summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary
  plot(b,est.b[,1])
  cor(b,est.b[,1])
  
  summary(stanfit, pars = c("mu_b"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("sigma_b"), probs = c(0.025, 0.975))$summary
  summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary

  
  est.th <- matrix(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary[,1],
                  nrow=200,ncol=2,byrow=T)

  plot(th[,1],est.th[,1])
  cor(th[,1],est.th[,1])
  
  mean(est.th[,1])
  sd(est.th[,1])
  
  mean(th[,1])
  sd(th[,1])
  
  plot(th[,2],est.th[,2])
  cor(th[,2],est.th[,2])
  
  mean(est.th[,2])
  sd(est.th[,2])
  
  mean(th[,2])
  sd(th[,2])
  
  par_name1 <- 'person[5,1]'
  par_name2 <- 'person[5,2]'
  
  mcmc_hist_by_chain(x = stanfit,pars=par_name1)
  mcmc_hist_by_chain(x = stanfit,pars=par_name2)
  
  par_name <- 'b[8]'
  mcmc_hist_by_chain(x = stanfit,pars=par_name)
  


  
  
  
  
  
  