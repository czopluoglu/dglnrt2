//This doesn't work
// have to figure out marginalizing the discrete parameters in Stan
// C ~ bernoulli(pr) can't happen in Stan
// Stan doesn't sample discrete parameters

// https://study.sagepub.com/sites/default/files/chapter16.pdf, page 10
// https://discourse.mc-stan.org/t/assigning-bernoulli-prior-to-missing-entries-in-covariate/10219/2
// https://mc-stan.org/docs/2_27/stan-users-guide/vectorizing-mixtures.html
// https://elevanth.org/blog/2018/01/29/algebra-and-missingness/

data{
  int <lower=1> I;                       // number of examinees          
  int <lower=1> J;                       // number of items
  int <lower=1> n_obs;                   // number of observations
  int <lower=1> ind_person_obs[n_obs];   // person position indicator
  int <lower=1> ind_item_obs[n_obs];     // item position indicator
  real RT[n_obs];                        // RT
  int <lower=0,upper=1> R[n_obs];        // item responses
}

parameters {
  vector[J] beta;                        // vector for time intensity parameters
  real<lower=0> sigma_beta;              // sd for time intensity parameters
  
  vector <lower=0> [J]  alpha;           // vector for time discrimination parameters
  real<lower=0> sigma_alpha;              // sd for time discrimination parameters

  vector[I] tau_t;                       // vector for latent speed for honest responses
  real<lower=0> sigma_taut;              // sd for tau_t

  vector[I] tau_c;                       // vector for latent speed for fraudulent responses
  real<lower=0> sigma_tauc;              // sd for tau_c
  
  vector[I] theta_t;                     // vector for latent ability for honest responses
  real<lower=0> sigma_thetat;            // sd for thetat

  vector[I] theta_c;                     // vector for latent ability for fraudulent responses
  real<lower=0> sigma_thetac;            // sd for thetac
  
  vector[J] b;                           // vector for item difficulty parameters
  real<lower=0> sigma_b;                 // sd for item difficulty parameters
  
  vector<lower = 0>[J] a;                // vector for item discrimination parameters
  real<lower=0> sigma_a;                 // sd for item discrimination parameters

  vector<lower=0,upper=1>[J] C;          // vector for the probability of 
                                         // item compromise status

  vector<lower=0,upper=1>[I] T;          // vector for the probability of examinee 
                                         // item peknowledge 
                                         
  vector[4] mu_P;                        // vector for mean vector of person parameters 
  vector[4] mu_I;                        // vector for mean vector of item parameters 
  
}

model{
  
  sigma_thetat = 1;
  
  sigma_taut ~ exponential(1);
  sigma_tauc ~ exponential(1);
  sigma_beta ~ exponential(1);
  sigma_alpha ~ exponential(1);
  sigma_a    ~ exponential(1);
  sigma_b  ~ exponential(1);
  
  
  tau_t    ~ normal(0,sigma_taut);
  tau_c    ~ normal(0,sigma_tauc);
  
  theta_t ~ normal(0,1);
  theta_c ~ normal(0,1);
  
  
  mu_beta      ~ normal(beta0,1);

  beta         ~ normal(mu_beta,sigma_beta);

  alpha    ~ inv_gamma(v1,v2);
      
  mu_b     ~ normal(0,1);
  sigma_b  ~ exponential(1);
  b        ~ normal(mu_b,sigma_b);
  
  mu_a     ~ lognormal(0,0.5);
  
  a        ~ lognormal(mu_a,sigma_a);
    
  C ~ beta(1,1);
  T ~ beta(1,1);
      
  for (i in 1:n_obs) {
    
    real p_t = beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]];
    real p_c = beta[ind_item_obs[i]]-tau_c[ind_person_obs[i]];
    
      real lprt1 = log1m(C[ind_item_obs[i]]) + log1m(T[ind_item_obs[i]]) + normal_lpdf(RT[i] | p_t, 1/alpha[ind_item_obs[i]]);  // T = 0, C=0
      real lprt2 = log1m(C[ind_item_obs[i]]) + log(T[ind_item_obs[i]])   + normal_lpdf(RT[i] | p_t, 1/alpha[ind_item_obs[i]]);  // T = 1, C=0
      real lprt3 = log(C[ind_item_obs[i]])   + log1m(T[ind_item_obs[i]]) + normal_lpdf(RT[i] | p_t, 1/alpha[ind_item_obs[i]]);  // T = 0, C=1
      real lprt4 = log(C[ind_item_obs[i]])   + log(T[ind_item_obs[i]])   + normal_lpdf(RT[i] | p_c, 1/alpha[ind_item_obs[i]]);  // T = 1, C=1 
      
      target += log_sum_exp([lprt1, lprt2, lprt3, lprt4]);
  
    real t_t =  inv_logit(a[ind_item_obs[i]]*(theta_t[ind_person_obs[i]] - b[ind_item_obs[i]]));
    real t_c =  inv_logit(a[ind_item_obs[i]]*(theta_c[ind_person_obs[i]] - b[ind_item_obs[i]]));
    
      real lpr1 = log1m(C[ind_item_obs[i]]) + log1m(T[ind_item_obs[i]]) + bernoulli_lpmf(R[i] | p_t, 1/alpha[ind_item_obs[i]]);  // T = 0, C=0
      real lpr2 = log1m(C[ind_item_obs[i]]) + log(T[ind_item_obs[i]])   + bernoulli_lpmf(R[i] | p_t, 1/alpha[ind_item_obs[i]]);  // T = 1, C=0
      real lpr3 = log(C[ind_item_obs[i]])   + log1m(T[ind_item_obs[i]]) + bernoulli_lpmf(R[i] | p_t, 1/alpha[ind_item_obs[i]]);  // T = 0, C=1
      real lpr4 = log(C[ind_item_obs[i]])   + log(T[ind_item_obs[i]])   + bernoulli_lpmf(R[i] | p_c, 1/alpha[ind_item_obs[i]]);  // T = 1, C=1 
      
      target += log_sum_exp([lpr1, lpr2, lpr3, lpr4]);
  }
}

