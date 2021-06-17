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
  real beta0;                            // prior for time intensity
  real v1;                               // hyperprior parameter for time disc.
  real v2;                               // hyperprior parameter for time disc.
}

parameters {
  vector[J] beta;                        // vector for time intensity parameters
  vector <lower=0> [J]  alpha;           // vector for time discrimination parameters
  vector[I] tau_t;                       // vector for latent speed for honest responses
  vector[I] tau_c;                       // vector for latent speed for fraudulent responses
  real mu_beta;                          // mean for time intensity parameters  
  real<lower=0> sigma_beta;              // sd for time intensity parameters
  real<lower=0> sigma_taut;              // sd for tau_t
  real<lower=0> sigma_tauc;              // sd for tau_c
  
  vector[J] b;                           // vector for item difficulty parameters
  vector[I] theta_t;                     // vector for latent ability for honest responses
  vector[I] theta_c;                     // vector for latent ability for fraudulent responses
  real mu_b;                             // mean for item difficulty parameters  
  real<lower=0> sigma_b;                 // sd for item difficulty parameters
  
  vector<lower=0,upper=1>[J] C;
 
}

transformed parameters {
  
  vector[I] T;

  for (i in 1:n_obs) {
    if((tau_t[ind_person_obs[i]]   < tau_c[ind_person_obs[i]]) &&
       (theta_t[ind_person_obs[i]] < theta_c[ind_person_obs[i]]))
      
      T[ind_person_obs[i]] = 1;
    
    else 
      
      T[ind_person_obs[i]] = 0;
  }
  
}


model{
  
  sigma_taut ~ exponential(1);
  sigma_tauc ~ exponential(1);
  
  tau_t    ~ normal(0,sigma_taut);
  tau_c    ~ normal(0,sigma_tauc);
  
  theta_t ~ normal(0,1);
  theta_c ~ normal(0,1);
  
  
  mu_beta      ~ normal(beta0,1);
  sigma_beta   ~ exponential(1);
      beta     ~ normal(mu_beta,sigma_beta);

      alpha    ~ inv_gamma(v1,v2);
      
  mu_b         ~ normal(0,1);
  sigma_b      ~ exponential(1);
      b        ~ normal(mu_b,sigma_b);
      
  C ~ uniform(0,1);
      
  for (i in 1:n_obs) {
    
    real p_t = beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]];
    real p_c = beta[ind_item_obs[i]]-tau_c[ind_person_obs[i]];
    
    real p = (p_t^(1-T[ind_person_obs[i]]))*
      (((1-C[ind_item_obs[i]])*p_t + 
          (C[ind_item_obs[i]])*p_c)^T[ind_person_obs[i]]);
    
    RT[i] ~ normal(p,1/(alpha[ind_item_obs[i]]));
    
    real t_t =  inv_logit(1.7*(theta_t[ind_person_obs[i]] - b[ind_item_obs[i]]));
    real t_c =  inv_logit(1.7*(theta_c[ind_person_obs[i]] - b[ind_item_obs[i]]));
    
    real t = (t_t^(1-T[ind_person_obs[i]]))*
      (((1-C[ind_item_obs[i]])*t_t + 
          (C[ind_item_obs[i]])*t_c)^T[ind_person_obs[i]]);
          
    R[i] ~ bernoulli(t);
  }
}

