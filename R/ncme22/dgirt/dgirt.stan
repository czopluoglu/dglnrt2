// Stan model syntax with modified deterministic gated item response (Rasch version)

data{
    int <lower=1> I;                       // number of examinees          
    int <lower=1> J;                       // number of items
    int <lower=1> n_obs;                   // number of observations (I xJ - missing responses)
    int <lower=1> p_loc[n_obs];            // person indicator   
    int <lower=1> i_loc[n_obs];            // item indicator
    int <lower=0,upper=1> Y[n_obs];       // vector of item responses
}

parameters {
  real mu_b;                  // mean for item difficulty parameters
  real<lower=0> sigma_b;      // sd for item difficulty parameters
  
 // real mu_thetat;             // mean for theta_t
 real<lower=0> sigma_thetat; // sd for theta_t
  
 // real mu_thetac;             // mean for theta_c
  real<lower=0> sigma_thetac; // sd for theta_c
  
  corr_matrix[2] omega_P;    // 2 x 2 correlation matrix for person parameters
  
  vector<lower=0,upper=1>[J] pC; // vector of length J for the probability of item compromise status
  
  vector<lower=0,upper=1>[I] pH; // vector of length I for the probability of examinee item peknowledge 
  
  ordered[2] person[I];          // an array with length I for person specific latent parameters
  // Each array has two elements
  // first element is theta_t
  // second element is theta_c
  // ordered vector assures that tau_c > tau_t for every person
  // to make sure chains are exploring the same mode and 
  // multiple chains do not go east and west leading multi-modal posteriors
 
  
  vector[J] b;    // vector of item difficulty parameters

}


transformed parameters{
  
  vector[2] mu_P;                        // vector for mean vector of person parameters 
  vector[2] scale_P;                     // vector of standard deviations for person parameters
  cov_matrix[2] Sigma_P;                 // covariance matrix for person parameters
  
  //vector[J] b2;                          // vector of item difficulty parameters

  mu_P[1] = 0;//mu_thetat;
  mu_P[2] = 0;//mu_thetac;
  
  scale_P[1] = sigma_thetat;               
  scale_P[2] = sigma_thetac;
  
  Sigma_P = quad_form_diag(omega_P, scale_P); 
  
  //b2 = b-mean(b);
}



model{
  
  //mu_thetat     ~ normal(0,1);
  sigma_thetat  ~ exponential(1);
  
  //mu_thetac     ~ normal(0,1);
  sigma_thetac  ~ exponential(1);
  
  omega_P       ~ lkj_corr(1);
  person        ~ multi_normal(mu_P,Sigma_P);
  
  mu_b      ~ normal(0,1);
  sigma_b   ~ exponential(1);
   b         ~ normal(mu_b,sigma_b);
  
  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
  for (i in 1:n_obs) {
    
      // b[i_loc[i]] represents b-parameter of the (i_loc[i])th item
      
      //person[p_loc[i],1] represents parameter theta_t of the (p_loc[i])th person
      //person[p_loc[i],2] represents parameter theta_c of the (p_loc[i])th person
      
      
      real p_t = person[p_loc[i],1] - b[i_loc[i]];  // non-cheating response
      real p_c = person[p_loc[i],2] - b[i_loc[i]];  // cheating response
      
      // log of probability densities for each combination of two discrete parameters
      // (C,T) = {(0,0),(0,1),(1,0),(1,1)}
      
      real lprt1 = log1m(pC[i_loc[i]]) + log1m(pH[p_loc[i]]) + bernoulli_logit_lpmf(Y[i] | p_t);  // T = 0, C=0
      real lprt2 = log1m(pC[i_loc[i]]) + log(pH[p_loc[i]])   + bernoulli_logit_lpmf(Y[i] | p_t);  // T = 1, C=0
      real lprt3 = log(pC[i_loc[i]])   + log1m(pH[p_loc[i]]) + bernoulli_logit_lpmf(Y[i] | p_t);  // T = 0, C=1
      real lprt4 = log(pC[i_loc[i]])   + log(pH[p_loc[i]])   + bernoulli_logit_lpmf(Y[i] | p_c);  // T = 1, C=1 
      
      target += log_sum_exp([lprt1, lprt2, lprt3, lprt4]);
  }
  
}


