// For tabular data, assumes no missing response
// only item response component

data{
  int <lower=1> I;          // number of examinees          
  int <lower=1> J;          // number of items
  int <lower=0,upper=1> R[I,J];              // Item response matrix
}

parameters {
  real mu_b;                 // mean for item difficulty parameters
  real<lower=0> sigma_b;     // sd for item difficulty parameters
  
  real mu_a;                // mean for log of item discrimination parameters
  real<lower=0> sigma_a;    // sd for item discrimination parameters
  
  corr_matrix[2] omega_P;       // 2 x 2 correlation matrix for person parameters
  corr_matrix[2] omega_I;       // 2 x 2 correlation matrix for item parameters
  
  vector<lower=0,upper=1>[J] pC; // vector of length J for the probability of item compromise status
  
  vector<lower=0,upper=1>[I] pH; // vector of length I for the probability of examinee item peknowledge 
  
  ordered[2] person[I];           // an array with length I for person specific latent parameters
  // Each array has two elements
  // first element is theta_t
  // second element is theta_c
  // ordered vector makes sure that theta_c > theta_t for every person
  // to make sure chains are exploring the same mode and 
  // do not go east and west leading multi-modal posteriors
  
  
  vector[2] item[J];           // an array with length J for item specific parameters
  // each array has two elements
  // first element is log of a
  // second element is log of b
}


transformed parameters{
  
  vector[2] mu_P;                        // vector for mean vector of person parameters 
  vector[2] mu_I;                        // vector for mean vector of item parameters
  
  vector[2] scale_P;                     // vector of standard deviations for person parameters
  vector[2] scale_I;                     // vector of standard deviations for item parameters
  
  cov_matrix[2] Sigma_P;                 // covariance matrix for person parameters
  cov_matrix[2] Sigma_I;                 // covariance matrix for person parameters
  
  mu_P[1] = 0;
  mu_P[2] = 0;
  
  scale_P[1] = 1;               
  scale_P[2] = 1;
  
  Sigma_P = quad_form_diag(omega_P, scale_P); 
  
  mu_I[1] = mu_a;
  mu_I[2] = mu_b;
  
  scale_I[1] = sigma_a;               
  scale_I[2] = sigma_b;
  
  Sigma_I = quad_form_diag(omega_I, scale_I); 
  
}



model{
  
  sigma_b ~ exponential(1);
  sigma_a ~ exponential(1);
  
  mu_b    ~ normal(0,1);
  mu_a    ~ lognormal(0,0.5);
  
  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
  omega_P   ~ lkj_corr(1);
  omega_I   ~ lkj_corr(1);
  
  person  ~ multi_normal(mu_P,Sigma_P);
  
  item    ~ multi_normal(mu_I,Sigma_I);

  for (i in 1:I) {
    for(j in 1:J) {

      // item[j,1] represents log of parameter a of the jth item
      // that's why we use exp(item[j,1]) below 
      // item[j,2] represents parameter b of the jth item
      
      //person[i,1] represents parameter theta_t of the ith person
      //person[i,2] represents parameter theta_c of the ith person

      real t_t =  exp(item[j,1])*(person[i,1] - item[j,2]);
      real t_c =  exp(item[j,1])*(person[i,2] - item[j,2]);

      // log of probability densities for each combination of two discrete parameters
      // (C,T) = {(0,0),(0,1),(1,0),(1,1)}
  
      real lp1 = log1m(pC[j]) + log1m(pH[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=0
      real lp2 = log1m(pC[j]) + log(pH[i])   + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 1, C=0
      real lp3 = log(pC[j])   + log1m(pH[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=1
      real lp4 = log(pC[j])   + log(pH[i])   + bernoulli_logit_lpmf(R[i,j] | t_c);  // T = 1, C=1 
      
      target += log_sum_exp([lp1, lp2, lp3, lp4]);
      
    }
  }
}

