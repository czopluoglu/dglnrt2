// https://study.sagepub.com/sites/default/files/chapter16.pdf, page 10
// https://discourse.mc-stan.org/t/assigning-bernoulli-prior-to-missing-entries-in-covariate/10219/2
// https://mc-stan.org/docs/2_27/stan-users-guide/vectorizing-mixtures.html
// https://elevanth.org/blog/2018/01/29/algebra-and-missingness/

data{
  int <lower=1> I;                       // number of examinees          
  int <lower=1> J;                       // number of items
  real RT[I,J];                          // RT
  int <lower=0,upper=1> R[I,J];         // item responses
}

parameters {
  
  vector<lower=0,upper=1>[J] pC; // vector of length J for the probability of item compromise status
  
  vector<lower=0,upper=1>[I] pH; // vector of length I for the probability of examinee item peknowledge 
  
  real mu_beta;                 // mean for time intensity parameters
  real<lower=0> sigma_beta;     // sd for time intensity parameters
  
  real mu_alpha;                // mean for log of time discrimination parameters
  real<lower=0> sigma_alpha;    // sd for time discrimination parameters
  
  real<lower=0> sigma_taut;     // sd for tau_t
  real<lower=0> sigma_tauc;     // sd for tau_c
  
  corr_matrix[2] omega_P;       // 2 x 2 correlation matrix for person parameters
  corr_matrix[2] omega_I;       // 2 x 2 correlation matrix for item parameters
  

  ordered[2] person[I];           // an array with length I for person specific latent parameters
                                  // Each array has two elements
                                  // first element is tau_t
                                  // second element is tau_c
                                  // ordered vector makes sure that tau_c > tau_t for every person
                                  // to make sure chains are exploring the same mode and 
                                  // do not go east and west leading multi-modal posteriors
  
  
  vector[2] item[J];           // an array with length J for item specific parameters
                               // each array has two elements
                               // first element is log(alpha)
                               // second element is beta
  
  
  // Parameters for response accuracy portion
  
  real mu_b;                 // mean for item difficulty parameters
  real<lower=0> sigma_b;     // sd for item difficulty parameters
  
  real mu_a;                // mean for log of item discrimination parameters
  real<lower=0> sigma_a;    // sd for item discrimination parameters
  
  corr_matrix[2] omega_P2;       // 2 x 2 correlation matrix for person parameters
  corr_matrix[2] omega_I2;       // 2 x 2 correlation matrix for item parameters
  
  
  ordered[2] person2[I];           // an array with length I for person specific latent parameters
                                   // Each array has two elements
                                   // first element is theta_t
                                   // second element is theta_c
                                   // ordered vector makes sure that theta_c > theta_t for every person
                                   // to make sure chains are exploring the same mode and 
                                   // do not go east and west leading multi-modal posteriors
  
  
  vector[2] item2[J];           // an array with length J for item specific parameters
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
  cov_matrix[2] Sigma_I;                 // covariance matrix for item parameters


  vector[2] mu_P2;                        // vector for mean vector of person parameters 
  vector[2] mu_I2;                        // vector for mean vector of item parameters

  vector[2] scale_P2;                     // vector of standard deviations for person parameters
  vector[2] scale_I2;                     // vector of standard deviations for item parameters
  
  cov_matrix[2] Sigma_P2;                 // covariance matrix for person parameters
  cov_matrix[2] Sigma_I2;                 // covariance matrix for item parameters
  
  mu_P[1] = 0;
  mu_P[2] = 0;
  scale_P[1] = sigma_taut;               
  scale_P[2] = sigma_tauc;
  Sigma_P = quad_form_diag(omega_P, scale_P); 
  
  mu_I[1] = mu_alpha;
  mu_I[2] = mu_beta;
  scale_I[1] = sigma_alpha;               
  scale_I[2] = sigma_beta;
  Sigma_I  = quad_form_diag(omega_I, scale_I); 
  
  
  
  mu_P2[1] = 0;
  mu_P2[2] = 0;
  scale_P2[1] = 1;                        // fix sd of theta_t to 1
  scale_P2[2] = 1;                        // fix sd of theta_c to 1
  Sigma_P2 = quad_form_diag(omega_P2, scale_P2); 
  
  
  mu_I2[1]  = mu_a;
  mu_I2[2]  = mu_b;
  scale_I2[1]  = sigma_a;               
  scale_I2[2]  = sigma_b;                     
  Sigma_I2 = quad_form_diag(omega_I2, scale_I2); 
  
}

model{
  
  sigma_taut  ~ exponential(1);
  sigma_tauc  ~ exponential(1);
  sigma_beta  ~ exponential(1);
  sigma_alpha ~ exponential(1);
  sigma_a     ~ exponential(1);
  sigma_b     ~ exponential(1);
  
  mu_b         ~ normal(0,1);
  mu_beta      ~ normal(4,1);
  mu_a         ~ lognormal(0,0.5);
  mu_alpha     ~ lognormal(0,0.5);

  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
  omega_P   ~ lkj_corr(1);
  omega_I   ~ lkj_corr(1);
  
  omega_P2   ~ lkj_corr(1);
  omega_I2   ~ lkj_corr(1);
  

  person  ~ multi_normal(mu_P,Sigma_P);

  item    ~ multi_normal(mu_I,Sigma_I);
  
  person2  ~ multi_normal(mu_P2,Sigma_P2);

  item2    ~ multi_normal(mu_I2,Sigma_I2);
  

  for (i in 1:I) {
    for(j in 1:J) {
      
      // item[j,1] is log of parameter alpha for item j, therefore use exp(item[j,1]) below
      // item[j,2] is parameter beta for item j
      // person[i,3] is parameter tau_t for person i
      // person[i,4] is parameter tau_c for person i
      
      
      // item2[j,1] is log of parameter a for item j, therefore use exp(item2[j,1]) below
      // item2[j,2] is parameter b for item j
      // person2[i,1] is parameter theta_t
      // person2[i,2] is parameter theta_c
      
      
    real p_t = item[j,2]-person[i,1]; 
    real p_c = item[j,2]-person[i,2];
    
      real lprt1 = log1m(pC[j]) + log1m(pH[i]) + normal_lpdf(RT[i,j] | p_t, 1/exp(item[j,1]));  // T = 0, C=0
      real lprt2 = log1m(pC[j]) + log(pH[i])   + normal_lpdf(RT[i,j] | p_t, 1/exp(item[j,1]));  // T = 1, C=0
      real lprt3 = log(pC[j])   + log1m(pH[i]) + normal_lpdf(RT[i,j] | p_t, 1/exp(item[j,1]));  // T = 0, C=1
      real lprt4 = log(pC[j])   + log(pH[i])   + normal_lpdf(RT[i,j] | p_c, 1/exp(item[j,1]));  // T = 1, C=1 
      
    real t_t =  exp(item2[j,1])*(person2[i,1] - item2[j,2]); 
    real t_c =  exp(item2[j,1])*(person2[i,2] - item2[j,2]); 
  
      real lp1 = log1m(pC[j]) + log1m(pH[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=0
      real lp2 = log1m(pC[j]) + log(pH[i])   + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 1, C=0
      real lp3 = log(pC[j])   + log1m(pH[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=1
      real lp4 = log(pC[j])   + log(pH[i])   + bernoulli_logit_lpmf(R[i,j] | t_c);  // T = 1, C=1 

      target += log_sum_exp([lprt1, lprt2, lprt3, lprt4]);

      target += log_sum_exp([lp1, lp2, lp3, lp4]);
  }
 }

}

