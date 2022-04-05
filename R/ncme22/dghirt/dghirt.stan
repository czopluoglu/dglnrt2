// Stan model syntax combining the Modified DG-IRT for response accuracy and
// DG-LNRT for response time

  
data{
    int <lower=1> I;                       // number of examinees          
    int <lower=1> J;                       // number of items
    int <lower=1> n_obs;                   // number of observations (I xJ - missing responses)
    int <lower=1> p_loc[n_obs];            // person indicator   
    int <lower=1> i_loc[n_obs];            // item indicator
    real RT[n_obs];                        // vector of log of responses
    int <lower=0,upper=1> Y[n_obs];       // vector of item responses
}

parameters {
  
// Parameters common for both response time and response accuracy component  
  
  vector<lower=0,upper=1>[J] pC; // vector of length J for the probability of item compromise status
  
  vector<lower=0,upper=1>[I] pH; // vector of length I for the probability of examinee item peknowledge 

// Parameters for the DG-LNRT - response time component  

  real mu_beta;                 // mean for time intensity parameters
  real<lower=0> sigma_beta;     // sd for time intensity parameters
  
  real mu_alpha;                // mean for log of time discrimination parameters
  real<lower=0> sigma_alpha;    // sd for log of time discrimination parameters
  
  real<lower=0> sigma_taut;     // sd for tau_t
  real<lower=0> sigma_tauc;     // sd for tau_c
  
  corr_matrix[2] omega_tau;     // 2 x 2 correlation matrix for person parameters
  corr_matrix[2] omega_item_rt; // 2 x 2 correlation matrix for item parameters
  
  ordered[2] tau[I];            // an array with length I for person specific 
                                // latent speed parameters
                                // Each array has two elements
                                // first element is tau_t
                                // second element is tau_c
                                // ordered vector assures that tau_c > tau_t for every person
                                // to make sure chains are exploring the same mode and 
                                // multiple chains do not go east and west leading multi-modal posteriors
  
  vector[2] item_rt[J];         // an array with length J for item specific parameters
                                // each array has two elements
                                // first element is alpha
                                // second element is beta
                        
// Parameters for the DG-IRT - response accuracy component  

  real mu_b;                    // mean for item difficulty parameters
  real<lower=0> sigma_b;        // sd for item difficulty parameters
  
  real mu_thetat;               // mean for theta_t
  real<lower=0> sigma_thetat;   // sd for theta_t
  
  real mu_thetac;               // mean for theta_c
  real<lower=0> sigma_thetac;   // sd for theta_c
  
  corr_matrix[2] omega_theta;   // 2 x 2 correlation matrix for person parameters
  
  
  ordered[2] theta[I];          // an array with length I for person specific latent parameters
                                // Each array has two elements
                                // first element is theta_t
                                // second element is theta_c
                                // ordered vector assures that theta_c > theta_t for every person
                                // to make sure chains are exploring the same mode and 
                                // multiple chains do not go east and west leading multi-modal posteriors
 
  
  vector[J] b;                  // vector of item difficulty parameters


}


transformed parameters{
  
// Transformed parameters for the DG-LNRT - response time component  

  vector[2] mu_tau;              // vector for mean of latent speed parameters
  vector[2] mu_item_rt;          // vector for mean of response time model item parameters
  
  vector[2] scale_tau;           // vector of standard deviations for latent speed parameters
  vector[2] scale_item_rt;       // vector of standard deviations for response time model item parameters
  
  cov_matrix[2] Sigma_tau;       // covariance matrix for latent speed parameters
  cov_matrix[2] Sigma_item_rt;   // covariance matrix for  response time model item parameters

// Transformed parameters for the DG-IRT - response accuracy component  
    
  vector[2] mu_theta;                // vector of mean for latent trait parameters 
  vector[2] scale_theta;             // vector of standard deviations for latent trait parameters
  cov_matrix[2] Sigma_theta;         // covariance matrix for latent trait parameters
  vector[J] b2;                      // vector of transformed item difficulty parameters (mean-centered)
                                     // We will put a constrain such that the average of b is zero
                                     // This is necessary for model identification of the Rasch model

// DG-LNRT component  

  mu_tau[1] = 0;
  mu_tau[2] = 0;
  scale_tau[1] = sigma_taut;               
  scale_tau[2] = sigma_tauc;
  Sigma_tau = quad_form_diag(omega_tau, scale_tau); 
  
  mu_item_rt[1] = mu_alpha;
  mu_item_rt[2] = mu_beta;
  scale_item_rt[1] = sigma_alpha;               
  scale_item_rt[2] = sigma_beta;
  Sigma_item_rt = quad_form_diag(omega_item_rt, scale_item_rt); 
  
// DG-IRT component
  

  mu_theta[1] = mu_thetat;
  mu_theta[2] = mu_thetac;
  scale_theta[1] = sigma_thetat;               
  scale_theta[2] = sigma_thetac;
  Sigma_theta = quad_form_diag(omega_theta, scale_theta); 
  
  b2 = b-mean(b);
  
  
}



model{

// Sampling the parameters of the DG-LNRT component  

  sigma_taut  ~ exponential(1);
  sigma_tauc  ~ exponential(1);
  sigma_beta  ~ exponential(1);
  sigma_alpha ~ exponential(1);
  
  mu_beta      ~ normal(4,1);
  mu_alpha     ~ lognormal(0,0.5);
  
  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
  omega_tau     ~ lkj_corr(1);
  omega_item_rt ~ lkj_corr(1);
  
  tau     ~  multi_normal(mu_tau,Sigma_tau);
  
  item_rt ~  multi_normal(mu_item_rt,Sigma_item_rt);

// Sampling the parameters of the DG-IRT component  

  mu_thetat     ~ normal(0,1);
  sigma_thetat  ~ exponential(1);
  
  mu_thetac     ~ normal(0,1);
  sigma_thetac  ~ exponential(1);
  
  omega_theta   ~ lkj_corr(1);
  theta         ~ multi_normal(mu_theta,Sigma_theta);
  
  mu_b      ~ normal(0,1);
  sigma_b   ~ exponential(1);
   b         ~ normal(mu_b,sigma_b);

// Sampling the probability of item compromise and probability of item preknowledge

  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
// Joint density of response time and response accuracy  
  
  for (i in 1:n_obs) {
    
  //Response time Component
    
    // item_rt[i_loc[i],1] represents log of parameter alpha of the (i_loc[i])th item
    // that's why we use exp(item[i_loc[i],1]) below 
      // item_rt[i_loc[i],1] represents parameter beta of the (i_loc[i])th item
      
      //tau[p_loc[i],1] represents parameter tau_t of the (p_loc[i])th person
      //tau[p_loc[i],2] represents parameter tau_c of the (p_loc[i])th person
      
      real p_taut = item_rt[i_loc[i],2] - tau[p_loc[i],1];   //expected response time for non-cheating response
      real p_tauc = item_rt[i_loc[i],2] - tau[p_loc[i],2];  //expected response time for cheating response
      
      // log of probability densities for each combination of two discrete parameters
      // (C,T) = {(0,0),(0,1),(1,0),(1,1)}
      
      real lprt1 = log1m(pC[i_loc[i]]) + log1m(pH[p_loc[i]]) + normal_lpdf(RT[i] | p_taut, 1/exp(item_rt[i_loc[i],1]));  // T = 0, C=0
      real lprt2 = log1m(pC[i_loc[i]]) + log(pH[p_loc[i]])   + normal_lpdf(RT[i] | p_taut, 1/exp(item_rt[i_loc[i],1]));  // T = 1, C=0
      real lprt3 = log(pC[i_loc[i]])   + log1m(pH[p_loc[i]]) + normal_lpdf(RT[i] | p_taut, 1/exp(item_rt[i_loc[i],1]));  // T = 0, C=1
      real lprt4 = log(pC[i_loc[i]])   + log(pH[p_loc[i]])   + normal_lpdf(RT[i] | p_tauc, 1/exp(item_rt[i_loc[i],1]));  // T = 1, C=1 
    
  // Response accuracy component
  
      // b[i_loc[i]] represents b-parameter of the (i_loc[i])th item
      
      //person[p_loc[i],1] represents parameter theta_t of the (p_loc[i])th person
      //person[p_loc[i],2] represents parameter theta_c of the (p_loc[i])th person
      
      
      real p_thetat = theta[p_loc[i],1] - b2[i_loc[i]];  // non-cheating response
      real p_thetac = theta[p_loc[i],2] - b2[i_loc[i]];  // cheating response
      
      // log of probability densities for each combination of two discrete parameters
      // (C,T) = {(0,0),(0,1),(1,0),(1,1)}
      
      real lpr1 = log1m(pC[i_loc[i]]) + log1m(pH[p_loc[i]]) + bernoulli_logit_lpmf(Y[i] | p_thetat);  // T = 0, C=0
      real lpr2 = log1m(pC[i_loc[i]]) + log(pH[p_loc[i]])   + bernoulli_logit_lpmf(Y[i] | p_thetat);  // T = 1, C=0
      real lpr3 = log(pC[i_loc[i]])   + log1m(pH[p_loc[i]]) + bernoulli_logit_lpmf(Y[i] | p_thetat);  // T = 0, C=1
      real lpr4 = log(pC[i_loc[i]])   + log(pH[p_loc[i]])   + bernoulli_logit_lpmf(Y[i] | p_thetac);  // T = 1, C=1 
  
  
      target += log_sum_exp([lprt1, lprt2, lprt3, lprt4]) + log_sum_exp([lpr1, lpr2, lprt3, lpr4]);
  }
  
}


