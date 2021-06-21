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
  real mu_beta;                 // mean for time intensity parameters
  real<lower=0> sigma_beta;              // sd for time intensity parameters
  
  real mu_alpha;                // mean for log of time discrimination parameters
  real<lower=0> sigma_alpha;             // sd for log of time discrimination parameters

  real mu_b;                    // mean for item difficulty parameters
  real<lower=0> sigma_b;                 // sd for item difficulty parameters
  
  real mu_a;                    // mean for log of item discrimination parameters
  real<lower=0> sigma_a;                 // sd for log of item discrimination parameters

  real<lower=0> sigma_taut;              // sd for tau_t

  real<lower=0> sigma_tauc;              // sd for tau_c
  
  corr_matrix[4] omega_P;                // correlation matrix for person parameters
	corr_matrix[4] omega_I;                // correlation matrix for item parameters

  vector<lower=0,upper=1>[J] pC;          // vector for the probability of 
                                         // item compromise status

  vector<lower=0,upper=1>[I] pH;          // vector for the probability of examinee 
                                         // item peknowledge 

  vector[4] person[I];         // person specific latent parameters
                               // first element is theta_t
                               // second element is theta_c
                               // third element tau_t
                               // fourth element tau_c

	vector[4] item[J];          // item specific parameters
	                            // first element is a
	                            // second element is b
	                            // third element ic alpha
	                            // fourth element is beta
 
}

transformed parameters{
  
  vector[4] mu_P;                        // vector for mean vector of person parameters 
  vector[4] mu_I;                        // vector for mean vector of item parameters
  
  vector[4] scale_P;                     // vector of standard deviations for person parameters
  vector[4] scale_I;                     // vector of standard deviations for item parameters
  
  cov_matrix[4] Sigma_P;                 // covariance matrix for person parameters
  cov_matrix[4] Sigma_I;                 // covariance matrix for person parameters
  
  mu_P[1] = 0;
  mu_P[2] = 0;
  mu_P[3] = 0;
  mu_P[4] = 0;
  
  scale_P[1] = 1;                        // fix sd of theta_t to 1
  scale_P[2] = 1;                        // fix sd of theta_c to 1
  scale_P[3] = sigma_taut;               
  scale_P[4] = sigma_tauc;
  
  Sigma_P = quad_form_diag(omega_P, scale_P); 
  
  mu_I[1] = mu_a;
  mu_I[2] = mu_b;
  mu_I[3] = mu_alpha;
  mu_I[4] = mu_beta;
  
  scale_I[1] = sigma_a;               
  scale_I[2] = sigma_b;                     
  scale_I[3] = sigma_alpha;               
  scale_I[4] = sigma_beta;
  
  Sigma_I = quad_form_diag(omega_I, scale_I); 
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
  mu_a         ~ normal(0,1);
  mu_alpha     ~ normal(0,1);

  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
  omega_P   ~ lkj_corr(1);
  omega_I   ~ lkj_corr(1);

  person  ~ multi_normal(mu_P,Sigma_P);

  item    ~ multi_normal(mu_I,Sigma_I);


  for (i in 1:I) {
    for(j in 1:J) {
      
      // item[j,1] is log of parameter a, therefore use exp(item[j,1]) below
      // item[j,2] is parameter b
      // item[j,3] is log of parameter alpha, therefore use exp(item[j,3]) below
      // item[j,4] is parameter beta
      
      //person[i,1] is parameter theta_t
      //person[i,2] is parameter theta_c
      //person[i,3] is parameter tau_t
      //person[i,4] is parameter tau_c
      
      
    real p_t = item[j,4]-person[i,3]; 
    real p_c = item[j,4]-person[i,4];
    
      real lprt1 = log1m(pC[j]) + log1m(pH[i]) + normal_lpdf(RT[i,j] | p_t, 1/exp(item[j,3]));  // T = 0, C=0
      real lprt2 = log1m(pC[j]) + log(pH[i])   + normal_lpdf(RT[i,j] | p_t, 1/exp(item[j,3]));  // T = 1, C=0
      real lprt3 = log(pC[j])   + log1m(pH[i]) + normal_lpdf(RT[i,j] | p_t, 1/exp(item[j,3]));  // T = 0, C=1
      real lprt4 = log(pC[j])   + log(pH[i])   + normal_lpdf(RT[i,j] | p_c, 1/exp(item[j,3]));  // T = 1, C=1 
      
    real t_t =  exp(item[j,1])*(person[i,1] - item[j,2]); 
    real t_c =  exp(item[j,1])*(person[i,2] - item[j,2]); 
  
      real lp1 = log1m(pC[j]) + log1m(pH[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=0
      real lp2 = log1m(pC[j]) + log(pH[i])   + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 1, C=0
      real lp3 = log(pC[j])   + log1m(pH[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=1
      real lp4 = log(pC[j])   + log(pH[i])   + bernoulli_logit_lpmf(R[i,j] | t_c);  // T = 1, C=1 

      target += log_sum_exp([lprt1, lprt2, lprt3, lprt4]);

      target += log_sum_exp([lp1, lp2, lp3, lp4]);
  }
 }

}

