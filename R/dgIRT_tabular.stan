// For tabular data, assumes no missing response
// only item response component

data{
  int <lower=1> I;          // number of examinees          
  int <lower=1> J;          // number of items
  int <lower=0,upper=1> R[I,J];              // Item response matrix
}

parameters {
  vector[J] b;                           // vector for item difficulty parameters
  vector<lower = 0>[J] a;                // vector for item discrimination parameters
  vector[I] theta_t;                     // vector for latent ability for honest responses
  vector[I] theta_c;                     // vector for latent ability for fraudulent responses
  real mu_b;                             // mean for item difficulty parameters  
  real<lower=0> sigma_b;                 // sd for item difficulty parameters
  real mu_a;                             // mean for item discrimination parameters  
  real<lower=0> sigma_a;                 // sd for item discrimination parameters
  vector<lower=0, upper=1>[I] T;        // vector for probability of preknowledge
  vector<lower=0, upper=1>[J] C;        // vector for probability of compromise
   
}

model{
  
  theta_t    ~ normal(0,1);             // prior for true theta
  theta_c    ~ normal(0,1);             // prior for cheating theta
  
  mu_b     ~ normal(0,1);               // hyperprior for mean of item difficulty
  sigma_b  ~ exponential(1);            // hyperprior for sd of item difficulty
  b        ~ normal(mu_b,sigma_b);      // prior for item difficulty
   
  mu_a     ~ lognormal(0,0.5);          // hyperprior for mean of item discrimination
  sigma_a  ~ exponential(1);            // hyperprior for sd of item discrimination
  a        ~ lognormal(mu_a,sigma_a);   // prior for item discrimination
  
  T  ~ beta(1,1);                       // prior for probability of item preknowledge
  C  ~ beta(1,1);                       // prior for probability of item compromise
  
  for (i in 1:I) {
    for(j in 1:J) {
      
      real t_t =  a[j]*(theta_t[i] - b[j]);
      real t_c =  a[j]*(theta_c[i] - b[j]);
  
      real lp1 = log1m(C[j]) + log1m(T[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=0
      real lp2 = log1m(C[j]) + log(T[i])   + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 1, C=0
      real lp3 = log(C[j])   + log1m(T[i]) + bernoulli_logit_lpmf(R[i,j] | t_t);  // T = 0, C=1
      real lp4 = log(C[j])   + log(T[i])   + bernoulli_logit_lpmf(R[i,j] | t_c);  // T = 1, C=1 
      
      target += log_sum_exp([lp1, lp2, lp3, lp4]);
      
    }
  }
}

