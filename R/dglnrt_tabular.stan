// For tabular data, assumes no missing response
// only response time component
// for Stan forums post

data{
  int <lower=1> I;          // number of examinees          
  int <lower=1> J;          // number of items
  real RT[I,J];             // Response time matrix
}

parameters {
  vector[J] beta;                        // vector for time intensity parameters
  vector <lower=0> [J]  alpha;           // vector for time discrimination parameters
  vector[I] tau_t;                       // vector for latent speed for honest responses
  vector[I] tau_c;                       // vector for latent speed for fraudulent responses
  real<lower=0> sigma_taut;              // sd for tau_t
  real<lower=0> sigma_tauc;              // sd for tau_c
  vector<lower=0, upper=1>[I] T;        // vector for probability of preknowledge
  vector<lower=0, upper=1>[J] C;        // vector for probability of compromise
   
}

model{
  
  sigma_taut ~ exponential(1);
  sigma_tauc ~ exponential(1);
  
  tau_t    ~ normal(0,sigma_taut);
  tau_c    ~ normal(0,sigma_tauc);
  
  beta     ~ normal(4,1);                // numbers are estimated from observed data
  alpha    ~ inv_gamma(100,128);         // numbers are estimated from observed data
  
  T  ~ beta(1,1);
  C  ~ beta(1,1);
  
  for (i in 1:I) {
    for(j in 1:J) {
      
      real p_t = beta[j]-tau_t[i];
      real p_c = beta[j]-tau_c[i];
  
      real lp1 = log1m(C[j]) + log1m(T[i]) + normal_lpdf(RT[i,j] | p_t, 1/alpha[j]);  // T = 0, C=0
      real lp2 = log1m(C[j]) + log(T[i])   + normal_lpdf(RT[i,j] | p_t, 1/alpha[j]);  // T = 1, C=0
      real lp3 = log(C[j])   + log1m(T[i]) + normal_lpdf(RT[i,j] | p_t, 1/alpha[j]);  // T = 0, C=1
      real lp4 = log(C[j])   + log(T[i])   + normal_lpdf(RT[i,j] | p_c, 1/alpha[j]);  // T = 1, C=1 
      
      target += log_sum_exp([lp1, lp2, lp3, lp4]);
      
    }
  }
}

