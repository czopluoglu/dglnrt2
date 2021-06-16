// For tabular data, assumes no missing response
// only response time component
// for Stan forums post

data{
  int <lower=1> I;          // number of examinees          
  int <lower=1> J;          // number of items
  real RT[I,J];             // Response time matrix
  int C[J];                 // vector for item compromise status (0s and 1s)
}

parameters {
  vector[J] beta;                        // vector for time intensity parameters
  vector <lower=0> [J]  alpha;           // vector for time discrimination parameters
  vector[I] tau_t;                       // vector for latent speed for honest responses
  vector[I] tau_c;                       // vector for latent speed for fraudulent responses
  real<lower=0> sigma_taut;              // sd for tau_t
  real<lower=0> sigma_tauc;              // sd for tau_c
}

transformed parameters {
  
  vector[I] T;

  for (i in 1:I) {
    for(j in 1:J) {
      
      if(tau_t[i] < tau_c[i])
        
        T[i] = 1;
      
      else 
        
        T[i] = 0;
    }
  }
  
}


model{
  
  sigma_taut ~ exponential(1);
  sigma_tauc ~ exponential(1);
  
  tau_t    ~ normal(0,sigma_taut);
  tau_c    ~ normal(0,sigma_tauc);
  
  beta     ~ normal(4,1);
  alpha    ~ inv_gamma(100,128);
      
  for (i in 1:I) {
    for(j in 1:J) {
      
      real p_t = beta[j]-tau_t[i];
      real p_c = beta[j]-tau_c[i];
      real p   = (p_t^(1-T[i]))*(((1-C[j])*p_t + C[j]*p_c)^T[i]);
    
       RT[i,j] ~ normal(p,1/alpha[j]);
    }
  }
}

