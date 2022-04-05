require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)

# This simulates the null condition using the exact same data structure
# as in the real dataset using the parameters estimated from the dataset


################################################################################
require(MASS)

set.seed(06202021)

################################################################################

N = 1000    # number of examinees
n = 171     # number of items

# Item parameters

  rho_i   <- 0.42
  sd_a    <- 0.08
  sd_b    <- 0.31
  Sigma_i <- diag(c(sd_a,sd_b)) %*% matrix(c(1,0.42,0.42,1),2,2) %*% diag(c(sd_a,sd_b))

  ipar <- mvrnorm(n,mu = c(0.79,4.57),Sigma=Sigma_i)

  
  alpha <- exp(ipar[,1])
  beta  <- ipar[,2]

# Tau_t and tau_c

tau <- mvrnorm(N,
               mu = c(-0.12,0.59),
               Sigma = matrix(c(0.29^2,0,0,0.57^2),2,2))

tau_t <- tau[,1]
tau_c <- tau[,2]

# Assign the examinees with item prekowledge

H <- rep(0,1000) # no examinee item preknowledge

  #H[c(95,135,157,169,174,184,206,207,229,252,287,326,343,347,354,370,442,448,481,517,520,566,568,584,586,603,641,
  #    655,656,660,713,719,750,752,753,832,854,855,944,946,999)] <- 1

# Assign the compromised items

C <- rep(0,171) # no compromised item

  # C[1:50] <- 1

# Generate observed response times according to the model

rt <- matrix(nrow=N,ncol=n)

for(i in 1:N){
  for(j in 1:n){
    
    p_t <- beta[j] - tau_t[i]
    p_c <- beta[j] - tau_c[i]
    
    if(H[i] == 1 & C[j] == 1){
      rt[i,j] = exp(rnorm(1,p_c,1/alpha[j]))
    } else {
      rt[i,j] = exp(rnorm(1,p_t,1/alpha[j]))
    }
    
  }
}

# Insert the missingness by design

  # every examinees gets the first 50 items
  # then gets a randomly selected 15 items from the remaining 121 items

  for(i in 1:N){
    
    rt[i,sample(51:n,106)] <- NA
    
  }



# Convert it to data frame and add group membership and a unique ID

rt       <- as.data.frame(rt)
rt$group <- H
rt$id    <- 1:nrow(rt)

# Check the data

head(rt)

# Reshape it to long format (for plotting purposes)

rt.long <- reshape(data        = rt,
                   idvar       = 'id',
                   varying     = list(colnames(rt)[1:n]),
                   timevar     = "Item",
                   times       = 1:n,
                   v.names      = "RT",
                   direction   = "long")

rt.long <- na.omit(rt.long)

d <- rt.long
################################################################################



data_rt <- list(
  J              = length(unique(d$Item)),
  I              = length(unique(d$id)),
  n_obs          = nrow(d),
  p_loc          = d$id,
  i_loc          = d$Item,
  Y              = log(d$RT)
)


# Compile the model syntax

mod <- cmdstan_model(here('R/ncme22/dglnrt/dglnrt.stan'))

# Fit the model

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup   = 250,
  iter_sampling = 750,
  refresh = 10,
  adapt_delta = 0.99)

# Compile the output files into an rstan object


fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

# Estimation Time 

get_elapsed_time(stanfit)

(sum(get_elapsed_time(stanfit))/4)/3600

# Analyze the parameter estimates


T <- as.numeric(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary[,1])
View(summary(stanfit, pars = c("pH"), probs = c(0.025, 0.975))$summary)
T

T[which(T>0.95)]


View(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary)

C <- as.numeric(summary(stanfit, pars = c("pC"), probs = c(0.025, 0.975))$summary[,1])
C

C[which(C>0.95)]

gr <- c(rep('operational',50),rep('pilot',121))

describeBy(C,gr)

plot(density(C[1:50]),type='l',xlim=c(0,1),ylim=c(0,15))
points(density(C[51:171]),type='l')


summary(stanfit, pars = c("mu_beta"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("sigma_beta"), probs = c(0.025, 0.975))$summary

summary(stanfit, pars = c("sigma_taut"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("sigma_tauc"), probs = c(0.025, 0.975))$summary

summary(stanfit, pars = c("omega_P"), probs = c(0.025, 0.975))$summary
summary(stanfit, pars = c("omega_I"), probs = c(0.025, 0.975))$summary


View(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary)

tau <- matrix(summary(stanfit, pars = c("person"), probs = c(0.025, 0.975))$summary[,1],
              1000,2,byrow=T)


View(summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary)


ipar <- summary(stanfit, pars = c("item"), probs = c(0.025, 0.975))$summary

ipar <- matrix(ipar[,1],ncol=2,byrow=TRUE)
ipar[,1] <- exp(ipar[,1])
  
cor(beta,ipar[,2])
plot(beta,ipar[,2])

cor(alpha,ipar[,1])
plot(alpha,ipar[,1])

describe(ipar)









