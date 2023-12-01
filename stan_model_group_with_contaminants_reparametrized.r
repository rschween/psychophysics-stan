  # Our data constrain the model less than those of Lee & Wagenmakers did. 
  # This leads to bad mixing (see trace plots), resulting in low R_hat.
  # Here, we address this by reparametrizing the model. 
  # cf. https://mc-stan.org/docs/stan-users-guide/reparameterizations.html
  # With our data, the result appears less severely impacted by bad mixing 
  # and becomes acceptable at a reasonable number of samples. 

  # Inputs: 
  # data: list of matrices, each size subject * comparison:
  # NAs should have been removed from comparison
  # x: comparison value
  # r: number of times responded "greater" 
  # n: number of times probed
  # rprop: rate; r/n
  # nstim: number of comparisons used by subject (in case they differ)

model_path = 'results_group_reparametrized'
parameters_to_watch = c("alphaout", "betaout", "z")  # Parameters to be monitored
 

model = "
  // Logistic Psychophysical Function with Contaminants
  data { 
    int nsubjs;
    int ncols;
    int nstim[nsubjs];
    int n[nsubjs,ncols];
    int r[nsubjs,ncols];
    int x[nsubjs,ncols];
    vector[nsubjs] xmean; 
  }
  parameters {
    real mua;
    real mub;
    real mup;
    real <lower=0,upper=1000> sigmaa;
    real <lower=0,upper=1000> sigmab;
    real <lower=0,upper=3> sigmap;
    vector[nsubjs] alpha;
    vector[nsubjs] beta;
    vector[nsubjs] probitphi;
    matrix<lower=0,upper=1>[nsubjs,ncols] pi;    
  } 
  transformed parameters {
    vector[2] lp_parts[nsubjs,ncols];
    vector<lower=0,upper=1>[nsubjs] phi;  
    
    for (i in 1:nsubjs)
      phi[i] = Phi(mup + probitphi[i] * sigmap);
    
    for (i in 1:nsubjs) {
      for (j in 1:nstim[i]) {  
        real theta; 
        theta = inv_logit( (mua + alpha[i] * sigmaa) + (mub + beta[i] * sigmab) * (x[i,j] - xmean[i]));
        
        lp_parts[i,j,1] = log1m(phi[i]) + binomial_lpmf(r[i,j]| n[i,j], theta);
        lp_parts[i,j,2] = log(phi[i]) + binomial_lpmf(r[i,j]| n[i,j], pi[i,j]);
      }
    }
  }
  model {
    // Priors
    mua ~ normal(0, inv_sqrt(.01)); 
    mub ~ normal(0, inv_sqrt(.01)); 
    mup ~ normal(0, 1); 
    sigmaa ~ gamma(1.1,0.01);
    sigmab ~ gamma(1.1,0.01);
    sigmap ~ gamma(1.1,3);
    for (i in 1:nsubjs) 
      pi[i] ~ beta(1, 1);  // can be removed
    
    alpha ~ normal(0, 1);
    beta ~ normal(0, 1);
    probitphi ~ normal(0, 1);
    
    for (i in 1:nsubjs)
      for (j in 1:nstim[i])
        increment_log_prob(log_sum_exp(lp_parts[i,j]));
  }
  generated quantities {
    vector[nsubjs] alphaout;
    vector[nsubjs] betaout;
    int<lower=0,upper=1> z[nsubjs,ncols];
    
    for (i in 1:nsubjs) {
      for (j in 1:nstim[i]) {  
        vector[2] prob;
        
        prob = softmax(lp_parts[i,j]);
        z[i,j] = bernoulli_rng(prob[2]);
      }
      alphaout[i] = alpha[i] * sigmaa + mua;
      betaout[i] = beta[i] * sigmab + mub;
    }
  }"
  
 
