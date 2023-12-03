### Licensing and copyright: see repository LICENSE note. ###

# Adapted from PsychophysicalFunction2_Stan.R in https://github.com/stan-dev/example-models/tree/master/Bayesian_Cognitive_Modeling/CaseStudies/PsychophysicalFunctions
# Figure 12.5 in Lee & Wagenmakers 2014.
# Inputs: size subjec * comparison:
# NAs should have been removed from comparison
# x: comparison value
# r: number of times responded "greater" 
# n: number of times probed
# rprop: rate; r/n
# nstim: number of comparisons used by subject (in case they differ)

model_path = 'results_group'
parameters_to_watch = c("alpha", "beta", "z")  # Parameters to be monitored

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
    real<lower=0,upper=1000> sigmaa;
    real<lower=0,upper=1000> sigmab;
    real<lower=0,upper=3> sigmap;
    vector[nsubjs] alpha;
    vector[nsubjs] beta;
    vector[nsubjs] probitphi;
    matrix<lower=0,upper=1>[nsubjs,ncols] pi;    
  } 
  transformed parameters {
    vector[2] lp_parts[nsubjs,ncols];
    vector<lower=0,upper=1>[nsubjs] phi;  
    
    for (i in 1:nsubjs)
      phi[i] = Phi(probitphi[i]);
    
    for (i in 1:nsubjs) {
      for (j in 1:nstim[i]) {  
        real theta; 
        theta = inv_logit(alpha[i] + beta[i] * (x[i,j] - xmean[i]));
        
        lp_parts[i,j,1] = log1m(phi[i]) + binomial_lpmf(r[i,j]| n[i,j], theta);
        lp_parts[i,j,2] = log(phi[i]) + binomial_lpmf(r[i,j]| n[i,j], pi[i,j]);
      }
    }
  }
  model {
    // Priors
    mua ~ normal(0, inv_sqrt(.001));
    mub ~ normal(0, inv_sqrt(.001));
    mup ~ normal(0, 1); 
    for (i in 1:nsubjs) 
      pi[i] ~ beta(1, 1);  // can be removed
    
    alpha ~ normal(mua, sigmaa);
    beta ~ normal(mub, sigmab);
    probitphi ~ normal(mup, sigmap);
    
    for (i in 1:nsubjs)
      for (j in 1:nstim[i])
        increment_log_prob(log_sum_exp(lp_parts[i,j]));
  }
  generated quantities {
    int<lower=0,upper=1> z[nsubjs,ncols];
    for (i in 1:nsubjs) {
      for (j in 1:nstim[i]) {  
        vector[2] prob;
        
        prob = softmax(lp_parts[i,j]);
        z[i,j] = bernoulli_rng(prob[2]);
      }
    }
  }"
  

