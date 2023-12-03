### Licensing and copyright: see repository LICENSE note. ###

model_path = paste0('results_participant_',subj)

samplePosterior <- function(outcome,stan_iterations = 1000) {
  r <- outcome[,2] # number of times comparison judged greater
  n <- r+outcome[,1] # number of times asked by comparison size
  
  nstim <- length(n) # number of comparison stimuli used
  x <- as.numeric(names(r)) # value of comparison stimuli used
  xmean <- mean(x) # mean value comparison
  rprop <- r/n # proportions (for plotting)
  
  ## define the model. This is as in the lecture, but we have 
  # removed the layer summarizing across participants 
  # i.e. we only analyze one participant at a time, for now. 
  model <- "
    // Logistic Psychophysical Function
    data { 
      int nstim;
      int n[nstim];
      int r[nstim];
      real x[nstim];
      real xmean; 
    }
    parameters {
      real mua;
      real mub;
      real<lower=0.01,upper=1000> sigmaa;
      real<lower=0.01,upper=1000> sigmab;
      real alpha;
      real beta;
    } 
    model {
      // Priors
      sigmaa ~ uniform(0.01,1000);
      sigmab ~ uniform(0.01,1000);
      
      mua ~ normal(0, inv_sqrt(.001));
      mub ~ normal(0, inv_sqrt(.001));
      
      alpha ~ normal(mua, sigmaa);
      beta ~ normal(mub, sigmab);
      
      for (i in 1:nstim) {
          real theta; 
          theta <- inv_logit(alpha + beta * (x[i] - xmean));
          r[i] ~ binomial(n[i], theta);
      }
    }"
  
  # Pack data to be passed on to Stan in a list (because that is how Stan wants it)
  data = list(x=x, xmean=xmean, n=n, r=r, nstim=nstim) 
  
  # Set initial values for parameters
  myinits <- list( 
    list(alpha=rep(0), beta=rep(0), 
         mua=0, mub=0, sigmaa=1, sigmab=1),
    list(alpha=rep(0), beta=rep(0), 
         mua=0, mub=0, sigmaa=1, sigmab=1))
  
  parameters = c("alpha", "beta")  # Parameters to be monitored
  
  # Run the stan model. This will take a while and return warnings. It is successful
  # if you have a "samples" variable in your environment, afterwards. 
  # The following command calls Stan with specific options.
  # For a detailed description type "?rstan".
  samples = stan(model_code=model,
                 data=data, 
                 init=myinits,
                 pars=parameters,
                 iter=stan_iterations, 
                 chains=2, 
                 cores =2
                 #thin=1,
                 # warmup = 100,  # Stands for burn-in; Default = iter/2
                 # seed = 123  # Setting seed; Default is random seed
  )
  
  # Collect relevant outputs in a list and return them. 
  out <- list(samples, xmean, r,n, nstim, x, rprop)
  return(out)
}