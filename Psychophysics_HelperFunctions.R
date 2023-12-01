### Helper functions for our psychophysics analysis

### Licensing ###
## Third-party code: 
# Contains code modified from https://github.com/stan-dev/example-models/tree/master/Bayesian_Cognitive_Modeling/CaseStudies/PsychophysicalFunctions 
# That code is subject to the following "new BSD license":
# Copyright 2014 martin-smira
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#    - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#    - Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Non-third-party code:
# Copyright 2021-23 Raphael Schween & Ruslan Spartakov - see repo copyright notice.

### Functions related to stan fit ### - require stan package to be loaded
run_stan_fit = function(model, data, myinits, parameters = NA, niter, nchains) {
    samples = stan(model_code=model,
                 data=data, 
                 init=myinits,
                 pars=parameters,
                 iter=niter, 
                 chains=nchains, 
                 #thin=1,
                 warmup=500, 
                 control = list( max_treedepth =15, adapt_delta = .99)
                 # seed = 123  # Default is random seed
  )
  return(samples)
}

# TODO: Generalize parameter extraction function. 
extract_parameters <- function(samples) {
  # Extract the parameters of interest for single participant
  alpha      = extract(samples)$alpha
  beta      = extract(samples)$beta
  
  (alphaMAP   <- map_estimate(alpha))
  (betaMAP    <- map_estimate(beta))
  (alpha_sel <- sample(alpha,20))
  (beta_sel  <- sample(beta,20))
  
  output <- list(alpha, beta, alphaMAP, betaMAP, alpha_sel, beta_sel)
  return(output)
}

extract_parameters_group = function(samples,nsubjs) {
  zmean = list()
  for (i in 1:nsubjs){
    ztemp = matrix(NA,nstim[i],1)
    for (j in 1:nstim[i]){
      ztemp[j,] = summary(samples)$summary[paste("z[",as.character(i),",",as.character(j),"]",sep=""),"mean"]
    }
    zmean[[i]] = ztemp
  }
  spl = extract(samples)
  alpha2    = spl$alphaout
  beta2     = spl$betaout
  alphaMAP2  = c(rep(0,nsubjs))
  betaMAP2  = c(rep(0,nsubjs))
  alpha_sel2 = matrix(NA,20,nsubjs) 
  beta_sel2 = matrix(NA,20,nsubjs) 
  
  # Constructing MAP-estimates and alpha/beta range
  for (i in 1:nsubjs)
  {
    alphaMAP2[i]   = density(alpha2[,i])$x[which(density(alpha2[,i])$y==max(density(alpha2[,i])$y))]
    betaMAP2[i]    = density(beta2[,i])$x[which(density(beta2[,i])$y==max(density(beta2[,i])$y))]
    alpha_sel2[,i] = sample(alpha2[,i],20)
    beta_sel2[,i]  = sample(beta2[,i],20)
  }
  
  return(list(alpha2,beta2,alphaMAP2,betaMAP2,alpha_sel2,beta_sel2,zmean))
}

## Legacy version - likely for group fit without contaminants.
# paramsGroup = function(samples,nsubjs) {
#   # Extracting the necessary parameters
  
#   alpha2    = extract(samples)$alpha
#   beta2     = extract(samples)$beta
#   alphaMAP2  = c(rep(0,nsubjs))
#   betaMAP2  = c(rep(0,nsubjs))
#   alpha_sel2 = matrix(NA,20,nsubjs) 
#   beta_sel2 = matrix(NA,20,nsubjs) 
  
#   # Constructing MAP-estimates and alpha/beta range
#   for (i in 1:nsubjs)
#   {
#     alphaMAP2[i]   = density(alpha2[,i])$x[which(density(alpha2[,i])$y==max(density(alpha2[,i])$y))]
#     betaMAP2[i]    = density(beta2[,i])$x[which(density(beta2[,i])$y==max(density(beta2[,i])$y))]
#     alpha_sel2[,i] = sample(alpha2[,i],20)
#     beta_sel2[,i]  = sample(beta2[,i],20)
#   }
  
#   return(list(alpha2,beta2,alphaMAP2,betaMAP2,alpha_sel2,beta_sel2))
# }

### JND & PSE from alpha & beta ###
# Modified from https://github.com/stan-dev/example-models/tree/master/Bayesian_Cognitive_Modeling/CaseStudies/PsychophysicalFunctions 
# See license note at top. 
# only the MAP estimate; use this to plot psychometric functions
F1 <- function(X, s = 1, alphaMAP, betaMAP, xmean) {
  out <- exp(alphaMAP[s] + betaMAP[s]*(X - xmean[s]))/
    (1+exp(alphaMAP[s] + betaMAP[s]*(X - xmean[s])))
  return(out)
}

F1inv <- function(Y, s=1, alphaMAP, betaMAP) {
  out <- (log(Y/(1-Y))-alphaMAP[s])/betaMAP[s] #this is called the logit
  return(out)
}

# function for all the posterior alpha/beta values; use this to calculate JND 
# posterior
F2 <- function(X, s=1, alpha, beta, xmean) {
  out <- exp(alpha[,s] + beta[,s]*(X - xmean[s]))/
    (1+exp(alpha[,s] + beta[,s]*(X - xmean[s])))
  return(out)
}

F2inv <- function(Y, s=1, alpha, beta) {
  out <- (log(Y/(1-Y))-alpha[,s])/beta[,s]
  return(out)
}

# function for 20 grabbed posterior alpha/beta values; use this to plot 
# overlapping sigmoids to visualize variance
F3 <- function(X, s=1, g, alpha_sel, beta_sel, xmean){
  out <- exp(alpha_sel[g,s] + beta_sel[g,s]*(X - xmean[s]))/
    (1+exp(alpha_sel[g,s] + beta_sel[g,s]*(X - xmean[s])))
  return(out)
}

### Predictions by maximum-likelihood estimate for Gaussian cue fusion ###
predict_bimodal_JND <- function(unimodal_JND_1, unimodal_JND_2) {
  predicted_bimodal_JND <- sqrt( unimodal_JND_1^2*unimodal_JND_2^2 / (unimodal_JND_1^2 + unimodal_JND_2^2) )
  return(predicted_bimodal_JND)
}

predict_bimodal_PSE <- function(unimodal_PSE_1, unimodal_PSE_2, unimodal_JND_1, unimodal_JND_2) {
  weight1 <- unimodal_JND_2^2 / (unimodal_JND_1^2 + unimodal_JND_2^2)
  weight2 <- unimodal_JND_1^2 / (unimodal_JND_1^2 + unimodal_JND_2^2)
  predicted_bimodal_JNE <- weight1*unimodal_PSE_1 + weight2*unimodal_PSE_2
  return(predicted_bimodal_JNE)
}

# ### Transforms ###
# logitTransform = function(y, s, alpha, beta) (psych::logit(y)-alpha[,s])/beta[,s] 
# logisticTransform = function(x, s, alpha, beta, xmean) psych::logistic(alpha[,s]+beta[,s]*(x-xmean[s]))

