### Code to analyze the psychophysical duration assessment experiment from ExPra ###
# The psychophysical experiment was based on Hartcher-O'Brien et al. 2014 https://doi.org/10.1371/journal.pone.0089339 
# It presented visual, auditory and combined audivisual cues of varying durations,
# where audiovisual cues could contain a mismatch between auditory and visual duration.
# Participants judged which of 2 cues on each presentation was perceived longer, 
# in a 2-alternative-forced choice task. 
# The hypothesis was that audiovisual cues would show Bayesian cue fusion, as 
# indicated by a) a mean perceived duration between the duration of the auditory and
# visual cue, as predicted by their respective uncertainties, b) an overall decraesed
# uncertainty of bimodal relative to unimodal cue presentation. 
# The following code:
# 1. preprocesses the average data for each participant
# 2. fits psychometric functions to the average data of each condition by a hierarchical 
#    Bayesian model, based on M. Lee & J.-E. Wagenmakers' 2014 Book: Bayesian Cognitive
#    Modelling, Chapter 12. Parameters of the psychometric functions are inferred 
#    by Hamiltonian Monte Carlo Sampling in Stan
# 3. predicts the theoretical maximum likelihood estimate of optimal cue fusion
# 4. plots the results
# Note: A desirable extension would be to implement cue fusion as an additional 
# layer of the hierarchical model in stan and doing fully bayesian analysis rather
# than MLE estimation. 

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

### Prepare environment ###
Sys.setenv(LANG = 'en')

rm(list=ls()) # Clear workspace.

# load required packages
library(Rlab)
library(jsonlite)
library(dplyr)
library(rstan)
# library(gdata)
# library(bayestestR)
library(stringi)
library(sjmisc)

rstan_options(auto_write = TRUE)
options(mc.cores = 2) # parallel::detectCores())# for RStan

# Set the working directory
# setwd("~/teaching/ExPra/Bayesian Integraton/DataAnalysis")
plot_to_pdf = TRUE 
data_path = 'data_example'
# Load functions from external scripts
# CAVE: The model sets the result path. Overwrites without prompting.
source('stan_model_group_with_contaminants_reparametrized.R')
result_path = paste(model_path,data_path, sep = '_')
dir.create(result_path) 
source('Psychophysics_HelperFunctions.R')
source('Psychophysics_ExtractData.R')
source('Psychophysics_PlotFunctions.R')

set.seed(1801) # this is to make "randomness" in the script reproducible
conditions = data.frame( names = c('vis',rep('aud',4)),
                        modality = c('vis','aud01','aud06','av01','av06'),
                        noiseLevel = c('None','01','06','01','06'),  # '01' and '06' indicate 10 and 60% background noise. 
                        type = c('visual','auditory','auditory','audiovisual','audiovisual')
                        )
experiment_standards <- c('450','550') 
# Note: In our design, the 'standards" in bimodal mean the two different conflicts and we use the auditory standard. 
# By  this  convention, bimodal450 would mean auditory: 450, visual 550 and vice versa.
# However, we use the true visual standard as identifier for unimodal visual trials!
experiment_comparison_values = c(100, 300, 400, 450, 550, 600, 700, 900)
experiment_n_modalities = length(conditions$modality)
experiment_n_standards = length(experiment_standards)
experiment_n_comparisons = length(experiment_comparison_values)

parameters_to_estimate <- c('PSE','JND')
n_parameters_to_estimate <- length(parameters_to_estimate)
point_estimates = c('_lower','_map','_upper')
n_point_estimates = length(point_estimates)

# names for all combinations of conditions and values
condition_names <- paste0(rep(conditions$modality,each = n_parameters_to_estimate*experiment_n_standards*n_point_estimates),rep(experiment_standards,each = n_parameters_to_estimate*n_point_estimates, experiment_n_modalities),rep(parameters_to_estimate, each = n_point_estimates), point_estimates)


### Prepare data
# 1) Check if aggregated data exist, based on 'prop.csv'
# Aggregated data are 'responses.csv', 'nresp.csv', 'comp.csv', 
# 'prop.csv' and 'wflabel.csv' according to the format used in 
# examples to Lee & Wagenmakers Ch. 12.:
# https://github.com/stan-dev/example-models/tree/master/Bayesian_Cognitive_Modeling/CaseStudies/PsychophysicalFunctions 
# Select file to be processed
if (!file.exists( file.path(data_path,'prop.csv'))) {
  # If no aggregated data are found, assume that 
  # data are provided as one .json file per experimental participant 
  # in custom/jsPsych format that contains the responses on each individual 
  # trial, including preparation, instructions, etc.
  # Extract the main experiment data and aggregate. 
  writeLines(paste('Could not find prop.csv in directory', getwd(), 
  '. \n  Will check for individual participant data and aggregate.'))
  n_practice_trials = 94 # 76 in Summer2021 version.
  outcome_dfs_list = aggregate_data(data_path, experiment_n_modalities, experiment_n_standards,
                                    experiment_n_comparisons, conditions, experiment_standards, 
                                    experiment_comparison_values, n_practice_trials, result_path)
}

### Fit stan model (loaded at top) to aggregated data
responses = read.csv(file.path(data_path,'responses.csv'))
nresp = read.csv(file.path(data_path,'nresp.csv'))
comp = read.csv(file.path(data_path,'comp.csv'))
prop = read.csv(file.path(data_path,'prop.csv'))

n_subjects = nrow(prop)

# initialize data structures to collect our processed values
output_values <- setNames( data.frame( matrix(
  NA, ncol = experiment_n_modalities*experiment_n_standards*n_parameters_to_estimate*n_point_estimates, nrow = n_subjects) ), condition_names)
collect_outputs = list()

# fit model to data - separately by modalities*standards
for( modality in conditions$modality ) { 
  for( standard in experiment_standards) { 
    # Extract data for current modality*standard combination.
    data_columns = grepl( paste0( modality, standard), names( responses))
    stopifnot( sum(data_columns) == experiment_n_comparisons)
    
    x <- as.matrix( comp[, data_columns]) # Note: if data contain NAs, transform to -99 and transform back, later.
    n <- as.matrix( nresp[ , data_columns])
    r <- as.matrix( responses[, data_columns])
    rprop <- as.matrix( prop[ , data_columns])
    
    # rows & columns appear switched - no idea why. Maybe a matrix thing? 
    xmean = apply(x, MARGIN =  1, mean)
    
    nstim <- rowSums(!is.na(x)) 
    ncols <- max(nstim)
    n_subjects <- dim( x)[1]
    nchains = 2
    niter = 2000
    
    # to be passed on to Stan
    data <- list(x=x, xmean=xmean, n=n, r=r, nsubjs=n_subjects, nstim=nstim, ncols = ncols) 
    
    # myinits <- rep( list( 
    #   list(alpha=rep(0, nsubjs), beta=rep(.05, nsubjs), probitphi=rep(-2, nsubjs),
    #        mua=0, mub=0, mup=-2, sigmaa=1, sigmab=1, sigmap=1, z = 1,  
    #        pi=matrix(runif(nsubjs * ncols), nsubjs, ncols) )), nchains)
    
    myinits = list()
    for (c in 1:nchains) {
      myinits[[c]] <- list(
            alpha = rnorm(n_subjects,0,1), probitphi = rnorm(n_subjects,0,1),
            beta=rep(.05, n_subjects), 
            mua=abs(rnorm(1, 0, 1)),  mup = rnorm(1, 0, .1),
            mub=0, 
            sigmaa=1,sigmab=1,
            sigmap=min(rgamma(1,1.1,3),3),
            z = matrix(rbern(n_subjects * ncols,.5), n_subjects, ncols  ),
            pi=matrix(runif(n_subjects * ncols), n_subjects, ncols) )      
    }

    # Note: to avoid recomputing each time, we save the result and load it automatically.
    # To recompute existing samples, delete corresponding files in folder.
    parameters_to_watch = NA
    samplefile_name = paste('StanSamples',modality,standard,'.RStan',sep = '_')
    if (file.exists( file.path(result_path,samplefile_name))) {
      load( file.path(result_path,samplefile_name))
      print(paste('Found file',samplefile_name,'in', result_path,'. Loading...'))
    } else {
      samples <- run_stan_fit(model, data, myinits, parameters_to_watch, niter, nchains) # using function defined above
      save(samples, data, rprop, file  = file.path(result_path,samplefile_name))  
    }
    
    # ## Diagnose the stan fit - exemplified for mus, should also check other parameters. 
    # # Note: fixes not to be applied in order (in fact, rather from bottom to top, if anything).
    # print(samples, pars = c('mua','mub','mup')) 
    # # Rhat should not be smaller than 1.0, n_eff should not be too low. 
    # # Possible fixes: longer warmup, higher adapt_delt & max_treedepth, more samples.
    # traceplot(samples, pars = c('mua','mub','mup'), nrow = 3) 
    # # Should resemble "fat hairy caterpillar" as opposed to "wiggly snake". 
    # # Different chains should cover same region (but not be identical!)
    # # Possible fix for snake is reparametrization. 
    # pairs(samples, pars = c('mua','mub','mup')) 
    # # Posterior histograms and bivariate correlations. 
    # # Helps identify sampling issues in specific regions.
    # # Cf. https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html
    # # Further, strong correlations may point to identifiability issues between 2 parameters 
    # # Possible fix would be defining stronger priors.

    stopifnot(xmean == data$xmean,r == data$r, n == data$n, x == data$x, ncols == data$ncols) 
    
    # TODO: Automate dealing with parameter names.
    # inferred_parameters <- extract_parameters(samples)
    inferred_parameters <- extract_parameters_group(samples, n_subjects)
    
    # assign to variables named as in external code. 
    alpha2 <- inferred_parameters[[1]]
    beta2 <- inferred_parameters[[2]]
    alphaMAP2 <- inferred_parameters[[3]]
    betaMAP2 <- inferred_parameters[[4]]
    alpha_sel2 <- inferred_parameters[[5]]
    beta_sel2 <- inferred_parameters[[6]]
    zmean <- inferred_parameters[[7]]
    
    JND 	  <- F2inv(0.84,c(1:n_subjects),alpha2,beta2)-F2inv(0.5,c(1:n_subjects),alpha2,beta2)
    JNDmap <- F1inv(0.84,c(1:n_subjects),alphaMAP2, betaMAP2)-F1inv(0.5,c(1:n_subjects),alphaMAP2, betaMAP2)								  				             
    PSE 	  <- F2inv(0.5,c(1:n_subjects),alpha2,beta2)+xmean
    PSEmap <- F1inv(0.5,c(1:n_subjects),alphaMAP2,betaMAP2)+xmean
    
    collect_outputs[[paste(modality,standard,sep='')]] = list( JND = JND, JNDmap = JNDmap, PSE = PSE, PSEmap = PSEmap,
                                                       alphaMAP = alphaMAP2, betaMAP = betaMAP2,
                                                       alpha_sel = alpha_sel2, beta_sel = beta_sel2,
                                                       rprop = rprop, x = x, xmean = xmean, zmean = zmean,
                                                       nstim = nstim)
    
    # assign output values to relevant columns in our output data frame (initialized before loop)
    output_values[paste(modality,standard,'JND','_map',sep = '')] <- JNDmap
    output_values[paste(modality,standard,'PSE','_map',sep = '')] <- PSEmap
    JNDq <- apply(JND, MARGIN = 2, quantile,probs = c(.05, .95))
    PSEq <- apply(PSE, MARGIN = 2, quantile, probs = c(.05, .95))
    output_values[paste(modality,standard,'JND','_lower',sep = '')] <- JNDq[1,]
    output_values[paste(modality,standard,'PSE','_lower',sep = '')] <- PSEq[1,]
    output_values[paste(modality,standard,'JND','_upper',sep = '')] <- JNDq[2,]
    output_values[paste(modality,standard,'PSE','_upper',sep = '')] <- PSEq[2,]
  }
}

### Process results ###
# At this point, the psychometric functions have been fit. 
# What follows is specific to the Bayesian cue fusion experiment,
# and largely follows Hartcher-O'Brien et al. 2014. 
# TODO: Add more commentary. 
weber_fractions = output_values[, grepl('JND',names(output_values),fixed = TRUE)] / cbind( 
  matrix(rep(matrix(c(450,450,450,550,550,550), nrow = n_subjects, ncol = 6, byrow = TRUE),3), nrow = n_subjects),
  matrix(rep(500, 12*n_subjects), nrow = n_subjects))

# Plot fitted curves
output_path_stump = file.path(result_path,paste0('plot_unimodal_Subj', sep=''))
plot_call = call('plot_group_contamin',substitute(collect_outputs),substitute(subj), experiment_standards, conditions$modality)
plot_participants_indexed_as_subj(1:n_subjects, plot_call, plot_to_pdf, output_path_stump)

## Predictions
# We have estimated the PSE and JND for the unimodal values corresponding to the
# standard (500ms) +/- 0.5*delta (50ms). Therefore, we can predict the JND/PSE
# for each biomdal delta from the actual PSEs and JNDs of the unimodal stimuli used.
other_standards = c("550", "450")
noise_levels = c("01","06")
for (n in noise_levels) {
  for (s in 1:experiment_n_standards) {
    for (p in point_estimates) {
      # predict audiovisual
      output_values[paste('pauv',n,experiment_standards[s],'JND',p,sep = '')] <- predict_bimodal_JND( 
        output_values[paste0('aud',n,experiment_standards[s],'JND',p,sep = '')], 
        output_values[paste0('vis',other_standards[s],'JND',p,sep = '')])
    }
    for (p in point_estimates) {
      # predict audiovisual / note: the time in predicted audivisual belongs to audio. 
      output_values[paste('pauv',n,experiment_standards[s],'PSE',p,sep = '')] <- predict_bimodal_PSE( 
        output_values[paste0('aud',n,experiment_standards[s],'PSE',p,sep = '')], 
        output_values[paste0('vis',other_standards[s],'PSE',p,sep = '')],
        output_values[paste0('aud',n,experiment_standards[s],'JND',p,sep = '')], 
        output_values[paste0('vis',other_standards[s],'JND',p,sep = '')])
    }
  }
}

# Plot Gaussian for demonstrating cue fusion
plot_call = call('plot_gaussian',substitute(output_values[subj,]),substitute(subj), experiment_standards)
output_path_stump = file.path(result_path,paste0('plot_bimodal_Gaussians_Subj', sep=''))
plot_participants_indexed_as_subj(1:n_subjects, plot_call, plot_to_pdf, output_path_stump)

# Plots summarizing parameters. You can select 'JND' or 'PSE' as value.
# Individual participants. MAP parameters + Credible interval
val = 'PSE' # which value to plot
if (val=='JND') {ylim = c(0,500)} else { ylim = c(200,800)}
plot_call = call('plot_predictions',substitute(output_values[subj,]),substitute(subj), value = val, ylim = ylim, ylab = val, experiment_standards)
output_path_stump = file.path(result_path, paste0( 'plot_bimodal_parameters_',val,'_','Subj', sep=''))
plot_participants_indexed_as_subj(1:n_subjects, plot_call, plot_to_pdf, output_path_stump)
 
# Mean across participants. MAP paramters +/- SD of MAP parameters. 
plot_call = call('plot_predictions_mean',output_values, value = 'JND', ylim = c(0,500), ylab = 'JND', experiment_standards)
plot_to_output(plot_call, plot_to_pdf, file.path(result_path, 'plot_bimodal_mean_JNDs.pdf'))

# indexing sets for values of interest
pse_map_set = which(grepl('PSE_map',names(output_values), fixed =T))
jnd_map_set = which(grepl('JND_map',names(output_values), fixed =T))
jnd_lo_set = which(grepl('JND_lower',names(output_values), fixed =T))
jnd_hi_set = which(grepl('JND_upper',names(output_values), fixed =T))
pauv_set = which(grepl('pauv',names(output_values), fixed =T))
av_set = which(grepl('av',names(output_values), fixed =T))
vis_set = which(grepl('vis',names(output_values), fixed =T))
aud_set = which(grepl('aud',names(output_values), fixed =T))

# identify subjects to exclude
mean_unimodal_JNDs = apply( output_values[, intersect(jnd_map_set, union(vis_set,aud_set))], MARGIN = 2, mean)
sd_unimodal_JNDs = apply( output_values[, intersect(jnd_map_set, union(vis_set,aud_set))], MARGIN = 2, sd)
hist(output_values[, intersect(jnd_map_set, union(vis_set,aud_set))][,1])
mean_unimodal_JNDs + 3*sd_unimodal_JNDs

export = output_values[, intersect(jnd_map_set, union(pauv_set, av_set))]

export$subject = factor(paste0('Sub',1:dim(export)[1]))
export_long = reshape(export, direction = 'long', varying = 1:(length(export)-1), v.names = 'outcome')#,sep = '')#
export_long$condition = factor(rep(c('1','2'), each = n_subjects*4)) # 1: av, 2: pauv
export_long$noiseLevel = rep(c('01','06'), each = n_subjects*2,2)
export_long$auditoryStandard = rep(c('450','550'), each = n_subjects, 4)

write.csv(export, file.path(result_path,'JNDdata.csv'))
write.csv(export_long, file.path(result_path,'JNDdata_long.csv'))

# export tables
# write.csv(kennWerte[, intersect(jndmapSet, union(visSet,audSet))],'unimodal_JNDs.csv')
# write.csv(kennWerte, 'all_values.csv')

### Empirical Weights - cf. Hartcher-O'Brien et al. 2014
wRawVis_01neg100 = output_values[, "aud01450PSE_map"] - output_values[, "av01450PSE_map"]
wRawVis_01pos100 = output_values[, "aud01550PSE_map"] - output_values[, "av01550PSE_map"]
wRawVis_06neg100 = output_values[, "aud06450PSE_map"] - output_values[, "av06450PSE_map"]
wRawVis_06pos100 = output_values[, "aud06550PSE_map"] - output_values[, "av06550PSE_map"]
# Check whether raw weights calculated from positive vs. negative differences are similar.
matplot(rbind( wRawVis_01neg100, wRawVis_01pos100), type = 'b')

## Summary Weights
# This is only reasonable if we can also summarize the prediction, i.e. if 
# unimodal uncertainty is independent of the unimodal standard
# Regress E2-Ebi on conflict to get w1 as slope. 
# Here, E2 is auditory, so we're calculating visual weight
# For 2 data points to regress, the slope is simply rise over run, i.e. (y2-y1)/(x2-x1)
wVis_01 = (wRawVis_01pos100 - wRawVis_01neg100) / (100--100)
wVis_06 = (wRawVis_06pos100 - wRawVis_06neg100) / (100--100)
output_values[,"wVis_01_map"] = wVis_01
output_values[,'wVis_06_map'] = wVis_06

## Conflict-specific weights
output_values[,"wVis_01deltaNeg100_map"] = wRawVis_01neg100 / -100
output_values[,"wVis_01deltaPos100_map"] = wRawVis_01pos100 / 100
output_values[,'wVis_06deltaNeg100_map'] = wRawVis_06neg100 / -100
output_values[,'wVis_06deltaPos100_map'] = wRawVis_06pos100 / 100

# Predicted weights
# Our convention for bimodal ist: Name = auditory standard = standard + delta/2
# e.g. for delta = -100: av450 = auditory450 = 500+(-100/2); visual here is 550
output_values[,"pwVis_01deltaNeg100_map"] = output_values[,"aud01450JND_map"]^2 / (output_values[,"aud01450JND_map"]^2 + output_values[,"vis550JND_map"]^2)
output_values[,"pwVis_01deltaPos100_map"] = output_values[,"aud01550JND_map"]^2 / (output_values[,"aud01550JND_map"]^2 + output_values[,"vis450JND_map"]^2)
output_values[,'pwVis_06deltaNeg100_map'] = output_values[,"aud06450JND_map"]^2 / (output_values[,"aud06450JND_map"]^2 + output_values[,"vis550JND_map"]^2)
output_values[,'pwVis_06deltaPos100_map'] = output_values[,"aud06550JND_map"]^2 / (output_values[,"aud06550JND_map"]^2 + output_values[,"vis450JND_map"]^2)

map_set = which(grepl('_map',names(output_values), fixed =T))
write.csv(output_values[, map_set],file.path(result_path,'map_Values.csv'))



