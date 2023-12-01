# Basic data simulator
# TODO: Adapt it to the parameter transformations of the model
#       then, same priors can be used and we can test inference on simulated data. 

# TODO: This badly needs refactoring!
rm(list=ls())
source('Psychophysics_HelperFunctions.R')
require(rlang)
require(Rlab)

n_subjects = 40
n_repetitions_per_stimulus = 16
contamination_probability = 0.05

visual_450_PSEs = rnorm(n_subjects, mean = 450, sd = 10)
visual_550_PSEs = rnorm(n_subjects, mean = 550, sd = 10)
visual_JND_mean = 350
visual_JND_sd = 10
visual_JND_rate = visual_JND_mean / visual_JND_sd^2
visual_JND_shape = visual_JND_mean * visual_JND_rate
visual_JNDs = rgamma(n_subjects, visual_JND_shape, visual_JND_rate)
visual_450_JNDs = visual_JNDs + rnorm(n_subjects, 0, 1)
visual_550_JNDs = visual_JNDs + rnorm(n_subjects, 0, 1)

auditory01_450_PSEs = rnorm(n_subjects, mean = 450, sd = 10)
auditory01_550_PSEs = rnorm(n_subjects, mean = 550, sd = 10)
auditory01_JND_mean = 250
auditory01_JND_sd = 10
auditory01_JND_rate = auditory01_JND_mean / auditory01_JND_sd^2
auditory01_JND_shape = auditory01_JND_mean * auditory01_JND_rate
auditory01_JNDs = rgamma(n_subjects, auditory01_JND_shape, auditory01_JND_rate)
auditory01_450_JNDs = auditory01_JNDs + rnorm(n_subjects, 0, 1)
auditory01_550_JNDs = auditory01_JNDs + rnorm(n_subjects, 0, 1)

auditory06_450_PSEs = rnorm(n_subjects, mean = 450, sd = 10)
auditory06_550_PSEs = rnorm(n_subjects, mean = 550, sd = 10)
auditory06_JND_mean = 400
auditory06_JND_sd = 10
auditory06_JND_rate = auditory06_JND_mean / auditory06_JND_sd^2
auditory06_JND_shape = auditory06_JND_mean * auditory06_JND_rate
auditory06_JNDs = rgamma(n_subjects, auditory06_JND_shape, auditory06_JND_rate)
auditory06_450_JNDs = auditory06_JNDs + rnorm(n_subjects, 0, 1)
auditory06_550_JNDs = auditory06_JNDs + rnorm(n_subjects, 0, 1)

audiovisual06_450_PSEs = predict_bimodal_PSE(visual_550_PSEs, auditory06_450_PSEs, visual_550_JNDs, auditory06_450_JNDs)
audiovisual06_450_JNDs = predict_bimodal_JND(visual_550_JNDs, auditory06_450_JNDs)

audiovisual06_550_PSEs = predict_bimodal_PSE(visual_450_PSEs, auditory06_550_PSEs, visual_450_JNDs, auditory06_550_JNDs)
audiovisual06_550_JNDs = predict_bimodal_JND(visual_450_JNDs, auditory06_550_JNDs)

audiovisual01_450_PSEs = predict_bimodal_PSE(visual_550_PSEs, auditory01_450_PSEs, visual_550_JNDs, auditory01_450_JNDs)
audiovisual01_450_JNDs = predict_bimodal_JND(visual_550_JNDs, auditory01_450_JNDs)

audiovisual01_550_PSEs = predict_bimodal_PSE(visual_450_PSEs, auditory01_550_PSEs, visual_450_JNDs, auditory01_550_JNDs)
audiovisual01_550_JNDs = predict_bimodal_JND(visual_450_JNDs, auditory01_550_JNDs)

comparison_values_450 = matrix(rep(c(90, 270, 360, 405, 495, 540, 630, 810),n_subjects), nrow = n_subjects, byrow = T)
comparison_values_550 = matrix(rep(c(110, 330, 440, 495, 605, 660, 770, 990),n_subjects), nrow = n_subjects, byrow = T)
comparison_values_bimodal = matrix(rep(c(100,300,400,450,550,600,700,900),n_subjects), nrow = n_subjects, byrow = T)
n_comparisons = dim(comparison_values_450)[2]

x_mean_450 = apply(comparison_values_450, 1, mean)
x_mean_550 = apply(comparison_values_550, 1, mean)

my_rbinom = function(theta) rbinom(1, n_repetitions_per_stimulus, theta)

thetas = list(
    visual_450 = pnorm(comparison_values_450, mean = visual_450_PSEs, sd = visual_450_JNDs),
    visual_550 = pnorm(comparison_values_550, mean = visual_550_PSEs, sd = visual_550_JNDs),

    auditory01_450 = pnorm(comparison_values_450, mean = auditory01_450_PSEs, sd = auditory01_450_JNDs),
    auditory01_550 = pnorm(comparison_values_550, mean = auditory01_550_PSEs, sd = auditory01_550_JNDs),

    auditory06_450 = pnorm(comparison_values_450, mean = auditory06_450_PSEs, sd = auditory06_450_JNDs),
    auditory06_550 = pnorm(comparison_values_550, mean = auditory06_550_PSEs, sd = auditory06_550_JNDs),

    audiovisual01_450 = pnorm(comparison_values_bimodal, mean = audiovisual01_450_PSEs, sd = audiovisual01_450_JNDs),
    audiovisual01_550 = pnorm(comparison_values_bimodal, mean = audiovisual01_550_PSEs, sd = audiovisual01_550_JNDs),

    audiovisual06_450 = pnorm(comparison_values_bimodal, mean = audiovisual06_450_PSEs, sd = audiovisual06_450_JNDs),
    audiovisual06_550 = pnorm(comparison_values_bimodal, mean = audiovisual06_550_PSEs, sd = audiovisual06_550_JNDs)
)

for (name in names(thetas)) {
    contamination_indices = as.logical(rbern(n = n_subjects*n_comparisons, contamination_probability))
    n_contaminated = sum(contamination_indices)
    thetas[[name]][contamination_indices] = runif(n_contaminated)
}



vis450 = matrix(sapply(thetas[['visual_450']], my_rbinom), nrow = n_subjects)
vis550 = matrix(sapply(thetas[['visual_550']], my_rbinom), nrow = n_subjects)
aud01450 = matrix(sapply(thetas[['auditory01_450']], my_rbinom), nrow = n_subjects)
aud01550 = matrix(sapply(thetas[['auditory01_550']], my_rbinom), nrow = n_subjects)
aud06450 = matrix(sapply(thetas[['auditory06_450']], my_rbinom), nrow = n_subjects)
aud06550 = matrix(sapply(thetas[['auditory06_550']], my_rbinom), nrow = n_subjects)
av01450 = matrix(sapply(thetas[['audiovisual01_450']], my_rbinom), nrow = n_subjects)
av01550 = matrix(sapply(thetas[['audiovisual01_550']], my_rbinom), nrow = n_subjects)
av06450 = matrix(sapply(thetas[['audiovisual06_450']], my_rbinom), nrow = n_subjects)
av06550 = matrix(sapply(thetas[['audiovisual06_550']], my_rbinom), nrow = n_subjects)

responses_df = data.frame(row.names = 1:n_subjects)
nresp_df = data.frame(row.names = 1:n_subjects)
comp_df = data.frame(row.names = 1:n_subjects)

comparison_names = c(100,300,400,450,550,600,700,900)
for (variable in c(substitute(vis450), substitute(vis550), substitute(aud01450), 
substitute(aud01550), substitute(aud06450), substitute(aud06550), 
substitute(av01450), substitute(av01550), substitute(av06450),substitute(av06550))) {
    for (comparison_index in 1:n_comparisons) {
        full_name = paste0(variable, comparison_names[comparison_index], sep = '')
        responses_df[full_name] = eval(variable)[,comparison_index]
        nresp_df[full_name] = matrix(n_repetitions_per_stimulus, ncol = 1, nrow = n_subjects)
        variable_name = as_string(variable)
        if (substr(variable_name, 1,2) == 'av') {
            comp_df[full_name] = comparison_values_bimodal[,comparison_index]
        } 
        else if (grepl('550', variable_name)) {
            comp_df[full_name] = comparison_values_550[,comparison_index]
        } else if (grepl('450', variable_name)) {
            comp_df[full_name] = comparison_values_450[,comparison_index]
        }
        stopifnot(dim(comp_df)==dim(responses_df))
    }   
}

prop_df = responses_df/nresp_df

write.csv(responses_df, 'responses.csv')
write.csv(nresp_df,'nresp.csv')
write.csv(comp_df,'comp.csv')
write.csv(prop_df,'prop.csv')