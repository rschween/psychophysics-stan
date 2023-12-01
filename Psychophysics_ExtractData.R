# Copyright 2021-23 Raphael Schween & Ruslan Spartakov - see repo copyright notice.

# output <-  list()
check_randomizaion = function(main_data){
  
  psychophysics_trials = which(main_data$response == 'l' | main_data$response == 'a') # identified by their admissible responses
  dp = main_data[psychophysics_trials,]
  
  dp$trialDefinition = paste0( dp$type, dp$noiseAmp, dp$standrad, dp$comparison)
  
  for (b in 0:15) print(table(dp$trialDefinition[(b*80+1):(b*80+80)]))
  
}

plot_performance_curve = function(main_data, subject){
  # separate each block into clear and unclear cases and plot those over time to 
  # get an estimate of whether participant performance decays with time or if they guess randomly
  # Expectations for a well-performing participant are  
    # - clear cases should be greater or equal to unclear cases
    # - there should be no time trend in either one or the difference between them.
  psychophysics_trials = which(main_data$response == 'l' | main_data$response == 'a')
  dp = main_data[main_data$trial_type == 'psychophysics',]
  stopifnot(lengths(dp)[1]==1280)

  clear_cases = which( dp['comparison'] < (dp['standard'] * 0.75) |
                          dp['comparison'] > (dp['standard'] * 1.25))

  unclear_cases = which( dp['comparison'] > dp['standard'] * 0.75 &
                          dp['comparison'] < dp['standard'] * 1.25)

  all_cases = list('clear_cases' = clear_cases, 'unclear_cases' = unclear_cases)

  comparison_is_larger = which(dp['comparison'] > dp['standard'])
  comparison_is_smaller = which(dp['comparison'] < dp['standard'])

  lg = list('cmpSmaller' = comparison_is_smaller, 'cmpLarger' = comparison_is_larger)

  comparison_is_last = which( dp['whichFirst'] == 'standard')
  comparison_is_first = which( dp['whichFirst'] == 'comparison')

  fs = list('cmpLast' = comparison_is_last, 'cmpFirst' = comparison_is_first)
  n_correct = matrix(nrow = 2, ncol =16)
  case_names = c('clear_cases','unclear_cases')
  n_total = matrix(nrow=1, ncol = 16)

  for (b in 1:16) {
    for (c in 1:2) {
      cases = case_names[c]
      n_correct[c,b] = 0
      n_correct[c,b] = n_correct[c,b] + sum(dp[intersect(all_cases[[cases]][(40*(b-1)+1):(40*(b-1)+40)],intersect(comparison_is_smaller,comparison_is_first)) ,'key_press'] == 'l')
      n_correct[c,b] = n_correct[c,b] + sum(dp[intersect(all_cases[[cases]][(40*(b-1)+1):(40*(b-1)+40)],intersect(comparison_is_smaller,comparison_is_last)) ,'key_press'] == 'a')
      n_correct[c,b] = n_correct[c,b] + sum(dp[intersect(all_cases[[cases]][(40*(b-1)+1):(40*(b-1)+40)],intersect(comparison_is_larger,comparison_is_first)) ,'key_press'] == 'a')
      n_correct[c,b] = n_correct[c,b] + sum(dp[intersect(all_cases[[cases]][(40*(b-1)+1):(40*(b-1)+40)],intersect(comparison_is_larger,comparison_is_last)) ,'key_press'] == 'l')
    }
  }

  plot(n_correct[1,], ylim = c(0,40), pch = 3, main = paste0('Performance curve for subject ', subject))
  points(n_correct[2,], ylim = c(0,40), col = 'grey')
}

extract_data = function(main_data, condition_df, response_df, response_names, subject) {
  # fill relevant fields in response_df with number of times comparison
  # was judged greater by modality (incl. noiseLevel), standard and comparison,
  # where the latter is unchecked. 
  
  # Select type (visual, auditory, audiovisual) and noise level ('01', '06', or None)
  subset_index = (main_data$type == condition_df$type & main_data$noiseAmp == condition_df$noiseLevel)
  current_data <- subset(main_data, subset = subset_index) # this also gets rid of preload and pause trials.
  
  # it appears the elements stored in the data frame are wrapped as lists of length 1 (probably
  # a residual of the json import) - replace relevant columns by unlisted entries. There may
  # be a smarter way to do that in the import functions. 
  current_data$stim1aud <- unlist(current_data$stim1aud)
  current_data$stim2aud <- unlist(current_data$stim2aud)
  current_data$stim1vis <- unlist(current_data$stim1vis)
  current_data$stim2vis <- unlist(current_data$stim2vis)
  
  standards <- c(450, 550)
  standardIs1 <- current_data$whichFirst == 'standard'
  standardIs2 <- current_data$whichFirst == 'comparison'
  
  stopifnot( is_empty( intersect(which(standardIs1), which(standardIs2))))
  
  # make all standards & comparisons first stimulus
  current_data$standardaud <- current_data$stim1aud
  current_data$comparisonaud <- current_data$stim1aud
  current_data$standardvis <- current_data$stim1vis
  current_data$comparisonvis <- current_data$stim1vis
  # replace those standards where 2nd stimulus is standard by 2nd stimulus
  current_data$standardaud[standardIs2] <- current_data$stim2aud[standardIs2]
  current_data$standardvis[standardIs2] <- current_data$stim2vis[standardIs2]
  # replace those comparisons where 1st stimulus is standard (and therefore 2nd stimulus is comparison) by 2nd stimulus
  current_data$comparisonaud[standardIs1] <- current_data$stim2aud[standardIs1]
  current_data$comparisonvis[standardIs1] <- current_data$stim2vis[standardIs1]
  
  table( current_data$standardaud )
  table( current_data$standardvis )
  table( current_data$comparisonaud)
  table( current_data$comparisonvis)
  
  
  current_data$cPG <- rep(NaN, nrow(current_data)) # comparison perceived greater
  current_data$cPG[standardIs1 & current_data$key_press == 'l'] <- TRUE
  current_data$cPG[standardIs1 & current_data$key_press == 'a'] <- FALSE
  current_data$cPG[standardIs2 & current_data$key_press == 'a'] <- TRUE
  current_data$cPG[standardIs2 & current_data$key_press == 'l'] <- FALSE
  # did we cover all cases?
  stopifnot(sum(is.na(current_data$cPG)) == 0)
  
  # al$relativeAud <- (al$comparisonaud-al$standardaud) / al$standardaud
  # al$relativeVis <- (al$comparisonvis-al$standardvis) / al$standardvis
  
  for (s in standards){
    # here is where we reference different standards depending on condition.
    # Visual unimodal has visual standard! 
    # Auditory and audiovisual have auditory standard (visual is complementary)
    tmp = subset(current_data, subset = current_data[ ,paste('standard',condition_df$names,sep='')] == s)
    # responses
    response_df[[1]][subject, grepl( paste0( condition_df$modality, s), response_names, fixed = TRUE)] = 
      table(tmp[[paste('comparison',condition_df$names,sep='')]],tmp$cPG)[,2]
    # nresp
    response_df[[2]][subject, grepl( paste0( condition_df$modality, s), response_names, fixed = TRUE)] = 
      apply(table(tmp[[paste('comparison',condition_df$names,sep='')]],tmp$cPG),FUN = sum, MARGIN = 1)
    # comp
    response_df[[3]][subject, grepl( paste0( condition_df$modality, s), response_names, fixed = TRUE)] = 
      as.numeric(unlist(dimnames(table(tmp[[paste('comparison',condition_df$names,sep='')]],tmp$cPG))[1]))
    # wflabel
    response_df[[4]] = response_df[[3]] / s

  }
    
  # output[conditions$modality] = 
  # # note: we use the auditory standard+Delta for comparison here. The visual is complementary.
  # if (noiseAmp == "01") {
  #   wfav01 <- list( 'WF' = table(al$relativeAud,  al$cPG))
  #   av01 <- list( '450' = with( subset(al, subset = al$standardaud == 450), table(comparisonaud, cPG)),
  #                 '550' = with( subset(al, subset = al$standardaud == 550), table(comparisonaud, cPG))
  #   )} else if (noiseAmp == "06") {
  #     wfav06 <- list( 'WF' = table(al$relativeAud,  al$cPG))
  #     av06 <- list( '450' = with( subset(al, subset = al$standardaud == 450), table(comparisonaud, cPG)),
  #                   '550' = with( subset(al, subset = al$standardaud == 550), table(comparisonaud, cPG))
  #     )}
  return(response_df)
}

aggregate_data = function(data_path, experiment_n_modalities, experiment_n_standards, experiment_n_comparisons, 
conditions, experiment_standards, experiment_comparison_values, n_practice_trials, result_path) {
  file_names = dir(data_path, pattern = '*.json$') 
  n_subjects = length(file_names)

  n_fields = experiment_n_modalities * experiment_n_standards *  experiment_n_comparisons
  column_names = paste0( rep(conditions$modality, each = experiment_n_standards*experiment_n_comparisons),
                    rep(experiment_standards, each = experiment_n_comparisons, times = experiment_n_modalities),
                    rep(experiment_comparison_values, times = experiment_n_modalities*experiment_n_standards))

  responses = setNames(data.frame(matrix(data = NA, nrow = n_subjects, ncol = n_fields)), column_names)
  n_responses = setNames(data.frame(matrix(data = NA, nrow = n_subjects, ncol = n_fields)), column_names)
  comp = setNames(data.frame(matrix(data = NA, nrow = n_subjects, ncol = n_fields)), column_names)
  outcome_dfs_list = list(responses,n_responses,comp)

  ## extract responses for all participants
  for (n in 1:n_subjects) {
    # load data  
    filename <- file_names[n]
    d = fromJSON( file.path( data_path,filename))
    d = as.data.frame(d)
    for (i in 1:nrow(d)) {
      if (is_empty(d$noiseAmp[i])) d$noiseAmp[i] = "None"
    }
    dmain = d[-c(1:n_practice_trials),] # remove practice trials
    
    # plot performance curve
    output_path = file.path(result_path,paste0('plot_performance_curve_subj',n,'.pdf', sep=''))
    plot_call = call('plot_performance_curve',dmain, n)
    plot_to_output(plot_call, plot_to_pdf, output_path)
    
    # extract data
    # within the function, audiovisual trials are identified by their auditory standard.
    # the visual standard is assumed complementary (e.g. when auditory is 450, visual is 550)
    # CAVE: Unimodal visual trials are identified by their visual standard!
    for (j in 1:experiment_n_modalities) outcome_dfs_list = extract_data(dmain, conditions[j,], outcome_dfs_list, response_names = names(responses), subject = n) # function imported from script, before.
  }
  write.csv(outcome_dfs_list[[1]], file.path(data_path,'responses.csv'))
  write.csv(outcome_dfs_list[[2]], file.path(data_path,'nresp.csv'))
  write.csv(outcome_dfs_list[[3]], file.path(data_path,'comp.csv'))
  write.csv(outcome_dfs_list[[1]]/outcome_dfs_list[[2]], file.path(data_path,'prop.csv'))
  write.csv(outcome_dfs_list[[4]], file.path(data_path,'wflabel.csv'))

  return(outcome_dfs_list)
}

# #########################################################################
# ###--- get the data for the unimodal auditory, low noise condition ---###
# #########################################################################
# # Same for auditory condition - I did not adapt the comments, here. 
# #select all rows with auditory condition and noiseAmp == 01.
# for (noiseAmp in c("01","06")) {
#   sus = (dmain$type == 'auditory' & dmain$noiseAmp == noiseAmp)
#   al <- subset(dmain, subset = sus)
#   
#   # it appears the elements stored in the data frame are wrapped as lists of length 1 (probably
#   # a residual of the json import) - replace relevnat columns by unlisted entries. There may
#   # be a smarter way to do that in the import functions. 
#   al$stim1aud <- unlist(al$stim1aud)
#   al$stim2aud <- unlist(al$stim2aud)
#   
#   standards <- c(450, 550)
#   standardIs1 <- al$stim1aud %in% standards & !(al$stim2aud %in% standards)
#   standardIs2 <- al$stim2aud %in% standards & !(al$stim1aud %in% standards)
#   
#   same450 = which(al$stim1aud == 450 & al$stim2aud == 450)
#   same550 = which(al$stim1aud == 550 & al$stim2aud == 550)
#   
#   set.seed(745)
#   first450 = sample(same450, 8) 
#   second450 = same450[!(same450 %in% first450)] 
#   first550 = sample(same550, 8) 
#   second550 = same550[!(same550 %in% first550)] 
#   
#   same4t5 = which(al$stim1aud == 450 & al$stim2aud == 550)
#   same5t4 = which(al$stim1aud == 550 & al$stim2aud == 450)
#   
#   
#   first4t5 = sample(same4t5, 8) 
#   second4t5 = same4t5[!(same4t5 %in% first4t5)] 
#   first5t4 = sample(same5t4, 8) 
#   second5t4 = same5t4[!(same5t4 %in% first5t4)] 
#   
#   standardIs1[first450] = TRUE
#   standardIs2[second450] = TRUE
#   standardIs1[first550] = TRUE
#   standardIs2[second550] = TRUE
#   standardIs1[first4t5] = TRUE
#   standardIs2[second4t5] = TRUE
#   standardIs1[first5t4] = TRUE
#   standardIs2[second5t4] = TRUE
#   
#   intersect(which(standardIs1), which(standardIs2))
#   # make all standards & comparisons vis1
#   # TODO: this may be buggy in case of missing values. 
#   al$standard <- al$stim1aud
#   al$comparison <- al$stim1aud
#   # replace those standards where vis2 is standard by vis2
#   al$standard[standardIs2] <- al$stim2aud[standardIs2]
#   # replace those comparisons where vis1 is standard (and therefore comparison is vis2) by vis2
#   al$comparison[standardIs1] <- al$stim2aud[standardIs1]
#   
#   table( al$standard )
#   table( al$comparison)
#   
#   
#   al$cPG <- rep(NaN, nrow(v)) # comparison perceived greater
#   al$cPG[standardIs1 & al$key_press == 'j'] <- TRUE
#   al$cPG[standardIs1 & al$key_press == 'f'] <- FALSE
#   al$cPG[standardIs2 & al$key_press == 'f'] <- TRUE
#   al$cPG[standardIs2 & al$key_press == 'j'] <- FALSE
#   # did we cover all cases?
#   stopifnot(sum(is.na(al$cPG)) == 0)
#   
#   al$relative <- (al$comparison-al$standard) / al$standard
#   if (noiseAmp == "01") {
#     wfa01 <- list( 'WF' = table(al$relative,  al$cPG))
#     aud01 <- list( '450' = with( subset(al, subset = al$standard == 450), table(comparison, cPG)),
#                   '550' = with( subset(al, subset = al$standard == 550), table(comparison, cPG))
#     )} else if (noiseAmp == "06") {
#       wfa06 <- list( 'WF' = table(al$relative,  al$cPG))
#       aud06 <- list( '450' = with( subset(al, subset = al$standard == 450), table(comparison, cPG)),
#                     '550' = with( subset(al, subset = al$standard == 550), table(comparison, cPG))
#       )}
#   
# } # end loop for auditory
# 
# output["aud01" = aud01]
# output["aud06" = aud06]
# output["wfa01" = wfa01]
# output["wfa06" = wfa06]
# 
# for (noiseAmp in c("01","06")) {
#   sus = (dmain$type == 'audiovisual' & dmain$noiseAmp == noiseAmp)
#   al <- subset(dmain, subset = sus)
#   
#   # it appears the elements stored in the data frame are wrapped as lists of length 1 (probably
#   # a residual of the json import) - replace relevnat columns by unlisted entries. There may
#   # be a smarter way to do that in the import functions. 
#   al$stim1aud <- unlist(al$stim1aud)
#   al$stim2aud <- unlist(al$stim2aud)
#   al$stim1vis <- unlist(al$stim1vis)
#   al$stim2vis <- unlist(al$stim2vis)
#   
#   standards <- c(450, 550)
#   standardIs1 <- al$stim1aud %in% standards & !(al$stim2aud %in% standards)
#   standardIs2 <- al$stim2aud %in% standards & !(al$stim1aud %in% standards)
#   
#   same450 = which(al$stim1aud == 450 & al$stim2aud == 450)
#   same550 = which(al$stim1aud == 550 & al$stim2aud == 550)
#   
#   set.seed(831)
#   first450 = sample(same450, 8) 
#   second450 = same450[!(same450 %in% first450)] 
#   first550 = sample(same550, 8) 
#   second550 = same550[!(same550 %in% first550)] 
#   
#   same4t5 = which(al$stim1aud == 450 & al$stim2aud == 550)
#   same5t4 = which(al$stim1aud == 550 & al$stim2aud == 450)
#   
#   
#   first4t5 = sample(same4t5, 8) 
#   second4t5 = same4t5[!(same4t5 %in% first4t5)] 
#   first5t4 = sample(same5t4, 8) 
#   second5t4 = same5t4[!(same5t4 %in% first5t4)] 
#   
#   standardIs1[first450] = TRUE
#   standardIs2[second450] = TRUE
#   standardIs1[first550] = TRUE
#   standardIs2[second550] = TRUE
#   standardIs1[first4t5] = TRUE
#   standardIs2[second4t5] = TRUE
#   standardIs1[first5t4] = TRUE
#   standardIs2[second5t4] = TRUE
#   
#   intersect(which(standardIs1), which(standardIs2))
#   # make all standards & comparisons vis1
#   # TODO: this may be buggy in case of missing values. 
#   al$standardaud <- al$stim1aud
#   al$comparisonaud <- al$stim1aud
#   al$standardvis <- al$stim1vis
#   al$comparisonvis <- al$stim1vis
#   # replace those standards where vis2 is standard by vis2
#   al$standardaud[standardIs2] <- al$stim2aud[standardIs2]
#   al$standardvis[standardIs2] <- al$stim2vis[standardIs2]
#   # replace those comparisons where vis1 is standard (and therefore comparison is vis2) by vis2
#   al$comparisonaud[standardIs1] <- al$stim2aud[standardIs1]
#   al$comparisonvis[standardIs1] <- al$stim2vis[standardIs1]
#   
#   table( al$standardaud )
#   table( al$standardvis )
#   table( al$comparisonaud)
#   table( al$comparisonvis)
#   
#   
#   al$cPG <- rep(NaN, nrow(v)) # comparison perceived greater
#   al$cPG[standardIs1 & al$key_press == 'j'] <- TRUE
#   al$cPG[standardIs1 & al$key_press == 'f'] <- FALSE
#   al$cPG[standardIs2 & al$key_press == 'f'] <- TRUE
#   al$cPG[standardIs2 & al$key_press == 'j'] <- FALSE
#   # did we cover all cases?
#   stopifnot(sum(is.na(al$cPG)) == 0)
#   
#   al$relativeAud <- (al$comparisonaud-al$standardaud) / al$standardaud
#   al$relativeVis <- (al$comparisonvis-al$standardvis) / al$standardvis
#   
#   # note: we use the auditory standard+Delta for comparison here. The visual is complementary.
#   if (noiseAmp == "01") {
#     wfav01 <- list( 'WF' = table(al$relativeAud,  al$cPG))
#     av01 <- list( '450' = with( subset(al, subset = al$standardaud == 450), table(comparisonaud, cPG)),
#                  '550' = with( subset(al, subset = al$standardaud == 550), table(comparisonaud, cPG))
#     )} else if (noiseAmp == "06") {
#       wfav06 <- list( 'WF' = table(al$relativeAud,  al$cPG))
#       av06 <- list( '450' = with( subset(al, subset = al$standardaud == 450), table(comparisonaud, cPG)),
#                    '550' = with( subset(al, subset = al$standardaud == 550), table(comparisonaud, cPG))
#       )}
#   
# } # end loop for audiovisual
# 
# output['av01' = av01]
# output['av06' = av06]
# output['wfav01' = wfav01]
# output['wfav06' = wfav06]