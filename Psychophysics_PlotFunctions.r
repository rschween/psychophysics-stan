### Licensing and copyright: see repository LICENSE note. ###

### define functions to be used in plotting ###
initialize_plot <- function(xlim, ylim = c(0,1), title = "MAP Curve", ylab = "Proportion 'Long' Response", xlab = "Test Interval (ms)") {
  # initialize a new plot with specifications.
  par(mar=c(1,2,2,0),oma=c(5,5,1,1)) #specify the appearance of the margins
  plot(NULL, xlim = xlim, ylim = ylim, xaxt = 'n')
  title(title)
  mtext(ylab,side=2,line=2,outer=T,cex=1.4)
  mtext(xlab,side=1,outer=T,line=3,cex=1.4)
}

plot_to_output <- function (plot_call, plot_to_pdf = FALSE, output_path = '', subj = NA) {
  # Note: Input 'subj' is needed for plot_calls substituting 'subj'
  if (plot_to_pdf == TRUE) pdf( output_path, paper = 'a4r', family = 'sans', pointsize = 12, bg = 'transparent', fonts = 'sans') else dev.new(width=10,height=5)
  eval(plot_call)
  if (plot_to_pdf == TRUE) dev.off()  
}

plot_participants_indexed_as_subj <- function ( participant_indices = 1, plot_call, plot_to_pdf = FALSE, output_path_stump) {
  # works if plot_call is an expression that indexes participants as 'subj'
  for (subj in participant_indices) {
    output_path = paste0( output_path_stump, subj,'.pdf', sep='')
    plot_to_output(plot_call, plot_to_pdf, output_path, subj)
  }
}

plot_predictions <- function(participantData,i,value = 'JND',ylim, ylab, standards) {
  
  par(mar=c(1,2,2,0),oma=c(5,5,1,1))
  cols = matrix(c('cyan3','deeppink','darkgoldenrod2','black', 'cyan3','deeppink4','darkgoldenrod2','black' ), nrow = 2, byrow = T )
  types = c('vis','aud','pauv','av')
  noises = matrix(c('','01','01','01','','06','06','06'), nrow = 2, byrow = T)
  lins = c(2,1,1)
  
  initialize_plot(xlim = c(0.5,2.5), ylim = ylim, ylab = ylab,  xlab = 'Condition', title = paste0('Participant ',i))
  offsets = seq(-0.2,0.2,length.out = experiment_n_standards)
  jitters = seq(-.05,.05, length.out = length(types))

  valueSet = which(grepl(value,names(participantData), fixed = T))
  mapSet = which(grepl('_map',names(participantData), fixed = T))
  lowerSet = which(grepl('_lower',names(participantData), fixed = T))
  upperSet = which(grepl('_upper',names(participantData), fixed = T))  
  for (n in 1:2) {
    noise = noises[n,]
    for (s in 1:experiment_n_standards) {
        for (t in 1:length(types)) {
          type = types[t]
          
          # use auditory standard as reference but pair it with other standard vor visual. 
          if (type != 'vis') {
            standardSet = which(grepl(experiment_standards[s],names(participantData), fixed = T))
          } else {
            standardSet = which(grepl(other_standards[s],names(participantData), fixed = T))
          }
          
          # define indices
          typeSet = which(grepl(paste0(type,noise[t]),names(participantData), fixed = T))

          
          xs = n+offsets[s]+jitters[t]
          mainSet = intersect(valueSet,intersect(typeSet, standardSet))
          if (type != 'pauv') {
            arrows( xs, participantData[, intersect(lowerSet,mainSet)], xs, participantData[, intersect(upperSet,mainSet)],
                    length = 0.01, col = cols[n,t], angle = 90, code = 3)  
          }
          points( xs, participantData[, intersect(mapSet,mainSet)], 
                  col = cols[n,t], pch = t+14)
      }
    }
  }
  axis(1, at = c(1,2), labels = c('Noise 0.1','Noise 0.6'), tick = F)
  axis(1, at = c(1,2)+rep(offsets,each=2), labels = rep(c('aud: 450','aud: 550'),each = 2), padj = -1)
  legend('topleft', c('visual','auditory','predicted','audiovisual'),col = cols[1,], lty = 1)
}

plot_predictions_mean <- function(participantData, value = 'JND', ylim, ylab, standards) {
  
  ylab = paste(value,' (ms, +/-SD)')
  par(mar=c(1,2,2,0),oma=c(5,5,1,1))
  cols = matrix(c('cyan3','deeppink','darkgoldenrod2','black', 'cyan3','deeppink4','darkgoldenrod2','black' ), nrow = 2, byrow = T )
  types = c('vis','aud','pauv','av')
  noises = matrix(c('','01','01','01','','06','06','06'), nrow = 2, byrow = T)
  lins = c(2,1,1)
  
  initialize_plot(xlim = c(0.5,2.5), ylim = ylim, ylab = ylab,  xlab = 'Condition', title = 'Mean across participants')
  offsets = seq(-0.2,0.2,length.out = experiment_n_standards)
  jitters = seq(-.05,.05, length.out = length(types))
  
  valueSet = which(grepl(value,names(participantData), fixed = T))
  mapSet = which(grepl('_map',names(participantData), fixed = T))
  lowerSet = which(grepl('_lower',names(participantData), fixed = T))
  upperSet = which(grepl('_upper',names(participantData), fixed = T))
  
  for (n in 1:2) {
    noise = noises[n,]
    for (s in 1:experiment_n_standards) {
      for (t in 1:length(types)) {
        type = types[t]
        
        # use auditory standard as reference but pair it with other standard vor visual. 
        if (type != 'vis') {
          standardSet = which(grepl(experiment_standards[s],names(participantData), fixed = T))
        } else {
          standardSet = which(grepl(other_standards[s],names(participantData), fixed = T))
        }
        
        # define indices
        typeSet = which(grepl(paste0(type,noise[t]),names(participantData), fixed = T))
        
        
        xs = n+offsets[s]+jitters[t]
        mainSet = intersect(valueSet,intersect(typeSet, standardSet))
      
        arrows( xs, mean(participantData[, intersect(mapSet,mainSet)])-sd(participantData[, intersect(mapSet,mainSet)]), 
                xs, mean(participantData[, intersect(mapSet,mainSet)])+sd(participantData[, intersect(mapSet,mainSet)]),
                length = 0.01, col = cols[n,t], angle = 90, code = 3)  
        points( xs, mean(participantData[, intersect(mapSet,mainSet)]), 
                col = cols[n,t], pch = t+14)
      }
    }
  }
  axis(1, at = c(1,2), labels = c('Noise 0.1','Noise 0.6'), tick = F)
  axis(1, at = c(1,2)+rep(offsets,each=2), labels = rep(c('aud: 450','aud: 550'),each = 2), padj = -1)
  legend('topleft', c('visual','auditory','predicted','audiovisual'),col = cols[1,], lty = 1)
}

plot_gaussian = function(participant_data,i = 1, standards) {
  scale <- seq(0,1100, by=.1)
  cols = matrix(c('cyan3','deeppink','darkgoldenrod2','black', 'cyan3','deeppink4','darkgoldenrod2','black' ), nrow = 2, byrow = T )
  types = c('vis','aud','pauv','av')
  noises = matrix(c('','01','01','01','','06','06','06'), nrow = 2, byrow = T)
  myStandards = matrix(c('450','550','550','550','550','450','450','450'), nrow = 2, byrow = T)
  lins = c(2,1,1,3)
  
  layout(matrix(c(1,2,3,4),nrow = 2))
  for (n in 1:2){
    noise = noises[n,]
    for (s in 1:experiment_n_standards){ 
      plot(NA, ylim = c(0,0.003),  xlim = range(scale), ylab = 'Likelihood', xlab = 'Stim. location', yaxt="n")
      
      for (t in 1:length(types)) {
        type = types[t]
        mod = paste(type,noise[t],myStandards[s,t],sep='')
        probs = dnorm(scale, output_values[[paste(mod,'PSE_map', sep ='')]][i], sqrt(2)*output_values[[paste(mod,'JND_map', sep ='')]][i] )
        points(scale, probs,
               type ='l', lty = lins[t], col = cols[s,t], lwd = 2)
        segments(y0=0,y1=max(probs),x0=output_values[[paste(mod,'PSE_map', sep ='')]][i],
                 lty=lins[t], col = cols[s,t])
      }
      title(main=paste0('Aud. std.: ', experiment_standards[s],' ms; Noise amp.: ', 0.1*as.numeric(noise[4])))
      title(paste0('Subject ',i,' MAP estimates.'), outer=T   )
      if ((n*s) == 4) legend('topleft',c('visual','auditive','predicted','audiovisual'), lty = lins, col = cols[1,], ncol = 2)
    }
  }
  
  
}

plot_group_contamin <- function(est, i, standards, modalities) {
	# layout(matrix(1:2,1,2,byrow=T))
  layout(matrix(c(1,1,2,3,3,4), ncol = 2))
	par(mar=c(1,2,2,0),oma=c(5,5,1,1))
  cols = c('cyan3','deeppink','deeppink4')
  lins = c(2,1,1)
  scale <- seq(0,1000, by=.1)
	for (s in standards) {
		plot(NA, main = paste0('Standard: ', s,' ms'), xlab="",ylab="", ylim=c(0,1), xlim = c(0,1000), yaxt="n",xaxt="n")
		
		for (m in 1:3) {
		  mod = modalities[m]
		  # scale <- seq(est[[paste(mod,s,sep='')]]$x[i,1],x[i,est[[paste(mod,s,sep='')]]$nstim[i]], by=.1)
		  
			points(est[[paste(mod,s,sep='')]]$x[i,],est[[paste(mod,s,sep='')]]$rprop[i,],pch=14+m, col=cols[m], bg=c(grey(1-est[[paste(mod,s,sep='')]]$zmean[[i]])))
			lines(scale,F1(scale,i,est[[paste(mod,s,sep='')]]$alphaMAP, est[[paste(mod,s,sep='')]]$betaMAP, est[[paste(mod,s,sep='')]]$xmean),
			      type="l", lty = lins[m], col = cols[m], lwd = 2)
			segments(x0=x[i,1],x1=est[[paste(mod,s,sep='')]]$PSEmap[i]+est[[paste(mod,s,sep='')]]$JNDmap[i],y0=0.84,
			         lty=3, col = cols[m])
			segments(x0=x[i,1],x1=est[[paste(mod,s,sep='')]]$PSEmap[i],y0=0.5,
			         lty=3, col = cols[m])
			segments(y0=0,y1=0.84,x0=est[[paste(mod,s,sep='')]]$PSEmap[i]+est[[paste(mod,s,sep='')]]$JNDmap[i],
			         lty=3, col = cols[m])
			segments(y0=0,y1=0.5,x0=est[[paste(mod,s,sep='')]]$PSEmap[i],
			         lty=3, col = cols[m])
		}
  	if (s=='450') {
  	  axis(2,las=1,yaxp=c(0,1,2))
  	  axis(2,at=0.84,las=1)
  	  legend('topleft',c('visual','audio 0.1', 'audio 0.6'),pch = c(15:17), lty = lins, col = cols)
  	}
  	if (s>0) axis(1)
	  
	  # Plot the Gaussians below the psychometric curves. The layout is defined by "layout()" funciton, above. 
	  # For 2AFC, the SD of the cue's likelihood function is sqrt(2)*JND of the psychometric function (cf. Rohde et al. 2015, p. 13)
    plot(NA, xlab = "Comparison (ms)", ylim = c(0,0.003),  xlim = range(scale), yaxt="n")
    scale <- seq(0,1000, by=.1)
	 for (m in 1:3) {
	   mod = modalities[m]
	   points(scale, dnorm(scale, output_values[[paste(mod,s,'PSE_map', sep ='')]][i], sqrt(2)*output_values[[paste(mod,s,'JND_map', sep ='')]][i] ),
	        type ='l', lty = lins[m], col = cols[m], lwd = 2)
	 }
	}
	
	mtext("Proportion 'Long' Response",side=2,line=2,outer=T,cex=1.4)
	mtext("Test Interval (ms)",side=1,outer=T,line=3,cex=1.4)
	title(main = paste("Subject",as.character(i)), outer = T)
}

# plotWeights <- function(participantData) {
#   
#   ylab =  'Auditory weight'
#   dev.new()
#   par(mar=c(1,2,2,0),oma=c(5,5,1,1))
#   cols = matrix(c('cyan3','deeppink','darkgoldenrod2','black', 'cyan3','deeppink4','darkgoldenrod2','black' ), nrow = 2, byrow = T )
#   types = c('vis','aud','pauv','av')
#   noises = matrix(c('','01','01','01','','06','06','06'), nrow = 2, byrow = T)
#   lins = c(2,1,1)
#   
#   initPlot(xlim = c(0.5,2.5), ylim = ylim, ylab = ylab,  xlab = 'Condition', title = 'Mean across participants')
#   offsets = seq(-0.2,0.2,length.out = nStandards)
#   jitters = seq(-.05,.05, length.out = length(types))
#   
#   mapSet = which(grepl('PSE_map',names(participantData), fixed = T))
#   lowerSet = which(grepl('PSE_lower',names(participantData), fixed = T))
#   upperSet = which(grepl('PSE_upper',names(participantData), fixed = T))
#   
#   wheres = c('_lower','_map','_upper')
#   for (n in 1:2) {
#     noise = noises[n,]
#     for (s in 1:nStandards) {
#       for (w in 1:length(wheres)) {
#         where = wheres[w]
#         if (where == '_lower') whereSet = lowerSet
#         if (where == '_map') whereSet = mapSet
#         if (where == '_upper') whereSet = upperSet
#         
#         auditoryStandard = standards[s]
#         if (s == 1) {
#           audStandSet = which(grepl(standards[1],names(participantData), fixed = T))
#           visStandSet = which(grepl(standards[2],names(participantData), fixed = T))
#         } else if (s == 2) {
#           audStandSet = which(grepl(standards[2],names(participantData), fixed = T))
#           visStandSet = which(grepl(standards[1],names(participantData), fixed = T))
#         }
#         noiseSet = which(grepl(noise[2],names(participantData), fixed = T))
#         audSet = which(grepl('aud',names(participantData), fixed = T))
#         visSet = which(grepl('vis',names(participantData), fixed = T))
#         avSet = which(grepl('av',names(participantData), fixed = T))
#         uniAudSet = intersect( audSet, intersect(audStandSet, noiseSet )) 
#         uniVisSet = intersect( visSet, visStandSet )
#         biSet = intersect( avSet, intersect(audStandSet, noiseSet))
#           
# 
#         participantData[ , paste('aw',noise[2],auditoryStandard,where, sep = '')] = ( participantData[ , intersect(mapSet, biSet)] - participantData[ , intersect(mapSet, uniVisSet)] ) /
#           ( participantData[ , intersect(mapSet, uniAudSet)] - participantData[ , intersect(mapSet, uniVisSet)] )
#         # define indices
#         # conditionSet = which(grepl(paste0(noise[t]),names(participantData), fixed = T))
#         # 
#         # 
#         # xs = n+offsets[s]+jitters[t]
#         # mainSet = intersect(valueSet,intersect(typeSet, standardSet))
#         # 
#         # arrows( xs, mean(participantData[, intersect(mapSet,mainSet)])-sd(participantData[, intersect(mapSet,mainSet)]),
#         #         xs, mean(participantData[, intersect(mapSet,mainSet)])+sd(participantData[, intersect(mapSet,mainSet)]),
#         #         length = 0.01, col = cols[n,t], angle = 90, code = 3)
#         # points( xs, mean(participantData[, intersect(mapSet,mainSet)]),
#         #         col = cols[n,t], pch = t+14)
#       }
# 
#     }
#   }
#   axis(1, at = c(1,2), labels = c('Noise 0.1','Noise 0.6'), tick = F)
#   axis(1, at = c(1,2)+rep(offsets,each=2), labels = rep(c('aud: 450','aud: 550'),each = 2), padj = -1)
#   legend('topleft', c('visual','auditory','predicted','audiovisual'),col = cols[1,], lty = 1)
# }

# plotAggregatedResult <- function(Ergebnis, outcomeType = 'JND', xlim = c(.5, 2.5), ylim = c(0,1000)) {
#   plot(NULL, xlim = xlim, ylim = ylim)
#   axis(1)
#   axis(2)
#   k = 1
  
#   points(c(1:2),Ergebnis[k,paste(c("vis450", "vis550"),outcomeType,sep = '')], type = 'b', lty = 2, col = k) # dotted is visual
#   # nOTE: everything on the x-axis is sorted by 'corresponding' visual sitmulus duration, i.e. auditory unimodal is 550, then 450!
#   points(c(1:2),Ergebnis[k,paste(c("aud01550", "aud01450"),outcomeType,sep = '')], type = 'b', lty = 3, col = k+2) # dotted is auditory
  
#   points(c(1:2),Ergebnis[k,paste(c("predv450a550n01", "predv550a450n01"),outcomeType,sep = '')], type = 'b', lty = 4,col = k+2) # dashdot predicted
#   # note: in first naming, i used auditory as standard - here visual. Therefore inversion. 
#   points(c(1:2),Ergebnis[k,paste(c("av01550", "av01450"),outcomeType,sep = '')], type = 'b', lty = 1, col = k+2) # solid is estimated
  
#   points(c(1:2),Ergebnis[k,paste(c("aud06550", "aud06450"),outcomeType,sep = '')], type = 'b', lty = 3, col = k+1) # dotted is auditory
#   points(c(1:2),Ergebnis[k,paste(c("predv450a550n06", "predv550a450n06"),outcomeType,sep = '')], type = 'b', lty = 4,col = k+1) # dashdot predicted
#   # note: in first naming, i used auditory as standard - here visual. Therefore inversion. 
#   points(c(1:2),Ergebnis[k,paste(c("av06550", "av06450"),outcomeType,sep = '')], type = 'b', lty = 1, col = k+1) # solid is estimated
#   title(paste(outcomeType ,'by standard', sep = ' '))
#   legend('topleft',paste0('Subj',1:n_subjects), fill = c(1:n_subjects))
# }

# addPlot <- function (PSEmap,JNDmap,x, nstim, rprop, c, WF, alphaMAP, betaMAP, xmean){
#   # add a psychometric curve to the plot. 
#   scale <- seq(x[1],x[nstim], by=.1)
#   points(x[],rprop[],
#          pch=c,col=c)
#   for (g in 1:20)
#   {
#     lines(scale, F3(scale,g, alpha_sel, beta_sel, xmean),type="l",col= alpha(c, 0.05))
#   }
#   lines(scale,F1(scale,alphaMAP,betaMAP,xmean),col = c,type="l")
#   if (WF == TRUE) {
#     segments(x0=x[1],x1=PSEmap+JNDmap,y0=0.84,lty=2, col = c)
#     segments(x0=x[1],x1=PSEmap,y0=0.5,lty=2, col = c)
#     segments(y0=0,y1=0.84,x0=PSEmap+JNDmap,lty=2, col = c)
#     segments(y0=0,y1=0.5,x0=PSEmap,lty=2, col = c)
    
#   }
#   axis(1)
# }
