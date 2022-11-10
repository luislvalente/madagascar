rm(list=ls())
library(DAISIE)
## the file below can be found in the Figures folder
pars <- read.csv("best_pars_M26_posterior_1000.csv", header = T, sep = ',')
ERTs <- NULL


### SETTINGS
## high impact or low impact
human_impact<-'high'
M_all<-1000
 
## Time settings for the future diversity simulation
timestep<-0.5   #decrease for a "higher-resolution" simulation
ages<-seq(0,60,timestep)

## return_time_intercept() Author: Olle Odijk, 2020
## Function to determine intercept of return time in plots
return_time_intercept<-function(DivPreHuman,d,timestep){
  i<-0
  for(i in 1:(length(d)/timestep)){
    if (d[i]>DivPreHuman){
      #using linear interpolation between timesteps: y=ax+b -> x=(y-b)/a
      y2=d[i]
      y1=d[i-1]
      dy<-y2-y1
      x2<-(i-1)*timestep
      x1<-(i-2)*timestep
      dx <- x2-x1
      a=dy/dx
      b=y2-a*x2
      x=(DivPreHuman-b)/a
      break}
  }
  return(x)
}



for (i in 1:nrow(pars)) {
   print(i)

  
  ## NON_VOLANT
  
  Mad_pars_nv<-c(pars[i,1],pars[i,2],pars[i,3],pars[i,4],pars[i,5])
  nmax=400
  ## Based on the proportion of total bats vs terrestrial
  M_terrestrial=M_all-(M_all*0.22)
  if(human_impact=='high'){
    prehuman_diversity<-203}
  
  if(human_impact=='low'){
    prehuman_diversity<-191}
  
  current_diversity_total<-175
  current_diversity_endemic<-175
  current_diversity_non_endemic<-0
  future_2010_total<-122
  future_2010_endemic<-122
  future_2010_non_endemic<-0
  future_2015_total<-69
  future_2015_endemic<-69
  future_2015_non_endemic<-0
  future_2021_total<-52
  future_2021_endemic<-52
  future_2021_non_endemic<-0 
  
  
  ## TERRESTRIAL
  start_diversity<-c(current_diversity_endemic,current_diversity_non_endemic)
  d1<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_nv, M = M_terrestrial, initEI= start_diversity)
  a_terrestrial<-return_time_intercept(prehuman_diversity,d1[[3]],timestep)
  a_terrestrial_95<-return_time_intercept(prehuman_diversity*0.95,d1[[3]],timestep)
  
  start_diversity<-c(future_2010_endemic,future_2010_non_endemic)
  d2<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_nv, M = M_terrestrial, initEI= start_diversity)
  b_terrestrial<-return_time_intercept(current_diversity_total,d2[[3]],timestep)
  b_terrestrial_95<-return_time_intercept(current_diversity_total*0.95,d2[[3]],timestep)
  
  start_diversity<-c(future_2015_endemic,future_2015_non_endemic)
  d3<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_nv, M = M_terrestrial, initEI= start_diversity)
  c_terrestrial<-return_time_intercept(current_diversity_total,d3[[3]],timestep)
  c_terrestrial_95<-return_time_intercept(current_diversity_total*0.95,d3[[3]],timestep)
  
  start_diversity<-c(future_2021_endemic,future_2021_non_endemic)
  d4<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_nv, M = M_terrestrial, initEI= start_diversity)
  d_terrestrial<-return_time_intercept(current_diversity_total,d4[[3]],timestep)
  d_terrestrial_95<-return_time_intercept(current_diversity_total*0.95,d4[[3]],timestep)
  
  
  
  ### BATS
  
  Mad_pars_b<-c(pars[i,6],pars[i,7],pars[i,8],pars[i,9],pars[i,10])
  nmax=200
  ## Based on the proportion of total bats vs terrestrial
  M_bats=M_all*0.22
  if(human_impact=='high'){
    prehuman_diversity<-46}
  
  if(human_impact=='low'){
    prehuman_diversity<-44}
  
  current_diversity_total<-44
  current_diversity_endemic<-35
  current_diversity_non_endemic<-9
  future_2010_total<-41
  future_2010_endemic<-32
  future_2010_non_endemic<-9
  future_2015_total<-40
  future_2015_endemic<-31
  future_2015_non_endemic<-9
  future_2021_total<-39
  future_2021_endemic<-30
  future_2021_non_endemic<-9
  
  ### BATS
  

  
  start_diversity<-c(current_diversity_endemic,current_diversity_non_endemic)
  d1<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_b, M = M_bats, initEI= start_diversity)
  a_bats<-tryCatch(return_time_intercept(prehuman_diversity,d1[[3]],timestep))
  if(length(a_bats)==0){a_bats<-0}
  a_bats_95<-tryCatch(return_time_intercept(prehuman_diversity*0.95,d1[[3]],timestep))
  if(length(a_bats_95)==0){a_bats_95<-0}
  
  start_diversity<-c(future_2010_endemic,future_2010_non_endemic)
  d2<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_b, M = M_bats, initEI= start_diversity)
  b_bats<-tryCatch(return_time_intercept(current_diversity_total,d2[[3]],timestep))
  b_bats_95<-tryCatch(return_time_intercept(current_diversity_total*0.95,d2[[3]],timestep))
  if(length(b_bats_95)==0){b_bats_95<-0}
  
  start_diversity<-c(future_2015_endemic,future_2015_non_endemic)
  d3<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_b, M = M_bats, initEI= start_diversity)
  c_bats<-tryCatch(return_time_intercept(current_diversity_total,d3[[3]],timestep))
  c_bats_95<-tryCatch(return_time_intercept(current_diversity_total*0.95,d3[[3]],timestep))
  if(length(c_bats_95)==0){c_bats_95<-0}
  
  start_diversity<-c(future_2021_endemic,future_2021_non_endemic)
  d4<-DAISIE_ExpEIN(t = ages, pars = Mad_pars_b, M = M_bats, initEI= start_diversity)
  d_bats<-tryCatch(return_time_intercept(current_diversity_total,d4[[3]],timestep))
  d_bats_95<-tryCatch(return_time_intercept(current_diversity_total*0.95,d4[[3]],timestep))
  if(length(d_bats_95)==0){d_bats_95<-0}
  
  
  # ## Equilibrium value bats:
  # DAISIE_ExpEIN(t = Inf,pars = Mad_pars_b, M = M_bats)
  
  ## ADD RETURN TIMES 
  ERTs <- rbind(ERTs, c(round(a_terrestrial,2),round(a_terrestrial_95,2),round(a_bats,2),round(a_bats_95,2),
    round(b_terrestrial,2),round(b_terrestrial_95,2),round(b_bats,2),round(b_bats_95,2),
    round(c_terrestrial,2),round(c_terrestrial_95,2),round(c_bats,2),round(c_bats_95,2),
    round(d_terrestrial,2),round(d_terrestrial_95,2),round(d_bats,2),round(d_bats_95,2)))
  
  
  }

colnames(ERTs) <- c('pre_human_terrestrial','95%_pre_human_terrestrial','pre_human_bats','95%_pre_human_bats'
                         ,'IUCN_2010_terrestrial','95%_IUCN_2010_terrestrial','IUCN_2010_bats','95%_IUCN_2010_bats'
                         ,'IUCN_2015_terrestrial','95%_IUCN_2015_terrestrial','IUCN_2015_bats','95%_IUCN_2015_bats'
                         ,'IUCN_2021_terrestrial','95%_IUCN_2021_terrestrial','IUCN_2021_bats','95%_IUCN_2021_bats')

ERTs <- as.data.frame(ERTs)


####### Adding error intervals in ERT table

ERT_errors <- data.frame(matrix(vector(), 1, 16))
rownames(ERT_errors) <- c('Mean & Quantiles')
colnames(ERT_errors) <- colnames(ERTs)


ERT_errors[1,]<-paste0(round(sapply(ERTs,mean),2),
                       ' (',
                       round(sapply(ERTs,quantile,0.025),2),
                       '-',
                       round(sapply(ERTs,quantile,0.975),2),
                       ')',
                       sep="")

ERT_errors

t(ERT_errors)

