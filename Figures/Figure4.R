rm(list=ls())
library(DAISIE)

plot_resolution<-10000


## ML PARS
## 1st five for M26 for terrestrial species
# 2nd five for M26 for bat species
pars<-c(0.334401932,0.28844913,Inf,0.00035675,1.464940908,
    0.334401932,0.464002713,Inf,0.034309807,1.464940908)

## separate the ML pars
terrestrial_pars<-pars[1:5]
bat_pars<-pars[6:10]

## Time settings for the future diversity simulation
timestep<-0.5   #decrease for a higher-resolution simulation
ages<-seq(0,100,timestep)


### RETURN TO PRE-HUMAN
#pre-human
total_diversity_terrestrial<-203
total_diversity_bats<-46



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


## Terrestrial
Mad_pars<-terrestrial_pars
M=780
total_diversity<-total_diversity_terrestrial
iert_results_terrestrial<-c()
for (i in 1:plot_resolution){
  start_div<-sample(1:total_diversity-1,1)
  start_endemic<-sample(0:start_div,1)
  start_non_endemic<-start_div-start_endemic
  end_div<-total_diversity
  div_difference<-end_div-start_div
  stt<-DAISIE_ExpEIN(t = ages, pars = Mad_pars, 
                     M = M, initEI= c(start_endemic,
                                      start_non_endemic))
  iert<-return_time_intercept(end_div,stt[[3]],timestep)
  
  iert_results_terrestrial<-rbind(iert_results_terrestrial,c(i,start_div,
                                                             start_endemic,
                                                             start_non_endemic,
                                                             end_div,
                                                             div_difference,iert))
  
}
colnames(iert_results_terrestrial)<-c('i','start_div','start_endemic',
                                      'start_non_endemic',
                                      'end_div','div_difference',
                                      'iert')

## BATS
Mad_pars<-bat_pars
M=220
total_diversity<-total_diversity_bats
iert_results_bats<-c()
for (i in 1:plot_resolution){
  start_div<-sample(1:total_diversity-1,1)
  start_endemic<-sample(0:start_div,1)
  start_non_endemic<-start_div-start_endemic
  end_div<-total_diversity
  div_difference<-end_div-start_div
  stt<-DAISIE_ExpEIN(t = ages, pars = Mad_pars, 
                     M = M, initEI= c(start_endemic,
                                      start_non_endemic))
  iert<-return_time_intercept(end_div,stt[[3]],timestep)
  
  iert_results_bats<-rbind(iert_results_bats,c(i,start_div,
                                               start_endemic,
                                               start_non_endemic,
                                               end_div,
                                               div_difference,iert))
  
}
colnames(iert_results_bats)<-c('i','start_div','start_endemic',
                               'start_non_endemic',
                               'end_div','div_difference',
                               'iert')



#########################
# RETURN TO CONTEMPORARY
# contemporary diversity
total_diversity_terrestrial<-175
total_diversity_bats<-44


## Terrestrial
Mad_pars<-terrestrial_pars
M=780
total_diversity<-total_diversity_terrestrial
iert_results_terrestrial_2<-c()
for (i in 1:plot_resolution){
  start_div<-sample(1:total_diversity-1,1)
  start_endemic<-sample(0:start_div,1)
  start_non_endemic<-start_div-start_endemic
  end_div<-total_diversity
  div_difference<-end_div-start_div
  stt<-DAISIE_ExpEIN(t = ages, pars = Mad_pars, 
                     M = M, initEI= c(start_endemic,
                                      start_non_endemic))
  iert<-return_time_intercept(end_div,stt[[3]],timestep)
  
  iert_results_terrestrial_2<-rbind(iert_results_terrestrial_2,c(i,start_div,
                                                             start_endemic,
                                                             start_non_endemic,
                                                             end_div,
                                                             div_difference,iert))
  
}
colnames(iert_results_terrestrial_2)<-c('i','start_div','start_endemic',
                                      'start_non_endemic',
                                      'end_div','div_difference',
                                      'iert')


## BATS
Mad_pars<-bat_pars
M=220
total_diversity<-total_diversity_bats
iert_results_bats_2<-c()
for (i in 1:plot_resolution){
  start_div<-sample(1:total_diversity-1,1)
  start_endemic<-sample(0:start_div,1)
  start_non_endemic<-start_div-start_endemic
  end_div<-total_diversity
  div_difference<-end_div-start_div
  stt<-DAISIE_ExpEIN(t = ages, pars = Mad_pars, 
                     M = M, initEI= c(start_endemic,
                                      start_non_endemic))
  iert<-return_time_intercept(end_div,stt[[3]],timestep)
  
  iert_results_bats_2<-rbind(iert_results_bats_2,c(i,start_div,
                                               start_endemic,
                                               start_non_endemic,
                                               end_div,
                                               div_difference,iert))
  
}
colnames(iert_results_bats_2)<-c('i','start_div','start_endemic',
                               'start_non_endemic',
                               'end_div','div_difference',
                               'iert')



#########################
### RETURN TO HALF OF CURRENT


## Terrestrial
Mad_pars<-terrestrial_pars
M=780
total_diversity<-total_diversity_terrestrial/2
iert_results_terrestrial_3<-c()
for (i in 1:plot_resolution){
  start_div<-sample(1:total_diversity-1,1)
  start_endemic<-sample(0:start_div,1)
  start_non_endemic<-start_div-start_endemic
  end_div<-total_diversity
  div_difference<-end_div-start_div
  stt<-DAISIE_ExpEIN(t = ages, pars = Mad_pars, 
                     M = M, initEI= c(start_endemic,
                                      start_non_endemic))
  iert<-return_time_intercept(end_div,stt[[3]],timestep)
  
  iert_results_terrestrial_3<-rbind(iert_results_terrestrial_3,c(i,start_div,
                                                                 start_endemic,
                                                                 start_non_endemic,
                                                                 end_div,
                                                                 div_difference,iert))
  
}
colnames(iert_results_terrestrial_3)<-c('i','start_div','start_endemic',
                                        'start_non_endemic',
                                        'end_div','div_difference',
                                        'iert')


## BATS
Mad_pars<-bat_pars
M=220
total_diversity<-total_diversity_bats/2
iert_results_bats_3<-c()
for (i in 1:plot_resolution){
  start_div<-sample(1:total_diversity-1,1)
  start_endemic<-sample(0:start_div,1)
  start_non_endemic<-start_div-start_endemic
  end_div<-total_diversity
  div_difference<-end_div-start_div
  stt<-DAISIE_ExpEIN(t = ages, pars = Mad_pars, 
                     M = M, initEI= c(start_endemic,
                                      start_non_endemic))
  iert<-return_time_intercept(end_div,stt[[3]],timestep)
  
  iert_results_bats_3<-rbind(iert_results_bats_3,c(i,start_div,
                                                   start_endemic,
                                                   start_non_endemic,
                                                   end_div,
                                                   div_difference,iert))
  
}
colnames(iert_results_bats_3)<-c('i','start_div','start_endemic',
                                 'start_non_endemic',
                                 'end_div','div_difference',
                                 'iert')


colours_ierts<-c("dodgerblue","#F0E442")


par(mfrow=c(1,2),mar=c(3.2,3,2,1),
    mgp=c(2, 1, 0))

plot(iert_results_terrestrial[,'div_difference'],iert_results_terrestrial[,'iert'],
     xlim=c(0,total_diversity_terrestrial),ylim=c(0,80),xlab='Number of species lost',
     ylab='ERT (Myr)',pch=16,cex=0.3,cex.lab=0.8,
     xaxs='i',yaxs='i',
     main='Non-volant mammals',cex.main=1)

points(iert_results_terrestrial_2[,'div_difference'],iert_results_terrestrial_2[,'iert'],pch=16,cex=0.3,col=colours_ierts[1])
points(iert_results_terrestrial_3[,'div_difference'],iert_results_terrestrial_3[,'iert'],pch=16,cex=0.3,col=colours_ierts[2])

points(28,3.12,pch=1,cex=1.1)
points(53,7.5,col=colours_ierts[1],pch=1,cex=1.1)
points(107,19.35,col=colours_ierts[1],pch=0,cex=1.1)
points(122,24.23,col=colours_ierts[1],pch=2,cex=1.1)

abline(1:250,1:250,lty=2)

plot(iert_results_bats[,'div_difference'],iert_results_bats[,'iert'],
     xlim=c(0,total_diversity_bats),ylim=c(0,15),xlab='Number of species lost',
     ylab="ERT (Myr)",pch=16,cex=0.3,cex.lab=0.8,
     xaxs='i',yaxs='i',
     main='Bats',cex.main=1)
abline(1:50,1:50,lty=2)
points(iert_results_bats_2[,'div_difference'],iert_results_bats_2[,'iert'],pch=16,cex=0.3,col=colours_ierts[1])
points(iert_results_bats_3[,'div_difference'],iert_results_bats_3[,'iert'],pch=16,cex=0.3,col=colours_ierts[2])

points(2,1.27,pch=1,cex=1.1)
points(2,1.1,col=colours_ierts[1],pch=1,cex=1.1)
points(3,1.59,col=colours_ierts[1],pch=0,cex=1.1)
points(5,2.48,col=colours_ierts[1],pch=2,cex=1.1)


legend(x =22.5,y=7, legend=c('x = y',
                            'Return to:',
                            'Pre-human diversity',
                            'Current diversity',
                            'Half of current diversity,',
                            'Values in the data:',
                            'Lost since human arrival',
                            'Threatened IUCN 2010',
                            'Threatened IUCN 2015',
                            'Threatened IUCN 2021'
                            ),
       lty=c(2,NA,NA,NA,NA,NA,NA,NA,NA,NA),
       pch=c(NA,NA,16,16,16,NA,1,1,0,2),box.col=NA,
       x.intersp = 0.3,
       cex=0.6,bg='transparent',
       col=c('black',NA,'black',colours_ierts[1],colours_ierts[2],NA,
                            'black',colours_ierts[1],colours_ierts[1],colours_ierts[1]))


