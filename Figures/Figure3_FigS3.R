rm(list=ls())

ert_cis<-read.csv("a_ERTS_CIs_1000.csv",header=T)
load('bat_ierts.rdata')
load('non_volant_ierts.rdata')

transparent_colour <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}





human_impact<-"high"
timestep<-0.5   #decrease for a "higher-resolution" simulation
ages<-seq(0,60,timestep)


 if(human_impact=='high'){
         prehuman_diversity_bats<-46
         prehuman_diversity_non_volant<-203}
 
 if(human_impact=='low'){
         prehuman_diversity_bats<-44
         prehuman_diversity_non_volant<-191}
 
 current_diversity_total_bats<-44
 future_2010_total_bats<-41
 future_2015_total_bats<-40
 future_2020_total_bats<-39
 
 current_diversity_total_non_volant<-175
 future_2010_total_non_volant<-122
 future_2015_total_non_volant<-69
 future_2020_total_non_volant<-52

mycol<-transparent_colour("orange",85)
 
 
colour_lines<-c('dodgerblue2','orange','red',mycol,"gray57")
 
## FIGURE 3
#######
### PLOTS Pre-Human and 2021

par(mfrow=c(2,2), mai = c(0.55, 0.5, 0.1, 0.2),
    oma=c(0.5,5.5,2,0),
mgp=c(2.1, 1, 0))


# 1 BATS Pre-human
## Back to pre-human
plot(NULL,NULL,xlim=c(0,10),
     ylim=c(0,max(bat_ierts$d1_vals)),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(bat_ierts$d1_vals)) {
lines(ages,bat_ierts$d1_vals[i,],type='l',
      lwd=0.1,col=colour_lines[5])
}

abline(h=prehuman_diversity_bats,col=colour_lines[1],lty=2)
abline(h=current_diversity_total_bats,col=colour_lines[2],lty=2)
rect(ert_cis$HPDL[5],-10,ert_cis$HPDH[5],100,col=mycol, border = NA )

mtext('Bats', side = 2, line = 3.5,las=2,font=2,cex=0.8  )
mtext('Return time from current to pre-human diversity', 
      side = 3, line = 1,cex=0.75,font=2 )



## 2 BATS IUCN 2021
plot(NULL,NULL,xlim=c(0,10),
     ylim=c(0,max(bat_ierts$d4_vals)),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(bat_ierts$d4_vals)) {
  lines(ages,bat_ierts$d4_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=current_diversity_total_bats,col=colour_lines[2],lty=2)
abline(h=future_2020_total_bats,col=colour_lines[3],lty=2)
rect(ert_cis$HPDL[8],-10,ert_cis$HPDH[8],100,col=mycol, border = NA )

mtext('Return time if threatened species go extinct 
     IUCN 2021', side = 3, line = 0.5,cex=0.75,font=2 )



### 3 TERRESTRIAL PRE HUMAN
## Back to pre-human
plot(NULL,NULL,xlim=c(0,30),
     ylim=c(0,500),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(non_volant_ierts$d1_vals)) {
  lines(ages,non_volant_ierts$d1_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=prehuman_diversity_non_volant,col=colour_lines[1],lty=2)
abline(h=current_diversity_total_non_volant,col=colour_lines[2],lty=2)
rect(ert_cis$HPDL[1],-20,ert_cis$HPDH[1],600,col=mycol, border = NA )

mtext('Non-volant
      mammals', side = 2, line = 3.5,
      las=2,font=2,  cex=0.8 )



## 4 TERRESTRIAL IUCN 2021
plot(NULL,NULL,xlim=c(0,30),
     ylim=c(0,500),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(non_volant_ierts$d4_vals)) {
  lines(ages,non_volant_ierts$d4_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=current_diversity_total_non_volant,col=colour_lines[2],lty=2)
abline(h=future_2020_total_non_volant,col=colour_lines[3],lty=2)
rect(ert_cis$HPDL[4],-20,ert_cis$HPDH[4],600,col=mycol, border = NA )


legend(x =1,y=500, 
       legend=c(' Pre-human diversity',
                ' Current diversity',
                ' Diversity if threatened species go extinct',
                ' ERT (confidence interval posterior dist.)'),
                lty=c(2,2,2,NA),
                pch=c(NA,NA,NA,15),
                box.col=NA,
                pt.cex = c(1,1,1,2.5),
                x.intersp = 0.3,
                cex=0.7,
                bg='transparent',
                col=colour_lines)





## FIGURE S3

########
### 2010 and 2015

par(mfrow=c(2,2), 
    mai = c(0.55, 0.5, 0.1, 0.2),
    oma=c(0.5,5.5,2,0),
    mgp=c(2.1, 1, 0))

## 5 BATS 2010
plot(NULL,NULL,xlim=c(0,10),
     ylim=c(0,max(bat_ierts$d2_vals)),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(bat_ierts$d2_vals)) {
  lines(ages,bat_ierts$d2_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=current_diversity_total_bats,col=colour_lines[2],lty=2)
abline(h=future_2010_total_bats,col=colour_lines[3],lty=2)
rect(ert_cis$HPDL[6],-20,ert_cis$HPDH[6],500,col=mycol, border = NA )

mtext('Bats', side = 2, line = 3.5,las=2,font=2,cex=0.8  )
mtext('Return time if threatened species go extinct 
     IUCN 2010', 
      side = 3, line = 0.5,cex=0.75,font=2 )




## 6 BATS 2015
plot(NULL,NULL,xlim=c(0,10),
     ylim=c(0,max(bat_ierts$d3_vals)),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(bat_ierts$d3_vals)) {
  lines(ages,bat_ierts$d3_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=current_diversity_total_bats,col=colour_lines[2],lty=2)
abline(h=future_2015_total_bats,col=colour_lines[3],lty=2)
rect(ert_cis$HPDL[7],-20,ert_cis$HPDH[7],500,col=mycol, border = NA )

mtext('Return time if threatened species go extinct 
     IUCN 2015', side = 3, line = 0.5,cex=0.75,font=2 )



### 7 TERRESTRIAL 2010
## Back to pre-human
plot(NULL,NULL,xlim=c(0,30),
     ylim=c(0,500),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(non_volant_ierts$d2_vals)) {
  lines(ages,non_volant_ierts$d2_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=current_diversity_total_non_volant,col=colour_lines[2],lty=2)
abline(h=future_2010_total_non_volant,col=colour_lines[3],lty=2)
rect(ert_cis$HPDL[2],-20,ert_cis$HPDH[2],600,col=mycol, border = NA )

mtext('Non-volant
      mammals', side = 2, line = 3.5,
      las=2,font=2, cex=0.8 )



## 8 TERRESTRIAL  2015
plot(NULL,NULL,xlim=c(0,30),
     ylim=c(0,500),
     xlab='Time (Myr)',ylab='Total species',
     xaxs = "i",cex.lab=0.9,
     main= NA,
     cex.main=1)
for(i in 1:nrow(non_volant_ierts$d3_vals)) {
  lines(ages,non_volant_ierts$d3_vals[i,],type='l',
        lwd=0.1,col=colour_lines[5])
}

abline(h=current_diversity_total_non_volant,col=colour_lines[2],lty=2)
abline(h=future_2015_total_non_volant,col=colour_lines[3],lty=2)
rect(ert_cis$HPDL[3],-20,ert_cis$HPDH[3],600,col=mycol, border = NA )


legend(x =1,y=500, 
       legend=c(' Current diversity',
                ' Diversity if threatened species go extinct',
                ' ERT (confidence interval posterior dist.)'),
       lty=c(2,2,NA),
       pch=c(NA,NA,15),
       pt.cex = c(1,1,2.5),
       box.col=NA,
       x.intersp = 0.3,
       cex=0.7,bg='transparent',
       col=c(colour_lines[2],colour_lines[3],colour_lines[4]))

 
