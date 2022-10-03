rm(list=ls())

library(tidyverse)
library(gapminder)
library(mapplots)
library(ggplot2)
library(gridExtra)

divs<-read.csv('a_colonisation_scenario1.csv',header=T)
madagascar<-read_csv("a_IUCN_status.csv")


colours_IUCN<-c("gray88","gray48","forestgreen","blue","yellow","gold2","tomato","black")


divs_maxage<-divs[-which(divs[,"MaxAge"]=='No'),]
divs_radiation<-divs[-which(divs[,"Total.species"]<2),]


shape_endemic<-c()
for (i in 1:nrow(divs)){
  if(divs[i,'Clade.status']=='Endemic') {shape_endemic[i]<-16}
  if(divs[i,'Clade.status']=='Non endemic'){shape_endemic[i]<-17}}

divs<-cbind(divs, shape_endemic)


# 1 Plot Colonisation times and table

par(mar=c(2,0,0,0),mgp=c(0,0.5,0))
plot(NULL, NULL,
     xlim=c(-100,12),ylim=c(-2,265),xlab="",ylab='',yaxs ='i',yaxt='n',xaxt='n',bty="n",cex.axis=0.6)
segments(-divs$Mean.age,-2,-divs$Mean.age,divs$Order,col='grey',lty=1,lwd=0.4)
axis(1,at=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0),cex.axis=0.6,line=0)

abline(v=-66,col='red',lty=3,lwd=0.5)
abline(v=-88,col='red',lty=3,lwd=0.5)
segments(-divs$MINIMUM_MAXAGE,divs$Order,-divs$Percentile.0.975,divs$Order,
         lty= ifelse(divs$MaxAge == 'Yes',3,1),
         col = ifelse(divs$Bat == 'Yes','dodgerblue1','orange'))
points(-divs$Mean.age,divs$Order,pch=divs$shape_endemic,col=ifelse(divs$Bat == 'Yes','dodgerblue1','orange'),cex=1.2)

rect(1,divs$Order-2.9,4,divs$Order+3,col='chartreuse4',border=NA)
text(0.5,divs$Order,label=divs$Total.species,col='black', cex = 0.6, pos = 4)

rect(4.4,divs$Order-2.9,10.5,divs$Order+3,col='yellow2',border=NA)
text(5.6,divs$Order,label=divs$Threatened.species,col='black', cex = 0.6, pos = 4)

rect(10.9,divs$Order-2.9,13.9,divs$Order+3,col='tomato',border=NA)
text(10.4,divs$Order,label=divs$Extinct.species,col='black', cex = 0.6, pos = 4)

text(((-divs_radiation$Mean.age)+1),((divs_radiation$Order)+1),
     expression(paste("*")), cex = 1.5)
text(-divs$Mean.age-5,divs$Order+4,
     label=divs$Clade_name,
     col = ifelse('black'),
     cex = 0.6,font=ifelse(divs$Clade_name == 'Lemurs' |
                           divs$Clade_name == 'Afrosoricida' |
                           divs$Clade_name == 'Nesomyinae' |
                           divs$Clade_name == 'Eupleridae' |
                           divs$Clade_name == 'Hippos', 1,3))
text(-89, 80, "Magadascar isolated from other landmasses",srt = 90,cex=0.8,col='darkgrey')
text(-67, 100, "K-T event",srt = 90,cex=0.8,col='darkgrey')

legend(-100, 250, c("Bat lineage",
                   'Non-volant terrestrial lineage',
                   "Endemic lineage","Non-endemic lineage",
                   'Lineage not sampled in DNA-only tree, Upper bound',
                   '"Radiation" (lineage with more than 1 species on Madagascar)',
                   'IUCN status:',
                   'Not Evaluated / Data deficient',
                   'Least concern',
                   'Near threatened',
                   'Vulnerable',
                   'Endangered',
                   'Critically endangered',
                   'Extinct' ),
       pt.cex = 1.25,
       lty=c(0,0,0,0,3,NA,NA,0,0,0,0,0,0,0),
       pch=c(16,16,1,2,NA,NA,NA,15,15,15,15,15,15,15),
       box.col = 'white',
        cex = 0.7,
       col = c('dodgerblue','orange','black','black','black',NA,NA,
               colours_IUCN)
       )



# 2 Plot IUCN barplots


madagascar$Status<-factor(madagascar$Status, levels = c("NE","DD","LC","NT","VU","EN","CR","EX"))


lemurs<-madagascar %>% filter(Taxon == "Lemurs")
afrosoricida<-madagascar %>% filter(Taxon == "Afrosoricida")
eupleridae<-madagascar %>% filter(Taxon == "Eupleridae")
nesomyinae<-madagascar %>% filter(Taxon == "Nesomyinae")
bats<-madagascar %>% filter(Taxon == "Chiroptera")


# Stacked + percent

gafrosoricida<-ggplot(afrosoricida, aes(fill=Status, y=Values, x=IUCN_year)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN)+
  #scale_fill_brewer(palette="Spectral",direction=-1)+

  #scale_fill_simpsons()+
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill="white"))

glemur<-ggplot(lemurs, aes(fill=Status, y=Values, x=IUCN_year)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN)+
 theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill="white"))

gnesomyinae<-ggplot(nesomyinae, aes(fill=Status, y=Values, x=IUCN_year)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN)+
 theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill="white"))

gbats<-ggplot(bats, aes(fill=Status, y=Values, x=IUCN_year)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN)+
 theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill="white"))


geupleridae<-ggplot(eupleridae, aes(fill=Status, y=Values, x=IUCN_year)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN)+
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill="white"))


grid.arrange(gbats,gafrosoricida,glemur,gnesomyinae,geupleridae, nrow=1)



