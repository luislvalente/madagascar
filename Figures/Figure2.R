rm(list=ls())
library(ggplot2)
library(tidyverse)
library(ggsci)
library(gridExtra)

df1<-read.csv("a_ERTS_CIs_1000.csv",header=T)
all_iucns<-read_csv("a_IUCN_for_barplots.csv")


### 1 ERTs
df1$Taxon <- factor(df1$Taxon,levels = c("Bats","Non-volant"))           
df1$Stage <- factor(df1$Stage,levels = c("Pre-human", "2010", "2015", "2021"))           

p1 <- ggplot(data=df1, aes(x=Stage, y=IERT, fill=Taxon)) +
  geom_bar(stat="identity", color="black", size=0,position=position_dodge(width = 0.85), width=0.8)+
 ylab("ERT (Myr)") +
  geom_errorbar(aes(ymin=HPDL, ymax=HPDH), width=0.2, lwd = 0.3,
                position=position_dodge(.9)) +
scale_fill_manual(values=c("#56B4E9","#E69F00"), labels = c("Bats","Non-volant mammals"))+
  scale_y_continuous(expand = c(0,0.1)) +
theme(axis.title=element_blank(),
      axis.title.y = element_text(),
      panel.background = element_rect(fill="white"))+
    geom_segment(aes(x = 0, y =0 , xend = 0, yend = 30),col="black")
p1




### 2 IUCN BARPLOTS
all_iucns$Status<-factor(all_iucns$Status, levels = c("NE+DD","Not_threatened","Threatened"))

bats<-all_iucns %>% filter(Taxon == "Bats")
non_volant<-all_iucns %>% filter(Taxon == "Non_volant")


colours_IUCN<-c('grey',"#F0E442","mediumorchid")

gbats<-ggplot(bats, aes(fill=Status, y=Values, x=IUCN_year)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN,
                    name = "IUCN red list status",
                    labels = c("Not evaluated + Data deficient", "Not Threatened",'Threatened'))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.text.x = element_text(size=7,angle = 90, vjust = 0.5, hjust=1),
        legend.position='left')

gnonvolant<-ggplot(non_volant, aes(fill=Status, y=Values, x=IUCN_year)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colours_IUCN)+
scale_y_continuous(expand = c(0,0)) +
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=7,angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        panel.background = element_rect(fill="white"))

lay <- c(1,2,3,3)

grid.arrange(gbats,gnonvolant,nrow=1)



