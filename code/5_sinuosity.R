# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code to calculate sinuosity of successional trajectories ####
##Recreate Figure 4 ab and Supplemental Fig 4 in MS####

library(plyr)
library(dplyr)
library(ggplot2)
library(ape)
library(devtools)
library(tidyr)
library(vegan) ##must be the most recent version from github##
library(viridis)
library(stringr)
library(here)
library(tidyverse)
library(tsvr)
library(ggpubr)
library(ecotraj)
library(investr)
library(data.table)
library(trajr)
library(gtable)
library(gridExtra)
library(codyn)
library(indicspecies)
library(BiodiversityR)

source(here("code/1_functions.r"))
source(here("code/2_upload_format_data.r"))
source(here("code/3_PCoA_Permanova.r"))
source(here("code/4_centroids.r"))


#### CALCULATE Sinuosity of trajectories over whole time series ####

#https://cran.rstudio.com/web/packages/trajr/vignettes/trajr-vignette.html#

#make each experiment, nutrient, plot, its own combo#
PCOAandplotsboth2<-PCOAandplotsboth%>% mutate(expntrtplot= paste(exp, ntrt, plot, sep = '_'))

PCOAandplotsboth3<-PCOAandplotsboth2 %>%  ungroup() %>% dplyr::select(expntrtplot, Axis.1, Axis.2, year, -exp, -disk)


#calculate sinuoisty for each plot using the traj from coords function#
trajectoriesplot<-lapply(X = split.data.frame(x = PCOAandplotsboth3, f=PCOAandplotsboth3$expntrtplot),function(t){
  PCOAscores<-t %>%   ungroup() %>% dplyr::select(Axis.1, Axis.2, year)
  trajs <- TrajFromCoords(PCOAscores)
  sinuosity<-TrajSinuosity2(trajs)
  return(data.frame(sinuosity))
}
)

#dataframe of results#
trajectoriesplotdf<-plyr::ldply(trajectoriesplot, data.frame)

trajectoriesplotdf_2<-separate(data = trajectoriesplotdf, col = .id, into = c("exp", "ntrt", "plot"))

#mean sinuoisty for each exp_ntrt##

sinuosityavg<-trajectoriesplotdf_2%>% group_by(exp, ntrt) %>% dplyr::summarise(sinuosityavg=mean(sinuosity, na.rm=TRUE),sdsinuosiy=sd(sinuosity, na.rm=TRUE),sesinuosity=sd(sinuosity, na.rm=TRUE)/sqrt(n()), uprconfint=sinuosityavg+(1.96*sesinuosity),lwrconfint=sinuosityavg-(1.96*sesinuosity), n=n())


#### CALCULATE Sinuosity across decades######

#split by decade##

PCOAandplotsfirstdecade<-subset(PCOAandplotsboth, as.numeric(as.character(year))<=1992)

PCOAandplotsseconddecade<-subset(PCOAandplotsboth, as.numeric(as.character(year))>1992)


#make each experiment, nutrient, plot, its own combo#
PCOAandplotsfirstdecade2<-PCOAandplotsfirstdecade%>% mutate(expntrtplot= paste(exp, ntrt, plot, sep = '_'))

PCOAandplotsfirstdecade3<-PCOAandplotsfirstdecade2 %>%  ungroup() %>% select(expntrtplot, Axis.1, Axis.2, year, -exp, -disk)

PCOAandplotsseconddecade2<-PCOAandplotsseconddecade%>% mutate(expntrtplot= paste(exp, ntrt, plot, sep = '_'))

PCOAandplotsseconddecade3<-PCOAandplotsseconddecade2 %>%  ungroup() %>% select(expntrtplot, Axis.1, Axis.2, year, -exp, -disk)

#calculate sinuoisty for each plot#

#
trajectoriesplotfirst<-lapply(X = split.data.frame(x = PCOAandplotsfirstdecade3, f=PCOAandplotsfirstdecade3$expntrtplot),function(t){
  PCOAscores<-t %>%   ungroup() %>% select(Axis.1, Axis.2, year)
  trajs <- TrajFromCoords(PCOAscores)
  sinuosity<-TrajSinuosity2(trajs)
  return(data.frame(sinuosity))
}
)

trajectoriesplotsecond<-lapply(X = split.data.frame(x = PCOAandplotsseconddecade3, f=PCOAandplotsseconddecade3$expntrtplot),function(t){
  PCOAscores<-t %>%   ungroup() %>% select(Axis.1, Axis.2, year)
  trajs <- TrajFromCoords(PCOAscores)
  sinuosity<-TrajSinuosity2(trajs)
  return(data.frame(sinuosity))
}
)


#dataframe of results#
trajectoriesplotdffirst<-plyr::ldply(trajectoriesplotfirst, data.frame)
trajectoriesplotdffirst_2<-separate(data = trajectoriesplotdffirst, col = .id, into = c("exp", "ntrt", "plot"))

trajectoriesplotdfsecond<-plyr::ldply(trajectoriesplotsecond, data.frame)

trajectoriesplotdfsecond_2<-separate(data = trajectoriesplotdfsecond, col = .id, into = c("exp", "ntrt", "plot"))

# add column for decades#

trajectoriesplotdffirst_2<-trajectoriesplotdffirst_2 %>% mutate(decade="1982-1992")

trajectoriesplotdfsecond_2<-trajectoriesplotdfsecond_2 %>% mutate(decade="1993-2004")

#combine them together#
trajectoriesboth<-rbind(trajectoriesplotdffirst_2, trajectoriesplotdfsecond_2)

#mean sinuoisty for each exp_ntrt for both decades#

sinuosityavgbydecade<-trajectoriesboth%>% group_by(exp, ntrt, decade) %>% dplyr::summarise(sinuosityavg=mean(sinuosity, na.rm=TRUE),sdsinuosiy=sd(sinuosity, na.rm=TRUE),sesinuosity=sd(sinuosity, na.rm=TRUE)/sqrt(n()), n=n())

###### RECREATE FIGURE 4 ab, Sinuosity over time series ##########

sinuosityavg2<-sinuosityavg %>% mutate(ntrt2=ntrt)

# plot specifications #
sinuosityavg2$ntrt2  <- mapvalues(sinuosityavg2$ntrt2, from=c("9", "1", "2", "4","6"), to=c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))
sinuosityavg2$ntrt2 <- factor(sinuosityavg2$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

sinuosityavg2$exp<-mapvalues(sinuosityavg2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

sinuosityavg2$exp<- factor(sinuosityavg2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

sinuosityavg2<-sinuosityavg2 %>% mutate(abclabel=exp)

levels(sinuosityavg2$abclabel) <- c( "(a) ", "(b) ")

#plot#
sinuosityplot<- ggplot(data=sinuosityavg2, aes(x=ntrt2, y=sinuosityavg, group=ntrt2, color=ntrt2)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lwrconfint, ymax=uprconfint), size=1)+
  facet_grid(~exp)+ 
  xlab("")+
  ylab(expression('Average sinuosity (1982 - 2004)'))+
  mytheme+
  theme(legend.position="none")+
  theme(legend.title = element_blank()) + 
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  theme(axis.text.x = element_text(angle = 25,vjust = 0.5))+
  theme(strip.text = element_text(face="bold"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  geom_text(aes(x = 0.8, y = 8.2, label = abclabel, group =     abclabel),size = 4, color= "black", check_overlap = T)+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ylim(4.8, 8.3)




##### RECREATE FIGURE S4, Sinuosity across decades ####

sinuosityavgbydecade2<-sinuosityavgbydecade %>% mutate(ntrt2=ntrt)

#plot specifications
sinuosityavgbydecade2$ntrt2  <- mapvalues(sinuosityavgbydecade2$ntrt2, from=c("9", "1", "2", "4","6"), to=c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

sinuosityavgbydecade2$ntrt2 <- factor(sinuosityavgbydecade2$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))


sinuosityavgbydecade2$exp<-mapvalues(sinuosityavgbydecade2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

sinuosityavgbydecade2$exp<- factor(sinuosityavgbydecade2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

#plot# 
sinuosityplotbydecade<- ggplot(data=sinuosityavgbydecade2, aes(x=ntrt2, y=sinuosityavg, group=ntrt2, color=ntrt2)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=sinuosityavg-sesinuosity, ymax=sinuosityavg+sesinuosity), size=1)+
  facet_grid(decade~exp)+ 
  xlab("")+
  ylab(expression('Average sinuosity'))+
  mytheme+
  theme(legend.position="none")+
  theme(legend.title = element_blank()) + 
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  theme(axis.text.x = element_text(angle = 25,vjust = 0.5))+
  theme(strip.text = element_text(face="bold"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())
