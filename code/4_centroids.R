# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code to calculate distance within and between centroids####
##Recreate Figure 3 abcd and Tables S1, S2, S3, S4 in MS####

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



##### Calculate the distance TO centroid within ntrt treatments####

#need to have the most recent verion of vegan from github for this to work#
#https://github.com/vegandevs/vegan/issues/331#

#for this analysis, getting rid of the years where not every field was sampeled#

exp12subset_2<-exp12subset_1%>%filter(year!=1995)%>%filter(year!=1998)%>%filter(year!=2001)%>%filter(year!=2003)

#just the veg data#

veg<-exp12subset_2%>% dplyr::select(`mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)`:`mass.above.VIOLA SP.`)

#log transformation#

vegmatrix<-as.matrix(veg)

logveg<-log(1+vegmatrix)

#using the apply function to calculation the dispersion (beta disp) for each nutrient treatment for each year#

bd_allboth<-betadisper(vegdist(logveg),exp12subset_2$expntrtyear, type = "centroid")    ####Note that you have to specify the centroid ##  ##takes a bit##

#bd_allbothEUC<-betadisper(vegdist(logveg, method="euclidean"),exp12subset_2$expntrtyear, type = "centroid")    ##same thing but in Euclidean distance###

##some squared distances are negative and changed to zero##


#making group distances into a dataframe with rownames#

bd_allboth_2<-as.data.frame(bd_allboth$group.distances, row.names = rownames(bd_allboth$group.distances))

setDT(bd_allboth_2, keep.rownames = TRUE)[]


bd_allboth_df<-separate(data = bd_allboth_2, col = rn, into = c("exp", "ntrt", "year"))

bd_allboth_df$year<-as.numeric(bd_allboth_df$year)
bd_allboth_df$ntrt<-as.factor(bd_allboth_df$ntrt)

bd_allboth_df$distances<-as.numeric(bd_allboth_df$`bd_allboth$group.distances`)

bd_allboth_df2<-bd_allboth_df%>%mutate(ntrt2 = ntrt)

#this dataframe (bd_allboth_df2) has the average distance to the centroid for each disurbance_nutrient group for each year. Distances = the average Bray dist for 18 plots to their group centroids (what is plotted in Fig 3 a) ##





#### AIC table dist TO centroid, RECREATE TABLE S1 in MS####

#test whether the distance to centroid is best fit by a linear or saturating function#

bd_allboth_df2<-bd_allboth_df2%>% mutate (year2=year-1981)

#linear fit#

distlin<-lm(distances~year2, data=bd_allboth_df2)

#saturating function#
distnonlin<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = bd_allboth_df2)

#Null###
distNull<-lm(distances~1, data=bd_allboth_df2)


rawaic<-AIC(distlin, distnonlin, distNull)
nR<-dim(bd_allboth_df2)[1]  #Sample size 
aictable(rawaic,nR)

#nonlinear is the best fit#





#### Model parameters dist TO centroid, RECREATE TABLE S2 in MS ####

#https://www.statforbiology.com/nonlinearregression/usefulequations#asymptotic_regression_model###

groupdistancesboth2_df<-bd_allboth_df2%>% mutate(year2=year-1981)

groupdistcontrol<-bd_allboth_df2%>% filter(ntrt==9)


#E002 seperate by nutrient and experiment#

controlE002<-subset(groupdistancesboth2_df, ntrt==9 & exp==2)
E002distntrt1<-subset(groupdistancesboth2_df, ntrt==1 & exp==2)
E002distntrt2<-subset(groupdistancesboth2_df, ntrt==2 & exp==2)
E002distntrt4<-subset(groupdistancesboth2_df, ntrt==4 & exp==2)
E002distntrt6<-subset(groupdistancesboth2_df, ntrt==6 & exp==2)

#E002 #Non linear model fits##

controlE002fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = controlE002)
E002distntrt1fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E002distntrt1)
E002distntrt2fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E002distntrt2)
E002distntrt4fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E002distntrt4)
E002distntrt6fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E002distntrt6)


#Seperate experiment and N treatments, E001 #

controlE001<-subset(groupdistancesboth2_df, ntrt==9 & exp==1)
E001distntrt1<-subset(groupdistancesboth2_df, ntrt==1 & exp==1)
E001distntrt2<-subset(groupdistancesboth2_df, ntrt==2 & exp==1)
E001distntrt4<-subset(groupdistancesboth2_df, ntrt==4 & exp==1)
E001distntrt6<-subset(groupdistancesboth2_df, ntrt==6 & exp==1)


#Nonlinear model fits E001#

controlE001fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = controlE001)
E001distntrt1fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E001distntrt1)
E001distntrt2fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E001distntrt2)
E001distntrt4fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E001distntrt4)
#E001distntrt6fit<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = E001distntrt6) ##doesn't work...see below###

#this (ntrt6) is weird (have to specify lrc and then it works...#

x6=E001distntrt6$year2
y6=E001distntrt6$distances

##first, specify the fit##
asympfunc <- function(x6, Asym, R0, lrc=coef(E002distntrt6fit)[3])
  Asym + (R0 - Asym) * exp(-exp(lrc) * x6)
nls(y6 ~ asympfunc(x6, Asym, R0),
    data = data.frame(x6, y6),
    start = list(Asym = 0.3,R0=0.5))


E001distntrt6fit<-nls(y6 ~ asympfunc(x6, Asym, R0),
                      data = data.frame(x6, y6),
                      start = as.list(coef(E002distntrt6fit)[c(1, 2)]))

#THESE ARE ALL THE VALUES IN TABLE S1#

summary(controlE001fit)
summary(E001distntrt1fit)
summary(E001distntrt2fit)
summary(E001distntrt4fit)
summary(E001distntrt6fit)

summary(controlE002fit)
summary(E002distntrt1fit)
summary(E002distntrt2fit)
summary(E002distntrt4fit)
summary(E002distntrt6fit)




#### RECREATE FIGURE 3 ab, dist TO centroid#####

#confidence intervals around the control only#

# confidence intervals E001#

confE001 <- as.data.frame(predFit(controlE001fit,interval = "confidence", level= 0.95))

confE001_1<-cbind(controlE001,confE001)

# confidence intervals E002#

confE002 <- as.data.frame(predFit(controlE002fit,interval = "confidence", level= 0.95))

confE002_2<-cbind(controlE002,confE002)

confboth<-rbind(confE001_1,confE002_2)

#combine and add label#
confboth<-rbind(confE001_1, confE002_2)

confboth_2<-confboth %>% mutate(abclabel=exp)
levels(confboth_2$abclabel) <- c( "(a) ", "(b) ")

#plot specifications#

groupdistcontrol<-groupdistcontrol %>% mutate(abclabel=exp)

levels(groupdistcontrol$abclabel) <- c("(a) ", "(b) ")

groupdistancesboth2_df<-groupdistancesboth2_df %>% mutate(abclabel=exp)

levels(groupdistancesboth2_df$abclabel)<- c("(a) ", "(b) ")

confboth_2$exp<-mapvalues(confboth_2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

confboth_2$exp<- factor(confboth_2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

groupdistcontrol$exp<-mapvalues(groupdistcontrol$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

groupdistcontrol$exp<- factor(groupdistcontrol$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

groupdistancesboth2_df$exp<-mapvalues(groupdistancesboth2_df$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

groupdistancesboth2_df$exp<- factor(groupdistancesboth2_df$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

groupdistancesboth2_df<-groupdistancesboth2_df%>%mutate(ntrt2 = as.factor(ntrt))

groupdistancesboth2_df$ntrt2  <- mapvalues(groupdistancesboth2_df$ntrt2, from=c("9", "1", "2", "3", "4","5","6"),
                                           to=c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))
groupdistancesboth2_df$ntrt2 <- factor(groupdistancesboth2_df$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

groupdistancesboth2_df<-groupdistancesboth2_df %>% mutate(abclabel=exp)
levels(groupdistancesboth2_df$abclabel) <- c( "(a) ", "(b) ")

ann_textdisk1 <- data.frame(year2 = c(12,12), distances =c(0.7, 0.7),exp= c("Intact in 1982 (E001)","Disturbed in 1982 (E002)"), ntrt2=c("None", "None"), abclabel = factor(c("(a) ","(b) " ),levels = c("(a) ", "(b) ")))

#plot###

disttocentroidplot<- ggplot(data=groupdistancesboth2_df) +
  geom_point(data=groupdistancesboth2_df, aes(x=as.numeric(as.character(year2)), y=distances, group=ntrt2, color=ntrt2)) +
  geom_point(data=groupdistcontrol, aes(x=as.numeric(as.character(year2)), y=distances),     color="gray45") +
  geom_smooth(data=groupdistancesboth2_df,aes(x=as.numeric(as.character(year2)), y=distances, color=ntrt2), method="nls",formula=y ~ a * x^b,method.args=list(start=c(a=0.5,b=-0.1)),se=FALSE, size = 0.5)+ 
  geom_smooth(data=groupdistcontrol,aes(x=as.numeric(as.character(year2)), y=distances), method="nls",formula=y ~ a * x^b,method.args=list(start=c(a=0.5,b=-0.1)),se=FALSE,colour = "gray45", size = 0.5)+
  geom_ribbon(data=confboth_2, aes(ymin=lwr, ymax=upr, x=year2), color=NA, alpha=0.2)+
  facet_grid(~exp)+ 
  ylab(expression('Average distance to group centroid'))+
  xlab("Time since experiment")+
  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  geom_text(aes(x = 1.5, y = 0.65, label = abclabel, group =     abclabel),size = 4,check_overlap = T)




#### Calculate the distance BETWEEN centroids##########

#pulling out the group centroids by year#

centroiddf<-as.data.frame(bd_allboth[["centroids"]])

centroiddf_2<-centroiddf[,1:2]

centroiddf_2 <- tibble::rownames_to_column(centroiddf_2, "expntrtyear")

centroiddf_2<-separate(data = centroiddf_2, col = expntrtyear, into = c("exp", "ntrt", "year"))

centroiddf_3<-centroiddf_2%>% mutate(expyear= paste(exp,year, sep = '_')) 


#grand mean centroid for all treatments each year#

grandcentroid<-centroiddf_2%>% group_by(exp, year) %>% dplyr::summarise(PcoA1centroid=mean(PCoA1, na.rm=TRUE),PcoA2centroid=mean(PCoA2, na.rm=TRUE))%>% mutate(expyear= paste(exp,year, sep = '_')) 


# merge the treatments means with the grand mean centroid for each year#

centroiddf_3<-merge(data.frame(centroiddf_3), data.frame(grandcentroid), by = "expyear", all = TRUE)

#calculate the distance from each centroid to grand centroid using pythagoreum theorem#

centroiddf_4<-centroiddf_3 %>% mutate(distgrandcent=sqrt((PcoA1centroid- PCoA1)^2+(PcoA2centroid-PCoA2)^2))

#average distance of each group to its grand mean#

centroiddists<-centroiddf_4%>% group_by(expyear) %>% dplyr::summarise(meandist=mean(distgrandcent, na.rm=TRUE))%>%separate(col = expyear, into = c("exp", "year"))


#this dataframe (centroidsdists) has the average distance between nutrient group centroids for each year, for each experiment#





#### AIC table distance BETWEEN centroids, RECREATE TABLE S3 in MS####

centroiddists2<-centroiddists%>% mutate (year2=as.numeric(as.character(year))-1981)

#linear fit#

distbtwnlin<-lm(meandist~year2, data=centroiddists2)

#saturating function#
distbtwnnonlin<- nls(meandist~ SSasymp(year2, Asym, R0, lrc), data = centroiddists2)

#NULL#
distbtwnNULL<- lm(meandist~1, data=centroiddists2)


rawaic<-AIC(distbtwnlin, distbtwnnonlin,distbtwnNULL)
nR<-dim(bd_allboth_df2)[1]  #Sample size 
aictable(rawaic,nR)

#Nonlinear is the best fit#

##### Model parameters dist BETWEEN centroids, RECREATE TABLE s4 in MS ###

#dist between centroids#
asympfitcentroidE001
asympfitcentroidE002

summary(asympfitcentroidE001)
summary(asympfitcentroidE002)






#### RECREATE FIGURE 3 cd, dist BETWEEN centroids#####

centroiddists$year<-as.numeric(as.character(centroiddists$year))

centroiddists_2<-centroiddists%>% mutate(year2=year-1981)

#subset by exp#

centroiddists_E001<-subset(centroiddists_2, exp==1)
centroiddists_E002<-subset(centroiddists_2, exp==2)


#https://www.r-bloggers.com/2021/07/asymptotic-confidence-intervals-for-nls-regression-in-r/#

asympfitcentroidE001<- nls(meandist~ SSasymp(year2, Asym, R0, lrc), data = centroiddists_E001)
asympfitcentroidE002 <- nls(meandist~ SSasymp(year2, Asym, R0, lrc), data = centroiddists_E002)

#confidence intervals around fit#

new.dataE001 <- data.frame(year2=centroiddists_E001$year2)
new.dataE002 <- data.frame(year2=centroiddists_E002$year2)

confE001_centroid <- as_tibble(predFit(asympfitcentroidE001, newdata = new.dataE001, interval = "confidence", level= 0.95)) %>% mutate(year2 = new.dataE001$year2)

confE001_centroiddf<-as.data.frame(confE001_centroid)

confE001_centroiddf2<-cbind(centroiddists_E001, confE001_centroiddf)

confE002_centroid <- as_tibble(predFit(asympfitcentroidE002, newdata = new.dataE002, interval = "confidence", level= 0.95)) %>% mutate(year2 = new.dataE002$year2)

confE002_centroiddf<-as.data.frame(confE002_centroid)

confE002_centroiddf2<-cbind(centroiddists_E002, confE002_centroiddf)

#combine togeter and add abc label#
confcentroidboth<-rbind(confE001_centroiddf2, confE002_centroiddf2)

#plot specifications#

centroiddists_2$exp<-mapvalues(centroiddists_2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

centroiddists_2$exp<- factor(centroiddists_2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


centroiddists_3<-centroiddists_2%>% mutate(abclabel=exp)

levels(centroiddists_3$abclabel) <- c( "(c) ", "(d) ")


confcentroidboth$exp<-mapvalues(confcentroidboth$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

confcentroidboth$exp<- factor(confcentroidboth$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


confcentroidboth_2<-confcentroidboth[,1:7] %>% mutate(abclabel=exp)

levels(confcentroidboth_2$abclabel) <- c( "(c) ", "(d) ")

#plot#

distbetweencentroids<-ggplot(data=centroiddists_3,aes(x=as.numeric(as.character(year2)),y=meandist)) +
  geom_point() +
  facet_grid(~exp) +
  ylab(expression('Average distance between group centroids'))+
  xlab("Time since experiment (years)")+
  geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), colour = "black", size = 0.5, se=F, fullrange=T)+
  geom_ribbon(data=confcentroidboth_2, aes(ymin=lwr, ymax=upr), color=NA, fill = "gray1", alpha=0.2)+
  mytheme+
  theme(legend.position="none")+
  theme(strip.text = element_blank())+
  theme(strip.background = element_blank())+
  geom_text(aes(x = 1.5, y = 0.18, label = abclabel, group =     abclabel),size = 4,check_overlap = T)


