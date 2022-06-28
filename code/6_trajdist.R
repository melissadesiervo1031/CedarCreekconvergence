# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code to calculate trajectory distances ####

#Recreate Figure 5 and Supplemental Tables 5 and 6 in MS#

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
source(here("code/5_sinuosity.r"))



#### CALCULATE TRAJ DIST BETWEEN YEARS########

#https://cran.r-project.org/web/packages/ecotraj/vignettes/IntroductionETA.html#

#get the PCoA points into a Euclidean distance matrix to calc. distance between points ###

PCOAandplotsboth_2 <- PCOAandplotsboth %>% arrange(exp, field, ntrt,plot, year)%>%mutate(exp_field_ntrt_plot= paste(exp, field, ntrt,plot, sep = '_'), exp_field_ntrt_plot2 = as.numeric(as.factor(exp_field_ntrt_plot))) %>% dplyr::select(field, exp, plot, subplot, year, disk, ntrt, nadd, expntrtyear, expntrtfieldyear, exp_field_ntrt_plot, exp_field_ntrt_plot2, Axis.1, Axis.2)

#split it by E001 and E002 since there are missing years in E002#
PCOAandplotsE001<-subset(PCOAandplotsboth_2, exp==1) 
PCOAandplotsE002<-subset(PCOAandplotsboth_2, exp==2) 

#drop the years for E002##
levels(PCOAandplotsE002$year)
PCOAandplotsE002$year <- as.factor(as.character(PCOAandplotsE002$year))

EucE001 = dist(PCOAandplotsE001[,13:14])
EucE002 = dist(PCOAandplotsE002[,13:14])

#Calculate traj lengths#
trajlengthE001<-trajectoryLengths(EucE001, PCOAandplotsE001$exp_field_ntrt_plot, PCOAandplotsE001$year)

trajlengthE002<-trajectoryLengths(EucE002, PCOAandplotsE002$exp_field_ntrt_plot, PCOAandplotsE002$year)

# formatting E001#
yearnamesE001<-levels(PCOAandplotsE001$year)
trajlengthE001_df=as.data.frame(trajlengthE001)

colnames(trajlengthE001_df)<-yearnamesE001
names(trajlengthE001_df)[names(trajlengthE001_df) == '2004'] <- 'Total'

trajlengthE001_df$names <- rownames(trajlengthE001_df)

trajlengthE001_df2<-trajlengthE001_df %>%separate(names, c("exp", "field", "ntrt", "plot"), "_")

trajlengthE001_df3<-trajlengthE001_df2%>%dplyr::select(exp, field, ntrt, plot, '1982':'Total')

# formatting E002#
yearnamesE002<-levels(PCOAandplotsE002$year)
trajlengthE002_df=as.data.frame(trajlengthE002)
colnames(trajlengthE002_df)<-yearnamesE002
names(trajlengthE002_df)[names(trajlengthE002_df) == '2004'] <- 'Total'

trajlengthE002_df$names <- rownames(trajlengthE002_df)

trajlengthE002_df2<-trajlengthE002_df %>%separate(names, c("exp", "field", "ntrt", "plot"), "_")%>%dplyr::select(exp, field, ntrt, plot, '1982':'Total')

#for this df, we need to drop the columns where we don't have back to back years#
#drop: 1994 (col 17(missing 1995), drop 2000 (missing 2001) (col 21), drop 2002 (missing 2003)(col 22), drop 2003 (missing most of 2004)(col 23)###

trajlengthE002_df3 <- subset(trajlengthE002_df2, select = -c(17,21:23))

#convert wide to long#

traj_longE001 <- gather(trajlengthE001_df3, startingyear, trajdist, '1982':'2003', factor_key=TRUE)

traj_longE002 <- gather(trajlengthE002_df3, startingyear, trajdist, '1982':'1999', factor_key=TRUE)

#combine the two dataframes back together#

traj_longboth<-rbind(traj_longE001, traj_longE002)

#this data frame (traj_long_both) contain the trajectory distance in Euclidean space for each plot between years. traj dist = the dist between starting year and the next year#



##summarize the average distance between years by treatment group##

avgTraj_byyear_E001E002<-traj_longboth%>% group_by(exp,ntrt,startingyear) %>% dplyr::summarise(meantrajdist=mean(trajdist, na.rm=TRUE),sdtrajdist=sd(trajdist, na.rm=TRUE),setrajdist=sd(trajdist, na.rm=TRUE)/sqrt(n()), n=n())

avgTraj_byyear_E001E002_2<- avgTraj_byyear_E001E002 %>% mutate (year2=as.numeric(as.character(startingyear))+0.5-1981)




#### AIC Table Traj Dist, Recreate Table S5 and S6####

#AIC test whether a linear or nonlinear fit is better#
#linear fit#

trajdistlin<-lm(meantrajdist~year2, data=avgTraj_byyear_E001E002_2)


#nonlinear#

#trajdistnonlin2 <- nls(meantrajdist ~ SSasymp(year2, Asym, R0, lrc), data=avgTraj_byyear_E001E002_2) ##not converging###

trajdistnonlin2_2<- nls(meantrajdist ~ SSmicmen(year2, a, b), data=avgTraj_byyear_E001E002_2)


#simple exponential decay model#

trajdistnonlin3 <- nls(meantrajdist ~ 1 / (1 + year2^c),start = list(c = 1), data=avgTraj_byyear_E001E002_2)

#null#

trajdistnull<-lm(meantrajdist~1, data=avgTraj_byyear_E001E002_2)

rawaic<-AIC(trajdistnonlin3,trajdistnonlin2_2, trajdistlin, trajdistnull)
nR<-dim(avgTraj_byyear_E001E002_2)[1]  #Sample size 
aictable(rawaic,nR)

#linear model performs better than null and better than nonlinear model#

#linear model equations#

#subset by experiment##
avgTraj_byyear_E001<-subset(avgTraj_byyear_E001E002_2, exp="Intact in 1982 (E001)")
avgTraj_byyear_E002<-subset(avgTraj_byyear_E001E002_2, exp="Disturbed in 1982 (E002)")


modelsE001 <- dlply(avgTraj_byyear_E001, "ntrt", function(avgTraj_byyear_E001) 
  lm(meantrajdist ~ year2, data = avgTraj_byyear_E001))

# Print the summary of each model
l_ply(modelsE001, summary, .print = TRUE)


modelsE002 <- dlply(avgTraj_byyear_E002, "ntrt", function(avgTraj_byyear_E002) 
  lm(meantrajdist ~ year2, data = avgTraj_byyear_E002))

# Print the summary of each model
l_ply(modelsE002, summary, .print = TRUE)

#These are the values in table s5#



### RECREATE FIG 5, Traj dist#####

# plot specifications#

avgTraj_byyear_E001E002$exp<-mapvalues(avgTraj_byyear_E001E002$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

avgTraj_byyear_E001E002$exp<- factor(avgTraj_byyear_E001E002$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


avgTraj_byyear_E001E002_2<- avgTraj_byyear_E001E002 %>% mutate (year2=as.numeric(as.character(startingyear))+0.5-1981)

avgTraj_byyear_E001E002_2<-avgTraj_byyear_E001E002_2 %>% mutate(ntrt2=ntrt)

avgTraj_byyear_E001E002_2$ntrt2  <- mapvalues(avgTraj_byyear_E001E002_2$ntrt2, from=c("9", "1", "2", "3", "4","5","6"),to=c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

avgTraj_byyear_E001E002_2$ntrt2 <- factor(avgTraj_byyear_E001E002_2$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

#seperate the control from the others...#

avgTraj_byyear_E001E002_2<-avgTraj_byyear_E001E002_2 %>% mutate(abclabel=exp)
levels(avgTraj_byyear_E001E002_2$abclabel) <- c( "(a) ", "(b) ")

avgTrajcontrol<-subset(avgTraj_byyear_E001E002_2, ntrt2=="None")

#plot#
avgTraj_byyearlinear<- ggplot(data=avgTraj_byyear_E001E002_2) +
  geom_point(data=avgTraj_byyear_E001E002_2, aes(x=as.numeric(as.character(year2)), y=meantrajdist, group=ntrt2, color=ntrt2)) +
  geom_smooth(data=avgTraj_byyear_E001E002_2, aes(x=as.numeric(as.character(year2)), y=meantrajdist, group=ntrt2, color=ntrt2), method = "lm", formula = y ~ x , size = 1, se=FALSE)+
  geom_point(data=avgTrajcontrol, aes(x=as.numeric(as.character(year2)), y=meantrajdist),     color="gray45") +
  geom_smooth(data = avgTrajcontrol, aes(x=as.numeric(as.character(year2)), y=meantrajdist), method = "lm", formula = y ~ x, colour = "gray45", size = 0.25)+
  facet_grid(~exp)+ 
  ylab(expression('Annual community trajectory distance'))+
  xlab("Time since experiment (years)")+
  mytheme+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  geom_text(aes(x = 2.0, y = 0.25, label = abclabel, group =     abclabel),size = 4, color= "black", check_overlap = T)


#add text and legend#

avgTraj_byyear_2_linear<-avgTraj_byyearlinear+ theme(legend.position = c(0.85, 0.85))+theme(legend.text = element_text(size = 8))+ guides(color = guide_legend(override.aes = list(size = 0.5)))+theme(legend.margin = margin(0, 0, 0, 0))+theme(legend.spacing.y = unit(0, "mm"))+ theme(legend.key.size = unit(1.5, "mm"))
