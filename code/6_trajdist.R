# MH DeSiervo, LG Shoemaker
# Disturbance alters transience but nutrients determine equilibria during grassland succession with multiple global change drivers submitted with DOI https://doi.org/10.5061/dryad.dbrv15f5t.  


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
source(here("code/5_traj_directionality.r"))



#### CALCULATE TRAJ DIST BETWEEN YEARS########

##https://cran.r-project.org/web/packages/ecotraj/vignettes/IntroductionETA.html##

### Start with the BC dist matrix ###


##arrange the dataframe so that that columns and rows are in right order##
exp12subsetsorted<-exp12subset%>%arrange(year)%>%arrange(plot)%>%arrange(field)%>%arrange(exp)


exp12subsetsorted$year<-as.numeric(as.character(exp12subsetsorted$year))


### subset by experiment, due to different numbers of years in E001 and E002 ####

##just the plot data##
plotE001<-exp12subsetsorted%>% dplyr::select('field':'ntrt2') %>% filter(exp==1) %>% mutate(exp_field_ntrt_plot= paste(exp, field, ntrt,plot, sep = '_'))

plotE002<-exp12subsetsorted%>% dplyr::select('field':'ntrt2') %>% filter(exp==2) %>% mutate(exp_field_ntrt_plot= paste(exp, field, ntrt,plot, sep = '_'))


##just the veg data##

vegE001<-exp12subsetsorted%>% filter(exp==1)%>% dplyr::select(`mass.above.ACER NEGUNDO`:`mass.above.VIOLA SP.`)

vegE002<-exp12subsetsorted%>% filter(exp==2)%>% dplyr::select(`mass.above.ACER NEGUNDO`:`mass.above.VIOLA SP.`)


###log transformation####

vegmatrixE001<-as.matrix(vegE001)

logvegE001<-log(1+vegmatrixE001)


vegmatrixE002<-as.matrix(vegE002)

logvegE002<-log(1+vegmatrixE002)


## Bray curtis dissimiliarty matrix for each experiment##

distE001<-vegdist(logvegE001, method="bray")
BCdistallE001<-as.matrix(distE001)  ##into a matrix#
BCdistallE001df<-as.data.frame(BCdistallE001)  #into a dataframe#

BCdistallE001df2<-cbind(plotE001,BCdistallE001df) #merge plot data w/ dissim matrix##


distE002<-vegdist(logvegE002, method="bray")
BCdistallE002<-as.matrix(distE002)  ##into a matrix#
BCdistallE002df<-as.data.frame(BCdistallE002)  #into a dataframe#

BCdistallE002df2<-cbind(plotE002,BCdistallE002df) #merge plot data w/ dissim matrix##





####
trajlengthE001_fullN<-trajectoryLengths(distE001, BCdistallE001df2$exp_field_ntrt_plot, BCdistallE001df2$year)

trajlengthE002_fullN<-trajectoryLengths(distE002, BCdistallE002df2$exp_field_ntrt_plot, BCdistallE002df2$year)



####
levels(BCdistallE001df2$year)
BCdistallE001df2$year <- as.factor(as.character(BCdistallE001df2$year))

yearnamesE001<-levels(BCdistallE001df2$year)
trajlengthE001_df_fullN=as.data.frame(trajlengthE001_fullN)

colnames(trajlengthE001_df_fullN)<-yearnamesE001
names(trajlengthE001_df_fullN)[names(trajlengthE001_df_fullN) == '2004'] <- 'Total'

trajlengthE001_df_fullN$names <- rownames(trajlengthE001_df_fullN)

trajlengthE001_df2_fullN<-trajlengthE001_df_fullN %>%separate(names, c("exp", "field", "ntrt", "plot"), "_")%>%dplyr::select(exp, field, ntrt, plot, '1982':'Total')

####

levels(BCdistallE002df2$year)
BCdistallE002df2$year <- as.factor(as.character(BCdistallE002df2$year))

yearnamesE002<-levels(BCdistallE002df2$year)


trajlengthE002_df_fullN=as.data.frame(trajlengthE002_fullN)
colnames(trajlengthE002_df_fullN)<-yearnamesE002
names(trajlengthE002_df_fullN)[names(trajlengthE002_df_fullN) == '2004'] <- 'Total'

trajlengthE002_df_fullN$names <- rownames(trajlengthE002_df_fullN)

trajlengthE002_df2_fullN<-trajlengthE002_df_fullN %>%separate(names, c("exp", "field", "ntrt", "plot"), "_")%>%dplyr::select(exp, field, ntrt, plot, '1982':'Total')

###for this df, we need to drop the columns where we don't have back to back years###
##drop: 1994 (col 17(missing 1995), drop 2000 (missing 2001) (col 21), drop 2002 (missing 2003)(col 22), drop 2003 (missing most of 2004)(col 23)###

trajlengthE002_df3_fullN <- subset(trajlengthE002_df2_fullN, select = -c(17,21:23))

##convert wide to long##

traj_longE001_fullN<- gather(trajlengthE001_df2_fullN, startingyear, trajdist, '1982':'2003', factor_key=TRUE)

traj_longE002_fullN <- gather(trajlengthE002_df3_fullN, startingyear, trajdist, '1982':'1999', factor_key=TRUE)

##combine the two dataframes back together##

traj_longboth_fullN<-rbind(traj_longE001_fullN, traj_longE002_fullN)

##this data frame (traj_long_both) contain the trajectory distance in Euclidean space for each plot between years. traj dist = the dist between starting year and the next year###



##summarize the average distance between years by treatment group##

avgTraj_byyear_E001E002_fullN<-traj_longboth_fullN%>% group_by(exp,ntrt,startingyear) %>% dplyr::summarise(meantrajdist=mean(trajdist, na.rm=TRUE),sdtrajdist=sd(trajdist, na.rm=TRUE),setrajdist=sd(trajdist, na.rm=TRUE)/sqrt(n()), n=n())

avgTraj_byyear_E001E002_2_fullN<- avgTraj_byyear_E001E002_fullN %>% mutate (year2=as.numeric(as.character(startingyear))+0.5-1981)


#### AIC Table Traj Dist, Recreate Table S5 and S6####

###AIC test whether a linear or nonlinear fit is better##
###linear fit###


trajdistlin<-lm(meantrajdist~year2, data=avgTraj_byyear_E001E002_2_fullN)


##nonlinear asymptotic##

trajdistnonlin2 <- nls(meantrajdist ~ SSasymp(year2, Asym, R0, lrc), data=avgTraj_byyear_E001E002_2_fullN) ##converging###

#trajdistnonlin2_2<- nls(meantrajdist ~ SSmicmen(year2, a, b), data=avgTraj_byyear_E001E002_2)

##simple exponential decay model####

#trajdistnonlin3 <- nls(meantrajdist ~ 1 / (1 + year2^c),start = list(c = 1), data=avgTraj_byyear_E001E002_2)


###quadratic function###

trajdistquad<-lm(meantrajdist~year2+I(year2^2), data=avgTraj_byyear_E001E002_2_fullN)


##null##

trajdistnull<-lm(meantrajdist~1, data=avgTraj_byyear_E001E002_2_fullN)


rawaic<-AIC(trajdistlin, trajdistnonlin2, trajdistquad, trajdistnull)
nR<-dim(traj_longboth_fullN)[1]  #Sample size 
aictable(rawaic,nR)

##NONLINEAR ASYMPTOTIC MODEL performs better than null and better than nonlinear model###

####non linear model equations#####

###E001 seperate by nutrient and experiment####

controlE001<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==9 & exp==1)
E001distntrt1<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==1 & exp==1)
E001distntrt2<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==2 & exp==1)
E001distntrt3<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==3 & exp==1)
E001distntrt4<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==4 & exp==1)
E001distntrt5<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==5 & exp==1)
E001distntrt6<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==6 & exp==1)
E001distntrt7<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==7 & exp==1)
E001distntrt8<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==8 & exp==1)




###E002 seperate by nutrient and experiment####

controlE002<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==9 & exp==2)
E002distntrt1<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==1 & exp==2)
E002distntrt2<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==2 & exp==2)
E002distntrt3<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==3 & exp==2)
E002distntrt4<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==4 & exp==2)
E002distntrt5<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==5 & exp==2)
E002distntrt6<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==6 & exp==2)
E002distntrt7<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==7 & exp==2)
E002distntrt8<-subset(avgTraj_byyear_E001E002_2_fullN, ntrt==8 & exp==2)

#
#E001 #Non linear model fits##

#controlE001fit<- nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data = controlE001) # not working, singular gradient#
E001ntrt1fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt1) 
E001ntrt2fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt2) 
E001ntrt3fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt3) 
E001ntrt4fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt4) 
E001ntrt5fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt5) 
#E001ntrt6fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt6) ## not working, singular gradient#
E001ntrt7fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt7) 
E001ntrt8fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E001distntrt8) 


#E002 #Non linear model fits##

controlE002fit<- nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data = controlE002) 
E002ntrt1fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt1) 
E002ntrt2fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt2) 
E002ntrt3fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt3) 
E002ntrt4fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt4) 
E002ntrt5fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt5) 
E002ntrt6fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt6) 
E002ntrt7fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt7) 
#E002ntrt8fit<-nls(meantrajdist~ SSasymp(year2, Asym, R0, lrc), data =E002distntrt8) ### not working singular gradient ###


### for the ones that didn't work above specify the LRC at the estimate...###

x=controlE001$year2
y=controlE001$meantrajdist

summary(E001ntrt1fit)

##first, specify lrc from the closest model##
asympfunc <- function(x, Asym, R0, lrc=0.194)
  Asym + (R0 - Asym) * exp(-exp(lrc) * x)

controlE001fit<-nls(y ~ asympfunc(x, Asym, R0),
                    data = data.frame(x, y),
                    start = as.list(coef(E001ntrt1fit)[c(1, 2)])) ## this one works now ###

###
x=E001distntrt6$year2
y=E001distntrt6$meantrajdist

summary(E001ntrt5fit)


##first, specify lrc from the closest model##
asympfunc <- function(x, Asym, R0, lrc=0.00531)
  Asym + (R0 - Asym) * exp(-exp(lrc) * x)

E001ntrt6fit<-nls(y ~ asympfunc(x, Asym, R0),
                  data = data.frame(x, y),
                  start = as.list(coef(E001ntrt5fit)[c(1, 2)])) ## this one works now ###


###
x=E002distntrt8$year2
y=E002distntrt8$meantrajdist

summary(E002ntrt7fit)

##first, specify lrc from the closest model##
asympfunc <- function(x, Asym, R0, lrc=0.24417 )
  Asym + (R0 - Asym) * exp(-exp(lrc) * x)

E002ntrt8fit<-nls(y ~ asympfunc(x, Asym, R0),
                  data = data.frame(x, y),
                  start = as.list(coef(E002ntrt7fit)[c(1, 2)])) ## this one works now ###




###THESE ARE ALL THE VALUES IN TABLE S5###

summary(controlE001fit)
summary(E001ntrt1fit)
summary(E001ntrt2fit)
summary(E001ntrt3fit)
summary(E001ntrt4fit)
summary(E001ntrt5fit)
summary(E001ntrt6fit)
summary(E001ntrt7fit)
summary(E001ntrt8fit)


summary(controlE002fit)
summary(E002ntrt1fit)
summary(E002ntrt2fit)
summary(E002ntrt3fit)
summary(E002ntrt4fit)
summary(E002ntrt5fit)
summary(E002ntrt6fit)
summary(E002ntrt7fit)
summary(E002ntrt8fit)


### add the predicted values from fit to dataframes for plotting###

controlE001$prednls = predict(controlE001fit)
E001distntrt1$prednls = predict(E001ntrt1fit)
E001distntrt2$prednls = predict(E001ntrt2fit)
E001distntrt3$prednls = predict(E001ntrt3fit)
E001distntrt4$prednls = predict(E001ntrt4fit)
E001distntrt5$prednls = predict(E001ntrt5fit)
E001distntrt6$prednls = predict(E001ntrt6fit)
E001distntrt7$prednls = predict(E001ntrt7fit)
E001distntrt8$prednls = predict(E001ntrt8fit)


controlE002$prednls = predict(controlE002fit)
E002distntrt1$prednls = predict(E002ntrt1fit)
E002distntrt2$prednls = predict(E002ntrt2fit)
E002distntrt3$prednls = predict(E002ntrt3fit)
E002distntrt4$prednls = predict(E002ntrt4fit)
E002distntrt5$prednls = predict(E002ntrt5fit)
E002distntrt6$prednls = predict(E002ntrt6fit)
E002distntrt7$prednls = predict(E002ntrt7fit)
E002distntrt8$prednls = predict(E002ntrt8fit)


avgTraj_byyear_E001E002_2_fullN_2<-rbind(controlE001, E001distntrt1,E001distntrt2,E001distntrt3,E001distntrt4,E001distntrt5,E001distntrt6,E001distntrt7,E001distntrt8, controlE002,E002distntrt1,E002distntrt2,E002distntrt3,E002distntrt4,E002distntrt5,E002distntrt6,E002distntrt7,E002distntrt8 )


####################FIGURES #########################

### RECREATE FIG 5, Traj dist#####
avgTraj_byyear_E001E002_2_fullN_2

avgTraj_byyear_E001E002_2_fullN_2$ntrt2=avgTraj_byyear_E001E002_2_fullN_2$ntrt




###seperate the control from the others and calculate confidence intervals for plotting...###

# confidence intervals E001#

newdataE001 <- data.frame(
  x= controlE001$year2
)

newdataE002 <- data.frame(
  year2= controlE002$year2
)


confcontrolE001<-as.data.frame(predFit(controlE001fit, newdata = newdataE001, interval = "confidence", level=0.95))

confE001df<-cbind(controlE001, confcontrolE001)


confcontrolE002<-as.data.frame(predFit(controlE002fit, newdata = newdataE002, interval = "confidence", level=0.95))

confE002df<-cbind(controlE002, confcontrolE002)


confdistboth<-rbind(confE001df, confE002df)



# plot specifications

avgTraj_byyear_E001E002_2_fullN_2$exp<-mapvalues(avgTraj_byyear_E001E002_2_fullN_2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

avgTraj_byyear_E001E002_2_fullN_2$exp<- factor(avgTraj_byyear_E001E002_2_fullN_2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))



avgTraj_byyear_E001E002_2_fullN_2$ntrt2  <- mapvalues(avgTraj_byyear_E001E002_2_fullN_2$ntrt2, from=c("0", "9", "1", "2", "3", "4","5","6", "7", "8"),
                                                      to=c("Never", "None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

avgTraj_byyear_E001E002_2_fullN_2$ntrt2 <- factor(avgTraj_byyear_E001E002_2_fullN_2$ntrt2, levels = c("Never", "None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))





avgTraj_byyear_E001E002_2_fullN_2<-avgTraj_byyear_E001E002_2_fullN_2%>% mutate(abclabel=if_else(exp=="Intact in 1982 (E001)", "(a)", "(b)"))


confdistboth$ntrt2<-confdistboth$ntrt


confdistboth$ntrt2 <- mapvalues(confdistboth$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

confdistboth$ntrt2 <- factor(confdistboth$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))



confdistboth$exp<-mapvalues(confdistboth$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

confdistboth$exp<- factor(confdistboth$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))





####



mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))





avgTraj_byyearnonlinear_fullN<- ggplot(data=avgTraj_byyear_E001E002_2_fullN_2) +
  geom_point(data=avgTraj_byyear_E001E002_2_fullN_2, aes(x=as.numeric(as.character(year2)), y=meantrajdist, group=ntrt2, color=ntrt2)) +geom_line(data=avgTraj_byyear_E001E002_2_fullN_2, aes(x=as.numeric(as.character(year2)), y=prednls, group=ntrt2, color=ntrt2))+
  facet_grid(~exp)+ 
  geom_ribbon(data=confdistboth, aes(x=year2, ymin=lwr, ymax=upr,group=ntrt2, color=ntrt2), color=NA, fill = "gray1", alpha=0.2)+
  ylab(expression('Annual community trajectory distance'))+
  xlab("Time since experiment (years)")+
  mytheme+
  theme(legend.title = element_blank()) + 
  scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  geom_text(aes(x = 2.5, y = 0.65, label = abclabel, group =     abclabel),size = 4, color= "black", check_overlap = T)

