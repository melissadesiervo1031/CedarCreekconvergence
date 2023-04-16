# MH DeSiervo, LG Shoemaker
# Disturbance alters transience but nutrients determine equilibria during grassland succession with multiple global change drivers ####submitted with DOI https://doi.org/10.5061/dryad.dbrv15f5t. 

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


#### CALCULATE directionality trajectories over whole time series ####

#https://cran.rstudio.com/web/packages/trajr/vignettes/trajr-vignette.html#
##START W/ the DISTANCE MATRIX instead of PCOA points ### 

head(BCdistalldf2)


##make each field, experiment, plot, its own combo###
plotE001E002<-BCdistalldf2 %>% select(field, exp, plot, subplot, year, disk, ntrt, expntrtyear)%>% mutate(fieldexpplot= paste(field, exp,plot, sep = '_'))

##
plotE001E002_2<-distinct(plotE001E002, fieldexpplot, .keep_all=TRUE)


trajdirection<-trajectoryDirectionality(BCdistall, plotE001E002$fieldexpplot, plotE001E002$year)


trajdirectdf<-cbind(plotE001E002_2, trajdirection)


##mean directionality for each exp_ntrt ###

trajdiravg<-trajdirectdf%>% group_by(exp, ntrt) %>% dplyr::summarise(trajdiravg=mean(trajdirection, na.rm=TRUE),sdtrajdir=sd(trajdirection, na.rm=TRUE),setrajdir=sd(trajdirection, na.rm=TRUE)/sqrt(n()), uprconfint=trajdiravg+(1.96*setrajdir),lwrconfint=trajdiravg-(1.96*setrajdir), n=n())

####################################################################
#### CALCULATE directionality across decades######

##arrange the dataframe so that that columns and rows are in right order##
exp12subsetsorted<-exp12subset%>%arrange(year)%>%arrange(plot)%>%arrange(field)%>%arrange(exp)


exp12subsetsorted$year<-as.numeric(as.character(exp12subsetsorted$year))

##just the plot data##
plotfirstdecade<-exp12subsetsorted%>% dplyr::select('field':'ntrt2') %>% filter(year<=1992)

plotseconddecade<-exp12subsetsorted%>% dplyr::select('field':'ntrt2') %>% filter(year>1992)


##just the veg data##

vegfirstdecade<-exp12subsetsorted%>% filter(year<=1992)%>% dplyr::select(`mass.above.ACER NEGUNDO`:`mass.above.VIOLA SP.`)

vegseconddecade<-exp12subsetsorted%>% filter(year>1992)%>% dplyr::select(`mass.above.ACER NEGUNDO`:`mass.above.VIOLA SP.`)


###log transformation####

vegmatrixfirst<-as.matrix(vegfirstdecade)

logvegfirst<-log(1+vegmatrixfirst)


vegmatrixsecond<-as.matrix(vegseconddecade)

logvegsecond<-log(1+vegmatrixsecond)


## Bray curtis dissimiliarty matrix for each decade##
distE1andE2first<-vegdist(logvegfirst, method="bray")
BCdistallfirst<-as.matrix(distE1andE2first)  ##into a matrix#
BCdistalldffirst<-as.data.frame(BCdistallfirst)  #into a dataframe#

BCdistalldffirst2<-cbind(plotfirstdecade,BCdistalldffirst) #merge plot data w/ dissim matrix##


distE1andE2second<-vegdist(logvegsecond, method="bray")
BCdistallsecond<-as.matrix(distE1andE2second)  ##into a matrix#
BCdistalldfsecond<-as.data.frame(BCdistallsecond)  #into a dataframe#

BCdistalldfsecond2<-cbind(plotseconddecade,BCdistalldfsecond) #merge plot data w/ dissim matrix##




##make each field, experiment, plot, its own combo###
plotE001E002first<-BCdistalldffirst2%>% select(field, exp, plot, subplot, year, disk, ntrt, expntrtyear)%>% mutate(fieldexpplot= paste(field, exp,plot, sep = '_'))

plotE001E002first_2<-distinct(plotE001E002first, fieldexpplot, .keep_all=TRUE)

plotE001E002second<-BCdistalldfsecond2%>% select(field, exp, plot, subplot, year, disk, ntrt, expntrtyear)%>% mutate(fieldexpplot= paste(field, exp,plot, sep = '_'))

plotE001E002second_2<-distinct(plotE001E002second, fieldexpplot, .keep_all=TRUE)


#### traj direction on each decade ####

trajdirectionfirst<-trajectoryDirectionality(BCdistallfirst, plotE001E002first$fieldexpplot, plotE001E002first$year)

trajdirectionsecond<-trajectoryDirectionality(BCdistallsecond, plotE001E002second$fieldexpplot, plotE001E002second$year)

####

trajdirectdffirst<-cbind(plotE001E002first_2, trajdirect=trajdirectionfirst)

trajdirectdfsecond<-cbind(plotE001E002second_2, trajdirect=trajdirectionsecond)

###

### add column for decades###

trajdirectdffirst_2<-trajdirectdffirst %>% mutate(decade="1982-1992")

trajdirectdfsecond_2<-trajdirectdfsecond %>% mutate(decade="1993-2004")

##combine them together##
trajectoriesboth<-rbind(trajdirectdffirst_2, trajdirectdfsecond_2)

##mean traj direction for each exp_ntrt###

trajdirectavgbydecade<-trajectoriesboth%>% group_by(exp, ntrt, decade) %>% dplyr::summarise(trajdiravg=mean(trajdirect, na.rm=TRUE),sdtrajdir=sd(trajdirect, na.rm=TRUE),setrajdir=sd(trajdirect, na.rm=TRUE)/sqrt(n()), n=n())

##########################################################################
#################FIGURES############################

###### RECREATE FIGURE 4 ab, directionality over time series ##########

##for plotting###
trajdiravg2<-trajdiravg %>% mutate(ntrt2=ntrt)


trajdiravg2$ntrt2  <- mapvalues(trajdiravg2$ntrt2, from=c("0", "9", "1", "2", "3", "4","5","6", "7", "8"),
                                to=c("Never", "None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

trajdiravg2$ntrt2 <- factor(trajdiravg2$ntrt2, levels = c("Never", "None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))


trajdiravg2$exp<-mapvalues(trajdiravg2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

trajdiravg2$exp<- factor(trajdiravg2$exp, levels = c("Remnant", "Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

trajdiravg2<-trajdiravg2 %>% mutate(abclabel=if_else(exp=="Intact in 1982 (E001)", "(a)", "(b)"))


mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title.x=element_text(size=12))+theme(axis.title.y=element_text(size=10))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))


directionalityplot<- ggplot(data=trajdiravg2, aes(x=ntrt2, y=trajdiravg, group=ntrt2, color=ntrt2)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lwrconfint, ymax=uprconfint), size=1)+
  facet_grid(~exp, scales = "free_x")+ 
  xlab("")+
  ylab(expression(atop("Average directionality", paste("(1982 - 2004)"))))+
  mytheme+
  theme(legend.position="none")+
  theme(legend.title = element_blank()) + 
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.text.x = element_text(angle = 25,vjust = 0.5))+
  theme(strip.text = element_text(face="bold"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+scale_y_continuous(labels = scales::number_format(accuracy = 0.005))+geom_text(aes(x = 0.8, y = 0.49, label = abclabel, group =     abclabel),size = 4, color= "black", check_overlap = T)

##### RECREATE FIGURE S4, directionality across decades ####

trajdirectavgbydecade<-trajectoriesboth%>% group_by(exp, ntrt, decade) %>% dplyr::summarise(trajdiravg=mean(trajdirect, na.rm=TRUE),sdtrajdir=sd(trajdirect, na.rm=TRUE),setrajdir=sd(trajdirect, na.rm=TRUE)/sqrt(n()), n=n())




##for plotting###
trajdirectavgbydecade2<-trajdirectavgbydecade %>% mutate(ntrt2=ntrt)

## plot specifications

trajdirectavgbydecade2$ntrt2 <- mapvalues(trajdirectavgbydecade2$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                          to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

trajdirectavgbydecade2$ntrt2<- factor(trajdirectavgbydecade2$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))



trajdirectavgbydecade2$exp<-mapvalues(trajdirectavgbydecade2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

trajdirectavgbydecade2$exp<- factor(trajdirectavgbydecade2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))



mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title.x=element_text(size=12))+theme(axis.title.y=element_text(size=10))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))



trajdirplotbydecade<- ggplot(data=trajdirectavgbydecade2, aes(x=ntrt2, y=trajdiravg, group=ntrt2, color=ntrt2)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=trajdiravg-setrajdir, ymax=trajdiravg+setrajdir), size=0.5)+
  facet_grid(decade~exp)+ 
  xlab("")+
  ylab(expression('Average directionality'))+
  mytheme+
  theme(legend.position="none")+
  theme(legend.title = element_blank()) + 
  theme(legend.title = element_blank()) + 
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.text.x = element_text(angle = 25,vjust = 0.5))+
  theme(strip.text = element_text(face="bold"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))

