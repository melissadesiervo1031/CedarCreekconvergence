# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code to diversity metrics and biomass ####

#Recreate Figure S1, S2, S3 in MS#

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
source(here("code/6_trajdist.r"))
source(here("code/7_Indicatorspp.r"))

#### CALCULATE BIOMASS OVER TIME SERIES ####
##starting with exp12subset, raw biomass data####


dim(exp12subset) ##6102 by 198, 1-10 plot data, the rest veg data###


##wide to long form##
exp12subset_long <- gather(exp12subset, species, biomass, 'mass.above.ACER NEGUNDO':'mass.above.VIOLA SP.', factor_key=TRUE)

##remove extra text##

exp12subset_long $species<-str_remove(exp12subset_long$species, "mass.above.")


##make a column for exp_ntrt_year##

exp12subset_long_2<-exp12subset_long  %>% mutate( expntrtyear_2= paste(exp,ntrt, year, sep = '_'))

exp12subset_long_3<-exp12subset_long_2 %>% mutate(expntrt= paste(exp,ntrt, sep = '_')) %>% dplyr::select(field, plot, disk,year, ntrt,  expntrtyear_2, expntrt, species, biomass)

##sum all species abundance across plots with the same treatment group and year ###includes zeroes!!!#####

#sum_allsp_exp12ABC<-exp12subset_long_3%>% group_by(species, disk,year, ntrt, expntrtyear_2, expntrt,) %>% dplyr::summarise(sumbiomass=sum(biomass, na.rm=TRUE),n=n())


###get rid of years where we don't have all 18 plots###

#sum_allsp_exp12ABC_2<-subset(sum_allsp_exp12ABC, n==18)

#sum_allsp_df<-as.data.frame(sum_allsp_exp12ABC_2)

### due to differing number of plots through time, use average instead###


mean_allsp_exp12ABC<-exp12subset_long_3%>% group_by(species, disk,year, ntrt, expntrtyear_2, expntrt,) %>% dplyr::summarise(meanbiomass=mean(biomass, na.rm=TRUE),n=n())



##this dataframe (sumallsp_df has the abundance of each plant in each treatment group (n = 18) for each year####

##calculating the total biomass##

totalbiomassplot<-exp12subset_long_3%>% group_by(plot,field, disk,year, ntrt, expntrtyear_2, expntrt,) %>% dplyr::summarise(totalbiomass=sum(biomass, na.rm=TRUE),sebiomass=sd(biomass, na.rm=TRUE)/sqrt(n()),n=n())

##average across plots##

totalbiomassgroups<-totalbiomassplot%>% group_by(disk,year, ntrt,  expntrtyear_2, expntrt,) %>% dplyr::summarise(totalplotbiomass=mean(totalbiomass, na.rm=TRUE),seplotbiomass=sd(totalbiomass, na.rm=TRUE)/sqrt(n()), n=n())


totalbiomassdf<-as.data.frame(totalbiomassgroups)

totalbiomassdf2<-totalbiomassdf%>%mutate(upper=totalplotbiomass+(2*seplotbiomass), lower=totalplotbiomass-(2*seplotbiomass))%>% dplyr::select(expntrtyear_2,year, disk, ntrt, totalplotbiomass, seplotbiomass,upper, lower)

####summary across treatments all years###

biomasssummary<-totalbiomassdf%>% group_by(disk) %>% dplyr::summarise(meanbiomass=mean(totalplotbiomass, na.rm=TRUE),sebiomass=sd(totalplotbiomass, na.rm=TRUE)/sqrt(n()),n=n())

############################################################

#### CALUCULATE BIOMASS OF FUNCTIONAL GROUPS OVER TIME SERIES####

##start with the mean_allsp_exp12ABC and totalbiomassgroups dataframes ####

## and the original Da (all plot and spp info from CC)

totalbiomassplot2<-as.data.frame(totalbiomassgroups) %>% select(expntrtyear_2, totalplotbiomass, seplotbiomass)


##merge the spp with total plot biomass##

sum_allsp_df_2 <- (merge(mean_allsp_exp12ABC,totalbiomassplot2, by = 'expntrtyear_2', sort=FALSE, all.x=TRUE))

sum_allsp_df_3<-sum_allsp_df_2 %>% mutate(propbiomass=meanbiomass/totalplotbiomass)

### this data set (sum_allsp_df_3 has the sum of all spp in each year in each plot and the proportion of biomass of that sp##

###merge the species info with functional group and other information##

specieslookup<-sp.df %>% distinct(species, .keep_all = TRUE)

specieslookup2<-specieslookup %>% dplyr::select(species,functional.group,duration,lifeform, pathway, origin)


mergedfuncgroup <- (merge(specieslookup2,sum_allsp_df_3, by = 'species', sort=FALSE, all.y=TRUE))

mergedfuncgroup_2<-mergedfuncgroup%>% mutate(functional.group2= paste(functional.group, origin, sep = '_'))

####now sum up by functional group + origin ##

sumfunctionalgroup<-mergedfuncgroup_2%>% group_by(disk,year, ntrt,expntrtyear_2, expntrt,functional.group, totalplotbiomass) %>% dplyr::summarise(sumfuncbiomass=sum(meanbiomass, na.rm=TRUE))%>% mutate(propbiomass=sumfuncbiomass/totalplotbiomass)

sumfunctionalgroupdf<-as.data.frame(sumfunctionalgroup)


###check that the proportions add to 1###

checksum<-sumfunctionalgroup%>% group_by(expntrtyear_2) %>% dplyr::summarise(sumpropbiomass=sum(propbiomass, na.rm=TRUE))



#################################################################

#### Calculate species ricchness over time series ####

##starting with exp12subset, raw biomass data####

dim(exp12subset) ##6102 by 198, 1-10 plot data, the rest veg data###

###summarize##

###Diversity metrics for plant data###
ShannonD_a <- diversity(exp12subset[,11:198], index = "shannon") ##shannon-weiner index##
Richness_a <- specnumber(exp12subset[,11:198])
Evenness_a <- ShannonD_a/log(Richness_a)

###combine##
plotdatarichness<-cbind(exp12subset[,1:10], Richness_a)

plotdatarichnes2<- plotdatarichness%>%   mutate(ntrt2=ntrt,expntrtyear= paste(exp, ntrt,year, sep = '_'))

##make year numeric##

plotdatarichnes2$year<-as.numeric(as.character(plotdatarichnes2$year))

##average across treatments##

richnessaverage<-plotdatarichnes2%>% group_by(exp, ntrt, ntrt2, year, expntrtyear) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())


##average across fields####

richnessaverageABC<-plotdatarichnes2%>% group_by(exp,field, year, ntrt, ntrt2, year, expntrtyear) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())

richnessABC19822004<-subset(richnessaverageABC,year==1982|year==2004)
richnessABC19822004<-subset(richnessaverageABC,year==1982|year==2004)


##richness at beginning###

richnessbeginning<-subset(plotdatarichnes2, year < 1985)

richnessbegsummary<-richnessbeginning%>% group_by(exp) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())

richnessbegsummaryboth<-richnessbeginning %>% group_by(ntrt)%>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())



richnessfirstyear<-subset(plotdatarichnes2, year == 1982)

richnessbegsummary<-richnessbeginning%>% group_by(exp) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())




##richness at end of experiment###
richnessend<-subset(plotdatarichnes2, year > 2000)

richnessendsummary<-richnessend%>% group_by(ntrt2) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())

######################################################################



##################FIGURES ############################


#### RECREATE FIGURE S1, BIOMASS OVER TIME SERIES#####


##plot specifications ##


totalbiomassdf2$disk<-mapvalues(totalbiomassdf2$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

totalbiomassdf2$disk<- factor(totalbiomassdf2$disk, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


totalbiomassdf22<-totalbiomassdf2 %>% mutate(ntrt2=ntrt)


totalbiomassdf22$ntrt2 <- mapvalues(totalbiomassdf22$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                    to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

totalbiomassdf22$ntrt2 <- factor(totalbiomassdf22$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

##



totalbiomassdf3<- totalbiomassdf22 %>% mutate(year2=as.numeric(as.character(year))-1981)


totalbiomassdf3subset<-subset(totalbiomassdf3, ntrt==9|ntrt==2|ntrt==5|ntrt==8)

###plotting###

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))

totalmeanbiomassplot<-totalbiomassdf3subset %>% ggplot(aes(x=year2, y=totalplotbiomass, color=ntrt2))+
  geom_line(size=1) + 
  #geom_ribbon(aes(ymin=lower, ymax=upper, fill=ntrt2), alpha=0.5, colour = NA) + 
  facet_grid(~disk) + 
  mytheme+
  ylab(expression(paste("Average aboveground biomass g / m  "^2)))+
  xlab("Time since experiment")+
  theme(legend.position = "right")+
  theme(strip.background = element_blank())+theme(strip.text.x = element_text(size = 14, colour = "black")) + scale_colour_manual(values = c("dark gray","#f98e09", "#8a226a","#000004"))+ scale_fill_manual(values = c("dark gray","#f98e09", "#8a226a","#000004"))


#######################################################################

#### RECREATE FIGURE S2, BIOMASS OF FUNCTION GROUPS OVER TIME####




###reordering so that other is last, first grasses, then forbs##

sumfunctionalgroupdf$functional.group<-factor(sumfunctionalgroupdf$functional.group, levels=c("C3", "C4" , "S", "W", "L", "F"))

#sumfuntionalgroupdf2<-sumfuntionalgroupdf[complete.cases(sumfuntionalgroupdf), ]

#drop the levels not sampeled
sumfunctionalgroupdf$functional.group <- as.factor(as.character(sumfunctionalgroupdf$functional.group))


##plot specifications####

sumfunctionalgroupdf$disk<-mapvalues(sumfunctionalgroupdf$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

sumfunctionalgroupdf$disk<- factor(sumfunctionalgroupdf$disk, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


sumfunctionalgroupdf$ntrt2<-sumfunctionalgroupdf$ntrt



sumfunctionalgroupdf$ntrt2 <- mapvalues(sumfunctionalgroupdf$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                        to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

sumfunctionalgroupdf$ntrt2 <- factor(sumfunctionalgroupdf$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))



sumfunctionalgroupdf$year<-as.numeric(as.character(sumfunctionalgroupdf$year))


sumfunctionalgroupdf2<- sumfunctionalgroupdf %>% mutate(year2=year-1981)



mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))


### subset for plotting###
sumfuntionalgroupdf3subset<-subset(sumfunctionalgroupdf2, ntrt==9|ntrt==2|ntrt==5|ntrt==8)



propfunctionalgrouptime<-sumfuntionalgroupdf3subset %>% ggplot(aes(x=year2, y=propbiomass, group=functional.group))+
  geom_area(position = 'stack', aes(fill=functional.group)) + 
  facet_grid(ntrt2~disk) + 
  mytheme+
  ylab("Proportion of Biomass")+
  xlab("Time since experiment")+
  theme(legend.title = element_blank())+
  theme(legend.key.size = unit(0.5, "cm"))+scale_fill_viridis_d(option = "D",direction=-1, na.value="grey72")+
  theme(strip.background = element_blank())+theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_y_continuous(breaks=c(0.25,0.50,0.75, 1.00))


##############################################################################

#### RECREATE FIGURE S3, SPECIES RICHNESS OVER TIME SERIES#####

##plotting specifications###

richnessaverage$exp<-mapvalues(richnessaverage$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

richnessaverage$exp<- factor(richnessaverage$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


richnessaverage$ntrt2 <- mapvalues(richnessaverage$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                   to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

richnessaverage$ntrt2 <- factor(richnessaverage$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))


richnessaverage$year2= richnessaverage$year - 1981



mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))




spprichnesstime<-richnessaverage %>% ggplot(aes(x=year2, y=meanrichness, group=ntrt2, color=ntrt2, fill=ntrt2))+
  geom_point()+
  geom_smooth(method="loess", se=FALSE)+
  facet_grid(~exp) + 
  mytheme+
  ylab("Plot species richness")+
  xlab("Time since experiment")+
  mytheme+
  theme(legend.title = element_blank()) + 
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+
  theme(strip.background = element_blank())
