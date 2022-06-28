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

#start with exp12subset_1# 
#wide to long form#
exp12subset_long <- gather(exp12subset_1, species, biomass, 'mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)':'mass.above.VIOLA SP.', factor_key=TRUE)

#remove extra text#

exp12subset_long $species<-str_remove(exp12subset_long$species, "mass.above.")


#make a column for exp_ntrt_year#

exp12subset_long_2<-exp12subset_long  %>% mutate( expntrtyear_2= paste(exp,ntrt, year, sep = '_'))

exp12subset_long_3<-exp12subset_long_2 %>% mutate(expntrt= paste(exp,ntrt, sep = '_')) %>% dplyr::select(field, plot, disk,year, ntrt, nadd, expntrtyear_2, expntrt, species, biomass)

#sum all species abundance across plots with the same treatment group and year #includes zeroes!!!#

sum_allsp_exp12ABC<-exp12subset_long_3%>% group_by(species, disk,year, ntrt, nadd, expntrtyear_2, expntrt,) %>% dplyr::summarise(sumbiomass=sum(biomass, na.rm=TRUE),n=n())


#get rid of years where we don't have all 18 plots#

sum_allsp_exp12ABC_2<-subset(sum_allsp_exp12ABC, n==18)

sum_allsp_df<-as.data.frame(sum_allsp_exp12ABC_2)


#this dataframe (sumallsp_df has the abundance of each plant in each treatment group (n = 18) for each year#

#calculating the total biomass#

totalbiomassplot<-exp12subset_long_3%>% group_by(plot,field, disk,year, ntrt, nadd, expntrtyear_2, expntrt,) %>% dplyr::summarise(totalbiomass=sum(biomass, na.rm=TRUE),sebiomass=sd(biomass, na.rm=TRUE)/sqrt(n()),n=n())

#average across plots#

totalbiomassgroups<-totalbiomassplot%>% group_by(disk,year, ntrt, nadd, expntrtyear_2, expntrt,) %>% dplyr::summarise(totalplotbiomass=mean(totalbiomass, na.rm=TRUE),seplotbiomass=sd(totalbiomass, na.rm=TRUE)/sqrt(n()), n=n())

totalbiomassgroups2<-subset(totalbiomassgroups, n==18)

totalbiomassdf<-as.data.frame(totalbiomassgroups2)

totalbiomassdf2<-totalbiomassdf%>%mutate(upper=totalplotbiomass+(2*seplotbiomass), lower=totalplotbiomass-(2*seplotbiomass))%>% dplyr::select(expntrtyear_2,year, disk, ntrt, totalplotbiomass, seplotbiomass,upper, lower)

#this dataframe (totalbiomass df2, has the total biomass of all 18 plots per treatment group for all  years#)

#summary across treatments all years#

biomasssummary<-totalbiomassdf%>% group_by(disk) %>% dplyr::summarise(meanbiomass=mean(totalplotbiomass, na.rm=TRUE),sebiomass=sd(totalplotbiomass, na.rm=TRUE)/sqrt(n()),n=n())




#### RECREATE FIGURE S1, BIOMASS OVER TIME SERIES#####

#plot specification#

totalbiomassdf2$disk<-mapvalues(totalbiomassdf2$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

totalbiomassdf2$disk<- factor(totalbiomassdf2$disk, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

totalbiomassdf22<-totalbiomassdf2 %>% mutate(ntrt2=ntrt)

totalbiomassdf22$ntrt2  <- mapvalues(totalbiomassdf22$ntrt2, from=c("9", "1", "2", "4","6"), to=c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

totalbiomassdf22$ntrt2 <- factor(totalbiomassdf22$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

totalbiomassdf3<- totalbiomassdf22 %>% mutate(year2=as.numeric(as.character(year))-1981)


#plotting#

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))

totalbiomass<-totalbiomassdf3 %>% ggplot(aes(x=year2, y=totalplotbiomass, color=ntrt2))+
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=ntrt2), alpha=0.5, colour = NA) + 
  facet_grid(ntrt2~disk) + 
  mytheme+
  ylab(expression(paste("Total aboveground biomass g / m  "^2)))+
  xlab("Time since experiment")+
  theme(legend.title = element_blank())+
  theme(legend.key.size = unit(0.5, "cm"))+theme(strip.background = element_blank())+theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  scale_fill_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+ theme(legend.position = "none") 



#### CALUCULATE BIOMASS OF FUNCTIONAL GROUPS OVER TIME SERIES####

#start with the sum_all_sp_df and totabiomassdf2 dataframes and the original Da (all plot and spp info from CC)

#merge the spp with total plot biomass#

sum_allsp_df_2 <- (merge(sum_allsp_df,totalbiomassdf2, by = 'expntrtyear_2', sort=FALSE, all.x=TRUE))

sum_allsp_df_3<-sum_allsp_df_2 %>% mutate(propbiomass=sumbiomass/totalplotbiomass)

#this data set (sum_allsp_df_3 has the sum of all spp in each year in each dist/ntrt treatmet and the proportion of biomass of that group#

#merge the species info with functional group and other information#

specieslookup<-da %>% distinct(species, .keep_all = TRUE)

specieslookup2<-specieslookup %>% dplyr::select(species,functional.group,duration,lifeform, pathway, origin)

mergedfuncgroup <- (merge(specieslookup2,sum_allsp_df_3, by = 'species', sort=FALSE, all.y=TRUE))

mergedfuncgroup_2<-mergedfuncgroup%>% mutate(functional.group2= paste(functional.group, origin, sep = '_'))

#now sum up by functional group + origin #

sumfunctionalgroup<-mergedfuncgroup_2%>% group_by(disk,year, ntrt, nadd, expntrtyear_2, expntrt,functional.group, totalplotbiomass) %>% dplyr::summarise(sumbiomass=sum(sumbiomass, na.rm=TRUE)) %>% mutate(propbiomass=sumbiomass/totalplotbiomass)

sumfuntionalgroupdf<-as.data.frame(sumfunctionalgroup)

#reordering so that other is last, first grasses, then forbs#

sumfuntionalgroupdf$functional.group<-factor(sumfuntionalgroupdf$functional.group, levels=c("C3", "C4" , "S", "W", "L", "F"))

sumfuntionalgroupdf2<-sumfuntionalgroupdf[complete.cases(sumfuntionalgroupdf), ]

#drop the levels not sampeled
sumfuntionalgroupdf2$functional.group <- as.factor(as.character(sumfuntionalgroupdf2$functional.group))



#### RECREATE FIGURE S2, BIOMASS OF FUNCTION GROUPS OVER TIME####

#plot specifications#

sumfuntionalgroupdf2$disk<-mapvalues(sumfuntionalgroupdf2$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

sumfuntionalgroupdf2$disk<- factor(sumfuntionalgroupdf2$disk, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


sumfuntionalgroupdf2<-sumfuntionalgroupdf2 %>% mutate(ntrt2=ntrt)

sumfuntionalgroupdf2$ntrt2  <- mapvalues(sumfuntionalgroupdf2$ntrt2, from=c("9", "1", "2", "4","6"), to=c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

sumfuntionalgroupdf2$ntrt2 <- factor(sumfuntionalgroupdf2$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))


sumfuntionalgroupdf2$year<-as.numeric(as.character(sumfuntionalgroupdf2$year))


sumfuntionalgroupdf3<- sumfuntionalgroupdf2 %>% mutate(year2=year-1981)


propfunctionalgrouptime<-sumfuntionalgroupdf3 %>% ggplot(aes(x=year2, y=propbiomass, group=functional.group))+
  geom_area(position = 'stack', aes(fill=functional.group)) + 
  facet_grid(ntrt2~disk) + 
  mytheme+
  ylab("Proportion of Biomass")+
  xlab("Time since experiment")+
  theme(legend.title = element_blank())+
  theme(legend.key.size = unit(0.5, "cm"))+scale_fill_viridis_d(option="C", na.value="grey72")+
  theme(strip.background = element_blank())+theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_y_continuous(breaks=c(0.25,0.50,0.75, 1.00))


#### CALCULATE SPECIES RICHNESS OVER TIME SERIES####

#Diversity metrics for plant data#
ShannonD_a <- diversity(exp12subset_1[,11:222], index = "shannon") ##shannon-weiner index##
Richness_a <- specnumber(exp12subset_1[,11:222])
Evenness_a <- ShannonD_a/log(Richness_a)

##combine#
plotdatarichness<-cbind(exp12subset_1[,1:10], Richness_a)

plotdatarichnes2<- plotdatarichness%>%   mutate(ntrt2=ntrt,expntrtyear= paste(exp, ntrt,year, sep = '_'))

plotdatarichnes2$year<-as.numeric(as.character(plotdatarichnes2$year))

#average across treatments#

richnessaverage<-plotdatarichnes2%>% group_by(exp, ntrt, ntrt2, year, expntrtyear) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())

#richness at beginning#

richnessbeginning<-subset(plotdatarichnes2, year < 1985)

richnessbegsummary<-richnessbeginning%>% group_by(exp) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())

#richness at end of experiment#
richnessend<-subset(plotdatarichnes2, year > 2000)

richnessendsummary<-richnessend%>% group_by(exp) %>% dplyr::summarise(meanrichness=mean(Richness_a, na.rm=TRUE),sdrichness=sd(Richness_a, na.rm=TRUE),serichness=sd(Richness_a, na.rm=TRUE)/sqrt(n()), n=n())





#### RECREATE FIGURE S3, SPECIES RICHNESS OVER TIME SERIES#####

#plotting specifications#

richnessaverage$exp<-mapvalues(richnessaverage$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

richnessaverage$exp<- factor(richnessaverage$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


richnessaverage_2<- richnessaverage%>% mutate (year2=as.numeric(as.character(year-1981)))


richnessaverage_3<-richnessaverage_2 %>% mutate(ntrt2=ntrt)

richnessaverage_3$ntrt2  <- mapvalues(richnessaverage_3$ntrt2, from=c("9", "1", "2", "4","6"), to=c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

richnessaverage_3$ntrt2 <- factor(richnessaverage_3$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "3.4 N + micro","9 N + micro"))

#plot#

spprichnesstime<-richnessaverage_3 %>% ggplot(aes(x=year2, y=meanrichness, group=ntrt2, color=ntrt2, fill=ntrt2))+
  geom_point()+
  geom_smooth(method="loess")+
  facet_grid(ntrt2~exp) + 
  mytheme+
  ylab("Plot species richness")+
  xlab("Time since experiment")+
  mytheme+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  scale_fill_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  theme(strip.background = element_blank())