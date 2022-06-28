# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code to run indicator spp analyses ####

#Recreate Table 1#

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

### INDICATOR SPECIES ANALYSIS#####

#just the plot data#

plot<-exp12subset_1_sorted%>% dplyr::select('field':'ntrt2')

#just the veg data#

veg<-exp12subset_1_sorted%>% dplyr::select(`mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)`:`mass.above.VIOLA SP.`)

#log transformation#

vegmatrix<-as.matrix(veg)

logveg<-log(1+vegmatrix)

#merge back together###

exp12subset_3<-data.frame(plot, logveg)

head(exp12subset_3) ##log transformed values###

#subset to the controls and high N treatments#

exp12subset_4<-subset(exp12subset_3,ntrt==9|ntrt==6)


# Look at first 3 years of experiment and last 3 years of experiment
vegearly <- subset(exp12subset_4, year=="1982"|year=="1983"|year=="1984")
veglate <- subset(exp12subset_4, year=="2000"|year=="2002"|year=="2004")

#skip years 2001, 2003 bc of imcomplete data#

#Put data in long format
vegearlylong <-reshape(vegearly,
                       v.names="mass",
                       timevar="species",
                       idvar=c("field", "exp", "plot", "year", "ntrt", "nadd", "disk"),  
                       times=grep("mass.above", names(vegearly), value=TRUE),
                       varying=list(grep("mass.above", names(vegearly), value=TRUE)),
                       direction="long")

vegearlylong_2<-vegearlylong %>% mutate(stage="early") %>% mutate(stage2=1)


veglatelong <-reshape(veglate,
                      v.names="mass",
                      timevar="species",
                      idvar=c("field", "exp", "plot", "year", "ntrt", "nadd", "disk"),  
                      times=grep("mass.above", names(veglate), value=TRUE),
                      varying=list(grep("mass.above", names(veglate), value=TRUE)),
                      direction="long")

veglatelong_2<-veglatelong %>% mutate(stage="late")%>% mutate(stage2=3)

#merge the early and late##

veglong2<-rbind(vegearlylong_2, veglatelong_2)

# Delete prefix from species names
veglong2$species <- with(veglong2, gsub("mass.above.","", species))

rownames(veglong2)<-NULL ##get rid of rownames##

#convert data back to wide form##

vegwide2 <- spread(veglong2, species, mass)

vegwideplotdata2<-vegwide2[,1:13]
vegwidematrix2<-vegwide2[,14:224]

#make a column for dist_nitrogen_earylate#

vegwideplotdata22<-vegwideplotdata2 %>% mutate(stagentrtdist= paste(stage, ntrt, exp, sep = '_'))

vegwideplotdata22$stagentrtdist<-as.factor(vegwideplotdata22$stagentrtdist)

vegwideplotdata22$stagentrtdist <- factor(vegwideplotdata22$stagentrtdist, levels = c("early_9_1", "early_6_1", "early_9_2", "early_6_2", "late_9_1", "late_6_1", "late_9_2", "late_6_2"))

vegwideplotdata33<-vegwideplotdata22 %>% mutate(stagentrtdist2= as.numeric(stagentrtdist))

vegwideplotdata33$stagentrtdist2<-as.factor(vegwideplotdata33$stagentrtdist2)

vegwideplotdata33$stagentrtdist2 <- factor(vegwideplotdata33$stagentrtdist2, levels = c("1", "2", "3", "4", "5", "6", "7", "8"))

#early  late ##

#indicator spp analyses#

indvalearlylate = multipatt(vegwidematrix2,vegwideplotdata33$stagentrtdist2,max.order = 3, control = how(nperm=99))## just looks at the individual groups, and twos  ##


indval<-data.frame(indvalearlylate[["sign"]]) ##RESULTS###

