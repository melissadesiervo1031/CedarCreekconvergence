# MH DeSiervo, LG Shoemaker
### Disturbance alters transience but nutrients determine equilibria during grassland succession with multiple global change drivers ######submitted with DOI https://doi.org/10.5061/dryad.dbrv15f5t. 
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

source(here("code/1_functions.r"))
source(here("code/2_upload_format_data.r"))
source(here("code/3_PCoA_Permanova.r"))


### INDICATOR SPECIES ANALYSIS#####
### 2 time periods (early late) and including field D ###

##starting with exp12subset, raw biomass data####


head(da_fieldD_3 ) #long form field D##
da_fieldD_3$exp<-0

head(d3woody)  ##long form fields ABC##

#######

d3E001E002<-d3woody %>% select(field, exp,ntrt, plot, year, species, biomass=mass.above)


da_fieldD_3$exp=0
da_fieldD_3$ntrt=0


da_fieldD_3<-da_fieldD_3 %>% select(field, exp,ntrt, plot, year, species, biomass=mass.above)


###
dim(d3E001E002)
dim(da_fieldD_3)

####

fieldABCDE001E002<-rbind(d3E001E002, da_fieldD_3)

###

###subset to the controls and high N treatments###

fieldABCDE001E002_subset<-subset(fieldABCDE001E002,ntrt==9|ntrt==0|ntrt==8) ###control treatment  (9) highest N 27.2 (8)



# Look at first 3 years of experiment and last 3 years of experiment
vegearlylong <- subset(fieldABCDE001E002_subset, year=="1982"|year=="1983"|year=="1984")
veglatelong <- subset(fieldABCDE001E002_subset, year=="2000"|year=="2002"|year=="2004")


##skip years 2001, 2003 bc of imcomplete data##

vegearlylong_2<-vegearly %>% mutate(stage="early") %>% mutate(stage2=1)

veglatelong_2<-veglatelong %>% mutate(stage="late")%>% mutate(stage2=3)

#####merge the early and late##

veglong2<-rbind(vegearlylong_2, veglatelong_2)


####this file for indicator species #####

dim(veglong2) ##dim 6316    9


#rownames(veglong2)<-NULL ##get rid of rownames##


##indicator species analysis##

##convert data back to wide form##

vegwide2 <- spread(veglong2, species, biomass) ##dim 648 X 161####


# Fill in NA's with zeros in wide data set
vegwide2 [is.na(vegwide2)] <- 0

vegwideplotdata2<-vegwide2[,1:7]
vegwidematrix2<-vegwide2[,8:161]


###make a column for dist_nitrogen_earylate##

vegwideplotdata22<-vegwideplotdata2 %>% mutate(stagentrtdist= paste(stage, ntrt, exp, sep = '_'))

vegwideplotdata22$stagentrtdist<-as.factor(vegwideplotdata22$stagentrtdist)

###
levels(vegwideplotdata22$stagentrtdist)

vegwideplotdata22$stagentrtdist <- factor(vegwideplotdata22$stagentrtdist, levels = c("early_0_0", "early_9_1", "early_8_1", "early_9_2", "early_8_2", "late_0_0", "late_9_1", "late_8_1", "late_9_2", "late_8_2"))

vegwideplotdata33<-vegwideplotdata22 %>% mutate(stagentrtdist2= as.numeric(stagentrtdist))

vegwideplotdata33$stagentrtdist2<-as.factor(vegwideplotdata33$stagentrtdist2)

vegwideplotdata33$stagentrtdist2 <- factor(vegwideplotdata33$stagentrtdist2, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

##early  late ##

##indicator spp analyses###

indvalearlylatefieldABCD = multipatt(vegwidematrix2,vegwideplotdata33$stagentrtdist2,max.order = 3, control = how(nperm=99))## just looks at the individual groups, and twos  ##

indvalearlylatefieldABCD2 = multipatt(vegwidematrix2,vegwideplotdata33$stagentrtdist2,max.order = 4, control = how(nperm=999))## just looks at the individual groups, and twos and threes  ##


indval<-data.frame(indvalearlylatefieldABCD2[["sign"]]) ##RESULTS###

indval

