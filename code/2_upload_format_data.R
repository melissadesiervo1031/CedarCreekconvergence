
# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code to upload the Cedar Creek biomass data and prep the data for analyses####


###PACKAGES####

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

############


# Reading in and cleaning Cedar Creek data###

#Read in e001 data## Plot data for the unplowed fields##

d1s <- read.csv(here("data/e001-soil-cn-2019-09-13.csv"), header=T, strip.white=T, skip=1)

# Delete some extra variables and add dummies to stack exp1 & exp2 data
d1s$natm.nadd <- NULL
d1s$fertadd <- NULL
d1s$burn.trt.1992 <- NA
d1s$ntrt.last.year <- 'CurrentYear'
d1s$ntrt.origin <- 1 # No cessation in E001
d1s$disk <- 0 # E001 was not disked before experiment (disturbed vegetation)
d1s$exp <- 1

#Read in e002 data ##Plot data for the plowed fields##

d2s <- read.csv(here("data/e002-soil-cn-2019-09-13.csv"), header=T, strip.white=T, skip=1)

## Delete some extra variables and add dummies to stack exp1 & exp2 data
d2s$X <- NULL
d2s$date <- NULL
d2s$BurnRecBeginning1992 <- NULL
d2s$NtrtRecBefore1992 <- NULL
d2s$TrtRecAfter1991 <- NULL
d2s$TrtRecAfter2013 <- NULL
d2s$burn <- NA
d2s$fence <- NA
d2s$fence.origin <- NA
d2s$depth <- NA
d2s$disk <- 1 # E002 was disked before experiment (disturbed vegetation)
d2s$exp <- 2


#Additioanl plot info for both experiments##

d3s <- read.csv(here("data/E001-E002-Soil-CN-2018.csv"), header=T, strip.white=T, skip=0, na.strings=c(""))


# Delete samples that have errors
d3s <- d3s[is.na(d3s$FLAGGED.FOR.RERUN..S.Spillage..E.Error.),]

# Delete some extra variables and add dummies to stack exp1 & exp2 data
d3s$FLAGGED.FOR.RERUN..S.Spillage..E.Error. <- NULL
d3s$Sample <- NULL
d3s$X <- NULL
d3s$burn <- NA
d3s$burn.origin <- NA
d3s$burn.trt.1992 <- NA
d3s$fence <- NA
d3s$fence.origin <- NA
d3s$depth <- NA
d3s$disk <- NA
d3s$nadd <- NA
d3s$ntrt <- NA
d3s$ntrt.last.year <- NA
d3s$ntrt.origin <- NA
d3s$trt.origin <- NA  

d3s$disk[d3s$exp==1] <- 0 # E001 was not disked before experiment (intact vegetation)
d3s$disk[d3s$exp==2] <- 1 # E002 was disked before experiment (disturbed vegetation)

# Check that data sets have the same variables
sort(names(d1s))
sort(names(d2s))
sort(names(d3s))

# Stack exp1 & 2 data
ds <- rbind(d1s,d2s,d3s)


#Read in e001 aboveground plant biomass data ##

d1a <- read.csv(here("data/e001-aboveground-mass-2019-09-13.csv"), header=T, strip.white=T, skip=1)



# Delete some extra variables and add dummies to stack exp1 & exp2 data
d1a$natm.nadd <- NULL
d1a$fertadd <- NULL
d1a$burn.trt.1992 <- NA
d1a$ntrt.last.year <- "CurrentYear" # No cessation in E001
d1a$ntrt.origin <- 1 # No cessation in E001
d1a$subplot <- NA
d1a$disk <- 0 # E001 was not disked before experiment (disturbed vegetation)
d1a$exp <- 1


#Read in e002 aboveground plant biomass data ##
d2a <- read.csv(here("data/e002-aboveground-mass-2019-09-13.csv"), header=T, strip.white=T, skip=1)


# Delete some extra variables and add dummies to stack exp1 & exp2 data
d2a$X <- NULL
d2a$date <- NULL
d2a$burn <- NA
d2a$fence <- NA
d2a$fence.origin <- NA
d2a$disk <- 1 # E002 was disked before experiment (intact vegetation)
d2a$exp <- 2


# Check that data sets have the same variables
sort(names(d1a))
sort(names(d2a))


# Stack exp1 & 2 data
da <- rbind(d1a,d2a)

# Delete records with key missing data
da <- da[!is.na(da$mass),]

# Disking not done if field D
da <- da[da$field != 'D',]
da$field <- as.factor(as.character(da$field))


# If subplot is missing specify it as a whole plot
da$subplot[is.na(da$subplot)] <- "Whole"

# Code other Nutrient additions
da$other.add <- 1
da$other.add[da$ntrt == 9] <- 0

# Capitalize species to get rid of capilization differences in spelling
da$species <- toupper(as.character(da$species))


# Capitalize species to get rid of capilization differences in spelling
da$species <- toupper(as.character(da$species))

# Don't include woody plants in biomass as there are some big trees in samples
da$live <- 1
da$sorted <- 1
da$wood <- 0
da$vasc <- 1

# Do some general substitutions

da$species <- gsub("APOCYNUM CANNABINUM", "APOCYNUM ANDROSAEMIFOLIUM", da$species) ## MD 11/1 based off email with Eric##
da$species <- gsub("MISC. FORB", "MISCELLANEOUS FORB", da$species)
da$species <- gsub("SEDGES", "CAREX SP.", da$species)
da$species <- gsub("QUERCUS RUBRUM", "QUERCUS RUBRA", da$species)

# Find litter and code as not alive
sel<-da$species == 'MISCELLANEOUS LITTER'
da$live[sel] <- 0

sel<-grep("PINE", da$species)
da$live[sel] <- 0

# Code unsorted material as not being sorted
sel<-grep("MISCELLANEOUS", da$species)
da$sorted[sel] <- 0

sel<-grep("FUNGI", da$species)
da$sorted[sel] <- 0

sel<-grep("MOSS", da$species)
da$sorted[sel] <- 0

sel<-grep("LICHEN", da$species)
da$sorted[sel] <- 0

# Woody stuff
sel<-grep("ACER NEGUNDO", da$species)
da$wood[sel] <- 1
sel<-grep("CEANOTHUS", da$species)
da$wood[sel] <- 1
sel<-grep("CORYLUS", da$species)
da$wood[sel] <- 1
sel<-grep("PINUS", da$species)
da$wood[sel] <- 1
sel<-grep("PINE", da$species)
da$wood[sel] <- 1
sel<-grep("POPULUS", da$species)
da$wood[sel] <- 1
sel<-grep("QUERCUS", da$species)
da$wood[sel] <- 1
sel<-grep("ULMUS", da$species)
da$wood[sel] <- 1

# Check for outlier masses
# High values are lorge amounts of litter or trees
da[da$mass.above > 1000 & da$live==1 & da$woo==0,c("year","field", "exp", "plot", "species", "mass.above" )]

# There are a couple species with more than 1500 g of mass that
# are very suspcious based on the species
da <- da[da$mass.above <= 1500, ]


# Read in species attribute look up table
sp.df <- read.csv(here("data/Cedar_Creek_Plant_Taxon_List.csv"),header=T)



names(sp.df)[] <- tolower(names(sp.df))

sp.df<-plyr::rename(sp.df, c("ï..species"="species"))
sp.df$species <- toupper(sp.df$species) 
sp.df$functional.group <- toupper(sp.df$functional.group)
sp.df$duration <- toupper(sp.df$duration)
sp.df$lifeform <- toupper(sp.df$lifeform)
sp.df$pathway <- toupper(sp.df$pathway)
sp.df$origin <- toupper(sp.df$origin)
sp.df$family <- toupper(sp.df$family)

da <- merge(da,sp.df, by='species', all.x=TRUE)

# Set Unknowns to missing
da$origin[da$origin == 'UNKNOWN'] <- NA
da$origin[da$origin == 'NATIVE AND/OR INTRODUCED'] <- NA


da$functional.group[da$functional.group == 'UNKNOW'] <- NA

da$duration[da$duration == 'UNKNOWN'] <- NA
da$duration[da$duration == 'BIENNIAL, PERENNIAL'] <- "BIENNIAL"

# Lump biennials in with annuals as they are pretty similar
# mostly weedy forbs 
unique(da$species[da$duration=="BIENNIAL"])
unique(da$species[da$duration=="BIENNIAL, PERENNIAL"])
unique(da$species[grep('BIENNIAL', da$duration)]) 

da$duration[grep('BIENNIAL', da$duration)] <- "ANNUAL"
unique(da$species[da$duration=="ANNUAL" & da$functional.group=="C4"])
unique(da$species[da$duration=="ANNUAL" & da$functional.group=="F"])

# Add in some species atrributes
sel <- da$species == "QUERCUS RUBRA"
da$functional.group[sel] <- "W"
da$lifeform[sel] <- "WOODY"
da$origin[sel] <- "NATIVE"

# Add in some species atrributes
sel <- da$species == "RUDBEKIA SEROTINA"
da$functional.group[sel] <- "F"
da$lifeform[sel] <- "FORB"
da$origin[sel] <- "NATIVE"

# Add in some species atrributes
sel <- da$species == "POA PRATENSIS"
da$origin[sel] <- "INTRODUCED"

sel <- da$species == "MISCELLANEOUS CAREX SP."
da$functional.group[sel] <- "S"
da$lifeform[sel] <- "SEDGE"
da$origin[sel] <- "NATIVE"


# Capitalize species to get rid of capilization differences in spelling
da$species <- toupper(as.character(da$species))

##deal with some entries where sp. were weighed twice##
summary(freq <- ddply(da[da$live==1 & da$sorted==1, ], .(year, field, exp, plot, subplot, disk, ntrt, nadd, species), colwise(length, .(mass.above)))) ##takes a while##

freq$freq <- freq$mass.above

freq$mass.above <- NULL
freq[freq$freq > 1,]

doubles.df <- merge(freq[freq$freq > 1,], da[c("field", "exp","plot", "year", "species", "mass.above")], by=c("field", "exp","plot", "year", "species"))
doubles.df[c("field", "exp","plot", "year", "species", "mass.above")]

da.mn <- ddply(da, .(field, exp, plot, subplot, year, disk, ntrt, nadd, species, live, sorted, wood, functional.group, lifeform, duration, origin), colwise(mean, .(mass.above)))
##takes a while##

# subset to live, sorted, herbaceous plants
d2 <- da.mn[da.mn$sorted ==1 & da.mn$live ==1 & da.mn$wood==0, c("field", "exp","plot", "subplot", "year", "disk", "ntrt", "nadd", "species", "mass.above")]


#Add some columns####
d3<-d2 %>% mutate(expntrtfieldyear= paste(exp,field, ntrt, year, sep = '_'), expntrtyear= paste(exp,ntrt, year, sep = '_'), ntrt2=ntrt) 

##This data frame (d3) has fields A, B, C, experiment 1 (intact) & exp 2 (disked in 1982) in a long data dataformat ##9 levels of nutrient addition#### years 1982 - 2019 (not every field / exp sampeled in every year## 


# Transpose data to be a site by species matrix
da.wide <- reshape(d3,
                   v.names="mass.above",
                   idvar=c("field", "exp","plot", "subplot", "year", "disk", "ntrt", "nadd", "expntrtyear"),
                   timevar="species",
                   direction="wide") 
# Fill in NA's with zeros in wide data set
da.wide[is.na(da.wide)] <- 0

##This data frame has fields A, B, C, experiment 1 (intact) & exp 2 (disked in 1982) ## years 1982 - 2019 for E002 and years 1982 - 2004 for exp E001 (not every field / exp sampeled in every year##
###wide version format###

# N Treatments (Ntrt) Details (from Table 1 Tilman 1987)
# Ntrt Trt g N/m2/yr Other nutrients 
# 1     A       0.0      All 
# 2     B       1.0      All 
# 3     C       2.0      All 
# 4     D       3.4      All 
# 5     E       5.4      All 
# 6     F       9.5      All 
# 7     G      17.0      All 
# 8     H      27.2      All 
# 9     I       0.0      None 

##make nutrient and years factors#

da.wide$ntrt<-as.factor(da.wide$ntrt)
da.wide$year<-as.factor(da.wide$year)

########Combine east / west subplots into "whole" to match up with rest of data#######

##make only the biomasses the numeric variables##
da.wide$exp<-as.factor(da.wide$exp)
da.wide$disk<-as.factor(da.wide$disk)
da.wide$nadd<-as.factor(as.character(da.wide$nadd))


#subset out the east/ west subplots plots and the whole plots ##

da.wideeastwest<-subset(da.wide, subplot=="East"|subplot=="West", row.names=NULL)
da.widewhole<-subset(da.wide, subplot=="Whole", row.names=NULL)

##Add the values from the subplots together##

EWaddall<-as.data.frame(da.wideeastwest%>%group_by(field,exp,plot,year,disk,ntrt,nadd,expntrtyear, expntrtfieldyear, ntrt2) %>% summarise_if(is.numeric,sum)%>% mutate(subplot="Whole"))


###merge the whole plots back together###
da.wide_allwhole<-rbind(da.widewhole,EWaddall)



dim(da.wide_allwhole) ## dimensions of the dataset = 8798 X 222##

##reorder columns so that all the rows are in the same order across fields etc. 
da.wide5<-da.wide_allwhole%>%arrange(year)%>%arrange(plot)%>%arrange(exp)%>%arrange(field) 


##### Check that the max number of rep within a treatment is 6, since 6 plots in each field and exp recieved the same N treatment##

max(ave(da.wide5$plot, da.wide5$expntrtyear,FUN = seq_along))





##subset to years before 2005## (bc of change in fire regime)

da.wide5$year<-as.numeric(as.character(da.wide5$year)) 

exp12subset<-subset(da.wide5,year<2005)


##subset to selected nutrient treatments below 9.5 g N##

exp12subset_1<-subset(exp12subset,ntrt==9|ntrt==1|ntrt==2|ntrt==4|ntrt==6)

### this data set (exp12subset_1 has fields ABC, E001 and E002 from 1982-2004) ## > 2004 was a change in the fire regime## ONLY the control and selected nutrient treatments below 9.5 g N) ##

###some notes about missing data##
## E001 data set is complete from 1982 - 2004. E002 dataset is missing years 1995, 1998, 2001, 2003 ALL fields##

