# MH DeSiervo, LG Shoemaker
# Disturbance alters transience but nutrients determine equilibria during grassland succession with multiple global change drivers submitted with DOI https://doi.org/10.5061/dryad.dbrv15f5t.

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


######################################################################################
##### Calculate the distance TO centroid within ntrt treatments####

###need to have the most recent verion of vegan from github for this to work##
#https://github.com/vegandevs/vegan/issues/331#

####datafrane with All nutrient treatments....####

exp12subset<-subset(da.wide5,year<2005)

##get rid of years not all fields were sampeled###

#exp12subset_full<-exp12subset%>%filter(year!=1995)%>%filter(year!=1998)%>%filter(year!=2001)%>%filter(year!=2003)

##just the veg data##

veg_full<-exp12subset%>% dplyr::select(`mass.above.ACER NEGUNDO`:`mass.above.VIOLA SP.`)



###log transformation####

vegmatrixfull<-as.matrix(veg_full)

logvegfull<-log(1+vegmatrixfull)



##using the apply function to calculation the dispersion (beta disp) for each nutrient treatment for each year#

bd_allbothfull<-betadisper(vegdist(logvegfull),exp12subset$expntrtyear, type = "centroid")    ####Note that you have to specify the centroid ##  ##takes a bit##

#bd_allbothEUC<-betadisper(vegdist(logveg, method="euclidean"),exp12subset_2$expntrtyear, type = "centroid")    ##same thing but in Euclidean distance###

##some squared distances are negative and changed to zero##


######making group distances into a dataframe with rownames###

bd_allbothfull_2<-as.data.frame(bd_allbothfull$group.distances, row.names = rownames(bd_allbothfull$group.distances))

setDT(bd_allbothfull_2, keep.rownames = TRUE)[]


bd_allbothfull_df<-separate(data = bd_allbothfull_2, col = rn, into = c("exp", "ntrt", "year"))

bd_allbothfull_df$year<-as.numeric(bd_allbothfull_df$year)
bd_allbothfull_df$ntrt<-as.factor(bd_allbothfull_df$ntrt)

bd_allbothfull_df$distances<-as.numeric(bd_allbothfull_df$`bd_allbothfull$group.distances`)

bd_allbothfull_df2<-bd_allbothfull_df%>%mutate(ntrt2 = ntrt)%>%mutate(year2 = year-1981)

###this dataframe (bd_allboth_fulldf2) has the average distance to the centroid for each disurbance_nutrient group for each year. Distances = the average Bray dist for 18 plots to their group centroids (what is plotted in Fig 3 a) ##

## INCLUDES ALL NUTRIENT GROUPS###


#### Model parameters dist TO centroid, RECREATE TABLE S1 in MS ####

head(bd_allbothfull_df2)

bd_allbothfull_df2<-bd_allbothfull_df2%>% mutate (year2=year-1981)


##null##

nulllin<-lm(distances~1, data=bd_allbothfull_df2)


###linear fit###

distlin<-lm(distances~year2, data=bd_allbothfull_df2)


##saturating function##
distnonlin<- nls(distances~ SSasymp(year2, Asym, R0, lrc), data = bd_allbothfull_df2)


###quadratic function###

distquad<-lm(distances~year2+I(year2^2), data=bd_allbothfull_df2)



rawaic<-AIC(nulllin, distlin, distnonlin, distquad)
nR<-dim(bd_allbothfull_df2)[1]  #Sample size 
aictable(rawaic,nR)

###QUADRACTIC is the best fit (now that we have the highest nutrient treatments###

##### parameters with error for table S2 ####

bd_allbothfull_df2<-bd_allbothfull_df2%>% mutate (year2=year-1981) %>% mutate(expntrt= paste(exp,ntrt, sep = '_'))



####
quadraticfits <- lmList(distances ~ year2+ I(year2^2) | expntrt, data=bd_allbothfull_df2)

quadraticfits_df <- data.frame(Subject=rownames(coef(fm1)),coef(fm1),check.names=FALSE)

######################################################################################

#### Calculate distance to nutrient centroid seperated by field ####

###need to have the most recent verion of vegan from github for this to work##
#https://github.com/vegandevs/vegan/issues/331#

##Used the datasubsetted to before 2004 ##

head(da.wide5)

####All nutrient treatments....####

exp12subset<-subset(da.wide5,year<2005)


##plot data##
plot_full<-exp12subset%>%dplyr::select(field:expntrtyear) %>% mutate(expfieldplotyear= paste(exp,field,plot,year, sep = '_'))

##just the veg data##

veg_full<-exp12subset%>% dplyr::select(`mass.above.ACER NEGUNDO` :`mass.above.VIOLA SP.`)



###log transformation####

vegmatrixfull<-as.matrix(veg_full)

logvegfull<-log(1+vegmatrixfull)



##using the apply function to calculation the dispersion (beta disp) for each nutrient treatment for each FIELD for each year#

bd_allbyfield<-betadisper(vegdist(logvegfull),plot_full$expntrtfieldyear, type = "centroid")    ####Note that you have to specify the centroid ##  ##takes a VERY LONG TIME ##

#bd_allbothEUC<-betadisper(vegdist(logveg, method="euclidean"),exp12subset_2$expntrtyear, type = "centroid")    ##same thing but in Euclidean distance###

##some squared distances are negative and changed to zero##


######making group distances into a dataframe with rownames###


bd_allbyfield_2<-as.data.frame(bd_allbyfield$group.distances, row.names = rownames(bd_allbyfield$group.distances))

setDT(bd_allbyfield_2, keep.rownames = TRUE)[]


bd_allbyfield_2_df<-separate(data = bd_allbyfield_2, col = rn, into = c("exp", "field", "ntrt", "year"))

bd_allbyfield_2_df$year<-as.numeric(bd_allbyfield_2_df$year)
bd_allbyfield_2_df$ntrt<-as.factor(bd_allbyfield_2_df$ntrt)

bd_allbyfield_2_df$distances<-as.numeric(bd_allbyfield_2_df$`bd_allbyfield$group.distances`)

###this dataframe (bd_allbyfield_2_df) has the average distance to the centroid for each disurbance_nutrient group in each field for each year. Distances = the average Bray dist for 6 plots to their group centroids (what is plotted in Fig 3 a) ##

######################################################################################

#### Calculate the distance BETWEEN nutrient group centroids##########


bd_allbothfull_df2<-bd_allbothfull_df%>%mutate(ntrt2 = ntrt)


##pulling out the group centroids by year####

centroiddf_full<-as.data.frame(bd_allbothfull[["centroids"]])

centroiddf_full_2<-centroiddf_full[,1:2]

centroiddf_full_2<- tibble::rownames_to_column(centroiddf_full_2, "expntrtyear")

centroiddf_full_2<-separate(data = centroiddf_full_2, col = expntrtyear, into = c("exp", "ntrt", "year"))

centroiddf_full_3<-centroiddf_full_2%>% mutate(expyear= paste(exp,year, sep = '_')) 


###grand mean centroid for all treatments each year###

grandcentroid_full<-centroiddf_full_3%>% group_by(exp, year) %>% dplyr::summarise(PcoA1centroid=mean(PCoA1, na.rm=TRUE),PcoA2centroid=mean(PCoA2, na.rm=TRUE))%>% mutate(expyear= paste(exp,year, sep = '_')) 


### merge the treatments means with the grand mean centroid for each year##

centroiddf_full<-merge(data.frame(centroiddf_full_3), data.frame(grandcentroid_full), by = "expyear", all = TRUE)

###calculate the distance from each centroid to grand centroid using pythagoreum theorem###

centroiddf_full_4<-centroiddf_full %>% mutate(distgrandcent=sqrt((PcoA1centroid- PCoA1)^2+(PcoA2centroid-PCoA2)^2))


###average distance of each group to its grand mean##

centroiddists_full<-centroiddf_full_4%>% group_by(expyear) %>% dplyr::summarise(meandist=mean(distgrandcent, na.rm=TRUE))%>%separate(col = expyear, into = c("exp", "year"))


#### AIC table distance BETWEEN centroids, RECREATE TABLE S3 in MS####

centroiddists_full2<-centroiddists_full%>% mutate (year2=as.numeric(as.character(year))-1981)

##null##
distnull<-lm(meandist~1, data=centroiddists_full2)



###linear fit###

distbtwnlin<-lm(meandist~year2, data=centroiddists_full2)


##saturating function##
distbtwnnonlin<- nls(meandist~ SSasymp(year2, Asym, R0, lrc), data = centroiddists_full2)

##quadratic##


distbtwnquad<-lm(meandist~year2+I(year2^2), data=centroiddists_full2)


rawaic<-AIC(distnull, distbtwnlin, distbtwnnonlin,distbtwnquad)
nR<-dim(bd_allboth_df2)[1]  #Sample size 
aictable(rawaic,nR)

###Nonlinear is the best fit###

######################################################################################

####### distance between disturbance group centroids ####

bd_allbothfull_df2<-bd_allbothfull_df%>%mutate(ntrt2 = ntrt)


##pulling out the group centroids by year####

centroiddf_full<-as.data.frame(bd_allbothfull[["centroids"]])

centroiddf_full_2<-centroiddf_full[,1:2]

centroiddf_full_2<- tibble::rownames_to_column(centroiddf_full_2, "expntrtyear")

centroiddf_full_2<-separate(data = centroiddf_full_2, col = expntrtyear, into = c("exp", "ntrt", "year"))

centroiddf_full_3<-centroiddf_full_2%>% mutate(ntrtyear= paste(ntrt,year, sep = '_')) 



###grand mean centroid E001 and E002 for all treatments each year###

grandcentroid_full_e001e002<-centroiddf_full_3%>% group_by(ntrt, year) %>% dplyr::summarise(PcoA1centroid=mean(PCoA1, na.rm=TRUE),PcoA2centroid=mean(PCoA2, na.rm=TRUE)) %>% mutate(ntrtyear= paste(ntrt,year, sep = '_'))


### merge the treatments means with the grand mean centroid for each year##

centroiddf_full_e001e002<-merge(data.frame(centroiddf_full_3), data.frame(grandcentroid_full_e001e002), by = "ntrtyear", all = TRUE)



###calculate the distance from each centroid to grand centroid using pythagoreum theorem###

centroiddf_full_4_e001e002<-centroiddf_full_e001e002 %>% mutate(distgrandcent=sqrt((PcoA1centroid- PCoA1)^2+(PcoA2centroid-PCoA2)^2))


###average distance of each group to its grand mean##

centroiddists_full_e001e002<-centroiddf_full_4_e001e002%>% group_by(ntrtyear) %>% dplyr::summarise(meandist=mean(distgrandcent, na.rm=TRUE))%>%separate(col = ntrtyear, into = c("ntrt", "year"))

######################################################################################

#######distance between disturbance group centroids by field ########

##START W/ DATA SET subsetted to 2004 to compare E001 and E002 ##

##subset to years before 2005 for BOTH E001 and E002 ## (bc of change in fire regime)

da.wide5$year<-as.numeric(as.character(da.wide5$year)) 

exp12subset<-subset(da.wide5,year<2005)


##just the veg data##

veg_full<-exp12subset%>% dplyr::select(`mass.above.ACER NEGUNDO` :`mass.above.VIOLA SP.`)


###log transformation####

vegmatrixfull<-as.matrix(veg_full)

logvegfull<-log(1+vegmatrixfull)



#### add on a exp field year column###


exp12subset_111 <- exp12subset %>%  mutate(expfieldyear= paste(exp, field, year, sep = '_'))

##using the apply function to calculation the dispersion (beta disp) for each nutrient treatment for each FIELD for each year#

bd_E001E002byfield<-betadisper(vegdist(logvegfull),exp12subset_111$expfieldyear, type = "centroid")    ####Note that you have to specify the centroid ##  ##takes a minute ##


##some squared distances are negative and changed to zero##


######making group distances into a dataframe with rownames###


bd_E001E002byfield_2<-as.data.frame(bd_E001E002byfield$group.distances, row.names = rownames(bd_E001E002byfield$group.distances))

setDT(bd_E001E002byfield_2, keep.rownames = TRUE)[]


bd_E001E002byfield_2_df<-separate(data = bd_E001E002byfield_2, col = rn, into = c("exp", "field", "year"))

bd_E001E002byfield_2_df$year<-as.numeric(bd_E001E002byfield_2_df$year)

bd_E001E002byfield_2_df$distances<-as.numeric(bd_E001E002byfield_2_df$`bd_E001E002byfield$group.distances`)

bd_E001E002byfield_2_df<-bd_E001E002byfield_2_df%>% mutate(year2=year-1981)

###this dataframe (bd_E001E002byfield_2_df) has the average distance to the centroid for each field in each year to its E001 or E002 group centroid. Distances = the average Bray dist for ALL 54 plots to their group centroids  ##



##pulling out the group centroids by year####

centroidE001E002df<-as.data.frame(bd_E001E002byfield[["centroids"]])

centroidE001E002df_2<-centroidE001E002df[,1:2]

centroidE001E002df_2 <- tibble::rownames_to_column(centroidE001E002df_2, "expfieldyear")

centroidE001E002df_2<-separate(data = centroidE001E002df_2, col = expfieldyear, into = c("exp", "field", "year"))

centroidE001E002df_2<-centroidE001E002df_2 %>%   mutate(fieldyear= paste(field, year, sep = '_')) %>% select(exp, field, year, fieldyear, PCoA1, PCoA2)
###

###calculate distance between E001 and E002 centroids###

centroidE001<-subset(centroidE001E002df_2, exp==1) #dim = 69 X 6 ##

centroidE002<-subset(centroidE001E002df_2, exp==2) #dim = 58 X 6 ##   ### three missing in E002####

##merge them...dropping the ones in E001 that we don't have E002 for####

mergedcentroids<-merge(centroidE002, centroidE001, by = "fieldyear", all.x = TRUE)

#use pythag theorem to get distance between E001 and E002 centroids###

mergedcentroids_dist<-mergedcentroids %>% mutate(distcentroids=sqrt((PCoA1.y - PCoA1.x)^2+(PCoA2.y-PCoA2.x)^2))

######################################################################################


################################# FIGURES ###############################

#### recreate figure 3AB distance to the nutrient centroid over time ###

head(bd_allbothfull_df2)


bd_allbothfull_df2<-bd_allbothfull_df2 %>% mutate(abclabel=exp)


bd_allbothfull_df2$abclabel<-mapvalues(bd_allbothfull_df2$abclabel, from=c("1", "2"), to=c("(a)", "(b)"))



##### plot specifications###

bd_allbothfull_df2$exp<-mapvalues(bd_allbothfull_df2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

bd_allbothfull_df2$exp<- factor(bd_allbothfull_df2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))



bd_allbothfull_df2$ntrt2  <- mapvalues(bd_allbothfull_df2$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                       to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

bd_allbothfull_df2$ntrt2 <- factor(bd_allbothfull_df2$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))


##just take out control for figure###


bd_allbothfull_control=subset(bd_allbothfull_df2, ntrt==9)



mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title.x=element_text(size=12))+theme(axis.title.y=element_text(size=10))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))


disttocentroidplot<- ggplot(bd_allbothfull_df2, aes(x=year2, y=distances, group=ntrt2, color=ntrt2)) +
  geom_point() +
  facet_grid(~exp)+ 
  ylab(expression('Average distance to nutrient centroid'))+
  xlab("Time since experiment")+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 0.5, se=FALSE)+
  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+
  stat_smooth(data=bd_allbothfull_control,aes(x=year2, y=distances), method = "lm", formula = y ~ x + I(x^2), size = 0.5)+ guides(color=guide_legend(override.aes=list(fill=NA)))+ geom_text(aes(x = 2, y = 0.6, label = abclabel, group = abclabel),size = 4, color= "black", check_overlap = T)


######################################################################################

#### RECREATE FIGURE 3 cd, dist BETWEEN centroids#####


centroiddists_full$year<-as.numeric(as.character(centroiddists_full$year))

centroiddists_full_2<-centroiddists_full%>% mutate(year2=year-1981)%>%mutate(abclabel=exp)



centroiddists_full_2$abclabel<-mapvalues(centroiddists_full_2$abclabel, from=c("1", "2"), to=c("(c)", "(d)"))

####subset by exp###

centroiddists_full_E001<-subset(centroiddists_full_2, exp==1)
centroiddists_full_E002<-subset(centroiddists_full_2, exp==2)




###https://www.r-bloggers.com/2021/07/asymptotic-confidence-intervals-for-nls-regression-in-r/##

asympfitcentroidE001_full<- nls(meandist~ SSasymp(year2, Asym, R0, lrc), data = centroiddists_full_E001)

asympfitcentroidE002_full <- nls(meandist~ SSasymp(year2, Asym, R0, lrc), data = centroiddists_full_E002)



###confidence intervals around fit###


new.dataE001 <- data.frame(year2=centroiddists_full_E001$year2)
new.dataE002 <- data.frame(year2=centroiddists_full_E002$year2)


confE001_centroid <- as_tibble(predFit(asympfitcentroidE001_full, newdata = new.dataE001, interval = "confidence", level= 0.95)) %>% mutate(year2 = new.dataE001$year2)

confE001_centroiddf<-as.data.frame(confE001_centroid)

confE001_centroiddf2<-cbind(centroiddists_full_E001, confE001_centroiddf)


confE002_centroid <- as_tibble(predFit(asympfitcentroidE002_full, newdata = new.dataE002, interval = "confidence", level= 0.95)) %>% mutate(year2 = new.dataE002$year2)


confE002_centroiddf<-as.data.frame(confE002_centroid)

confE002_centroiddf2<-cbind(centroiddists_full_E002, confE002_centroiddf)


##combine togeter and add abc label##
confcentroidboth<-rbind(confE001_centroiddf2, confE002_centroiddf2)

confcentroidboth<-confcentroidboth[,1:8]


###for plotting####


centroiddists_full_2$exp<-mapvalues(centroiddists_full_2$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

centroiddists_full_2$exp<- factor(centroiddists_full_2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

confcentroidboth$exp<-mapvalues(confcentroidboth$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

confcentroidboth$exp <- factor(centroiddists_full_2$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

distbetweencentroids_fullN<-ggplot(data=centroiddists_full_2,aes(x=as.numeric(as.character(year2)),y=meandist)) +
  geom_point() +
  facet_grid(~exp) +
  ylab(expression('Average distance between nutrient centroids'))+
  xlab("Time since experiment (years)")+
  geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), colour = "black", size = 0.5, se=F, fullrange=T)+
  geom_ribbon(data=confcentroidboth, aes(ymin=lwr, ymax=upr), color=NA, fill = "gray1", alpha=0.2)+
  mytheme+
  theme(legend.position="none")+
  theme(strip.text = element_blank())+
  theme(strip.background = element_blank())+
  mytheme+ geom_text(aes(x = 2, y = 0.22, label = abclabel, group = abclabel),size = 4, color= "black", check_overlap = T)+theme(strip.text = element_blank())



######################################################################################


################### SUPPLEMENTAL FIGURES ###################

## distance to centroid seperated by field ###

bd_allbyfield_2_df<-bd_allbyfield_2_df %>% mutate(abclabel=exp)%>% mutate(ntrt2=ntrt)%>% mutate(year2=year-1981)

bd_allbyfield_2_df$abclabel<-mapvalues(bd_allbyfield_2_df$abclabel, from=c("1", "2"), to=c("(a)", "(b)"))

##### plot specifications###

bd_allbyfield_2_df$exp<-mapvalues(bd_allbyfield_2_df$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

bd_allbyfield_2_df$exp<- factor(bd_allbyfield_2_df$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))



bd_allbyfield_2_df$ntrt2  <- mapvalues(bd_allbyfield_2_df$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                       to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

bd_allbyfield_2_df$ntrt2 <- factor(bd_allbyfield_2_df$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))


bd_allbyfield_2_df$field  <- mapvalues(bd_allbyfield_2_df$field, from=c("A", "B", "C"),
                                       to=c("Field A", "Field B", "Field C"))


mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title.x=element_text(size=12))+theme(axis.title.y=element_text(size=10))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))


disttocentroidplotBYFIELD<- ggplot(bd_allbyfield_2_df, aes(x=year2, y=distances, group=ntrt2, color=ntrt2)) +
  geom_smooth(method="lm", se=FALSE, size=0.5) +
  geom_point()  + 
  facet_grid(vars(field), vars(exp))+ 
  ylab(expression('Average distance to nutrient group centroid'))+
  xlab("Time since experiment")+
  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4))



disttocentroidplotBYFIELDloess<- ggplot(bd_allbyfield_2_df, aes(x=year2, y=distances, group=ntrt2, color=ntrt2)) +
  geom_smooth(method="loess", se=FALSE, size=0.5) +
  geom_point()  + 
  facet_grid(vars(field), vars(exp))+ 
  ylab(expression('Average distance to nutrient group centroid'))+
  xlab("Time since experiment")+
  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4))

######################################################################################

###### distance between nutrient centroids, seperated by field #####

##pulling out the group centroids by year####

centroidbyfield<-as.data.frame(bd_allbyfield[["centroids"]])

centroiddfbyfield_2<-centroidbyfield[,1:2]

centroiddfbyfield_2<- tibble::rownames_to_column(centroiddfbyfield_2, "expfieldntrtyear")

centroiddfbyfield_2<-separate(data = centroiddfbyfield_2, col = expfieldntrtyear, into = c("exp", "field", "ntrt", "year"))

centroiddfbyfield_3<-centroiddfbyfield_2%>% mutate(expfieldyear= paste(exp,field,year, sep = '_')) 
####



###grand mean centroid for all fields and treatments each year###

grandcentroidbyfield_full<-centroiddfbyfield_3%>% group_by(exp, field,year) %>% dplyr::summarise(PcoA1centroid=mean(PCoA1, na.rm=TRUE),PcoA2centroid=mean(PCoA2, na.rm=TRUE))%>% mutate(expfieldyear= paste(exp,field,year, sep = '_')) 


### merge the treatments means with the grand mean centroid for each year##

centroiddfbyfield_full<-merge(data.frame(centroiddfbyfield_3), data.frame(grandcentroidbyfield_full), by = "expfieldyear", all = TRUE)

###calculate the distance from each centroid to grand centroid using pythagoreum theorem###

centroiddfbyfield_full_4<-centroiddfbyfield_full %>% mutate(distgrandcent=sqrt((PcoA1centroid- PCoA1)^2+(PcoA2centroid-PCoA2)^2))


###average distance of each group to its grand mean##

centroiddistsbyfield_full<-centroiddfbyfield_full_4%>% group_by(expfieldyear) %>% dplyr::summarise(meandist=mean(distgrandcent, na.rm=TRUE))%>%separate(col = expfieldyear, into = c("exp","field", "year"))



centroiddistsbyfield_full$exp<-mapvalues(centroiddistsbyfield_full$exp, from=c("1", "2"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

centroiddistsbyfield_full$exp<- factor(centroiddistsbyfield_full$exp, levels = c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

centroiddistsbyfield_full_2<-centroiddistsbyfield_full%>% mutate(year2=as.numeric(year)-1981)


centroiddistsbyfield_full_2$field  <- mapvalues(centroiddistsbyfield_full_2$field, from=c("A", "B", "C"),
                                                to=c("Field A", "Field B", "Field C"))

distbetweencentroidsbyfield_fullN<-ggplot(data=centroiddistsbyfield_full_2,aes(x=as.numeric(as.character(year2)),y=meandist)) +
  geom_point() +
  facet_grid(vars(field), vars(exp))+ 
  ylab(expression('Average distance between nutrient group centroids'))+
  xlab("Time since experiment (years)")+
  geom_smooth(method="loess", colour = "black", size = 0.5, se=F, fullrange=T)+
  theme(legend.position="none")+
  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())

######################################################################################

#### supplemental figure, distance between disturbance group centroids ####

centroiddists_e001e002<-as.data.frame(centroiddists_full_e001e002)


centroiddists_e001e002_2<- centroiddists_e001e002 %>% mutate(year2=as.numeric(year)-1981)%>% mutate(ntrt2=ntrt)



centroiddists_e001e002_2$ntrt2  <- mapvalues(centroiddists_e001e002_2$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                             to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

centroiddists_e001e002_2$ntrt2 <- factor(centroiddists_e001e002_2$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))



distbetweencentroidse001e002<-ggplot(data=centroiddists_e001e002_2,aes(x=year2 ,y=meandist, group=ntrt2, color=ntrt2)) +
  geom_point() +
  ylab(expression('Average distance between \n disturbance group centroids'))+
  xlab("Time since experiment (years)")+
  geom_smooth(method="loess", size = 0.5, se=F, fullrange=T)+  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank()) + scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+theme(plot.margin = margin(10, 10, 10, 30))

######################################################################################
### supplemental figure, distance between disturbance group centroids, seperated by year ###

mergedcentroids_dist_2 <-mergedcentroids_dist%>% mutate(year2=as.numeric(year.x)-1981)

## for plotting##

mergedcentroids_dist_2$field.x<-mapvalues(mergedcentroids_dist_2$field.x, from=c("A", "B", "C"), to=c("Field A", "Field B", "Field C"))

distE001E002BYFIELD<- ggplot(mergedcentroids_dist_2, aes(x=year2, y=distcentroids)) +
  geom_point()  +
  geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), colour = "black", size = 0.5, se=F, fullrange=T)+
  facet_wrap(~field.x)+ 
  ylab(expression('Average distance between \n disturbace group centroids'))+
  xlab("Time since experiment")+
  mytheme+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) + 
  theme(strip.text = element_text(size=12, face="bold"))+
  theme(strip.background = element_blank())+
  theme(plot.margin = margin(10, 10, 10, 30))