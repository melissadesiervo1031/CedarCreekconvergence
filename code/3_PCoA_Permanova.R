# MH DeSiervo, LG Shoemaker
### Disturbance alters transience but nutrients determine equilibria during grassland succession with multiple global change drivers ####### submitted with DOI https://doi.org/10.5061/dryad.dbrv15f5t. 
# Code to run a PCOA on dataset and conduct permanova####
##Recreate Figures 1, 2, 4 (cd) in MS####


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

source(here("code/1_functions.r"))
source(here("code/2_upload_format_data.r"))


############   CONDUCT PCOA fields ABC (using all the veg data across years) ##############

exp12subset<-subset(da.wide5,year<2005)

##make nutrient and years factors#

exp12subset$ntrt<-as.factor(exp12subset$ntrt)
exp12subset$year<-as.factor(exp12subset$year)


##arrange the dataframe so that that columns and rows are in right order##
exp12subsetsorted<-exp12subset%>%arrange(year)%>%arrange(plot)%>%arrange(field)%>%arrange(exp)


##just the plot data##
plot<-exp12subsetsorted%>% dplyr::select('field':'ntrt2')

##just the veg data##

veg<-exp12subsetsorted%>% dplyr::select(`mass.above.ACER NEGUNDO`:`mass.above.VIOLA SP.`)



###log transformation####

vegmatrix<-as.matrix(veg)

logveg<-log(1+vegmatrix)


## Bray curtis dissimiliarty matrix##
distE1andE2<-vegdist(logveg, method="bray")
BCdistall<-as.matrix(distE1andE2)  ##into a matrix#
BCdistalldf<-as.data.frame(BCdistall)  #into a dataframe#

BCdistalldf2<-cbind(plot,BCdistalldf) #merge plot data w/ dissim matrix##

##BCdistalldf2 has the plot info and the BC dissimilarity matrix for the entire experiment (E001 and E002, all three fields combined, with a subsetted selection of the ntrt treatments)

###PcOA for both E011 and E002 together####

PCOA_allfieldE12_fullN <- pcoa(distE1andE2)  ## takes a while##

PCOA_allfieldE12_cmd <- cmdscale(distE1andE2, eig=TRUE, add=FALSE)  ## Another way of doing the pcoa that gives same results#


## extracting the axis scores for each exp##

PCOAaxesboth_fullN <- PCOA_allfieldE12_fullN$vectors[,c(1,2)]


##merge the axis scores with the plot data ###

PCOAandplotsboth_fullN<-cbind(plot,PCOAaxesboth_fullN )

##PCOAandplotsboth_full N has the plot info and the PCOA scores for the entire experiment (E001 and E002, all three fields combined, wiah ALL of the ntrt treatments)


####species scores#####

PCOA_allfieldE12_cmd_wscores<-add.spec.scores(PCOA_allfieldE12_cmd, veg,  method='cor.scores')

spscores<-as.data.frame(PCOA_allfieldE12_cmd_wscores$cproj)
names<-rownames(spscores)
rownames(spscores) <- NULL

spscores_df <- cbind(names,spscores)

spscores_df$names <- with(spscores_df, gsub("mass.above.","", names))

############# CONDUCT Permanova: nutrient, disturbance, field effect #########

exp12subset<-subset(da.wide5,year<2005)


####for this analysis, getting rid of the years where not every field was sampeled##

exp12subset_22<-exp12subset%>%filter(year!=1995)%>%filter(year!=1998)%>%filter(year!=2001)%>%filter(year!=2003)



###log transform the biomass data###

justveg<-as.matrix(exp12subset_22[,12:198])

logveg<-log(1+justveg)

exp12subset_22<-cbind(exp12subset_22[,1:11], logveg)



### Extract two  years year to get the apply function to work#####

da.wide_1990<-subset(exp12subset_22, year=="1990", row.names=NULL)

da.wide_1999<-subset(exp12subset_22, year=="1999", row.names=NULL)

###
exp12subset_19821991<-subset(exp12subset_22, as.numeric(as.character(year)) < 1992)
exp12subset_19922004<-subset(exp12subset_22, as.numeric(as.character(year)) >= 1992)

#drop the levels not sampeled

levels(exp12subset_19821991$year)
exp12subset_19821991$year <- as.factor(as.character(exp12subset_19821991$year))


levels(exp12subset_19922004$year)
exp12subset_19922004$year <- as.factor(as.character(exp12subset_19922004$year))


#######################PERMANOVA####################################

### investigate the blocking / strata issue ###

## Correct with strata
test<-with(da.wide_1990, adonis2(vegdist(da.wide_1990[,11:198]) ~ exp+ntrt+field, data = da.wide_1990,  strata = da.wide_1990$field))

test2<-adonis2(vegdist(da.wide_1990[,11:198])~exp+ntrt+field, data= da.wide_1990)

#### using strata rescricts permutations to within fields. Makes no difference in Sum of Sqs or R2 (what we report) -- only matters for the p value ####


##Do a PERMANOVA looping over years, looking at effect of exp (disturbed or not), ntrt treatment# and field (A, B, C)###

##need to do this in two batches because from 1992 on there are fewer plots..####

###
permallfields1<-lapply(X = split.data.frame(x = exp12subset_19821991, f=exp12subset_19821991$year),function(t){
  z<-adonis2(vegdist(t[,11:198])~exp+ntrt+field, data= da.wide_1990)
  return(data.frame(z))
}
)



permallfields2<-lapply(X = split.data.frame(x = exp12subset_19922004, f=exp12subset_19922004$year),function(t){
  z<-adonis2(vegdist(t[,11:198])~exp+ntrt+field, data= da.wide_1999)
  return(data.frame(z))
}
)




##dataframe of results##
permmaster1<-plyr::ldply(permallfields1, data.frame)
permmaster2<-plyr::ldply(permallfields2, data.frame)


names<-c("Disturbance", "Fertilization", "Field", "residual", "total")

sourcevariation<-rep(names,length(unique(permmaster1$.id)))

permmaster11<-cbind(sourcevariation,permmaster1)

sourcevariation<-rep(names,length(unique(permmaster2$.id)))

permmaster22<-cbind(sourcevariation,permmaster2)

###

permmasterall<-rbind(permmaster11, permmaster22)

permmaster3_allnut<-subset(permmasterall,sourcevariation=="Disturbance"|sourcevariation=="Fertilization"|sourcevariation=="Field") 


############ CONDUCT PCOA w/ field d ##########################

head(da_fieldD_2 ) #long form field D##
da_fieldD_2$exp<-0

head(d3woody)  ##long form fields ABC##

#######

##just subset out control plots##
d3E001E002<-d3woody %>% filter(ntrt==9) %>% select(field, exp, plot, year, species, mass.above)


da_fieldD_3<-da_fieldD_2 %>% filter(ntrt==9)%>% select(field, exp, plot, year, species, mass.above)

fieldABCDE001E002<-rbind(da_fieldD_3, d3E001E002)

# Transpose data to be a site by species matrix
fieldABCD_wide <- reshape(fieldABCDE001E002,
                          v.names="mass.above",
                          idvar=c("field", "plot", "year"),
                          timevar="species",
                          direction="wide") 


fieldABCD_wide_sorted<-fieldABCD_wide%>%arrange(year)%>%arrange(plot)%>%arrange(field)%>%arrange(exp)

###count of plots###

##there should be 3 or 6 plots per field A, B, C to 2004, and 5 in field D ###########
plotcountcheck2<-fieldABCD_wide_sorted%>%group_by(exp, field, year) %>%  tally()



##subset to fewer years###

fieldABCD_wide_sorted_1982<-subset(fieldABCD_wide_sorted, year==1982)


##just the plot data field ABCD##
plotABCD1982<-fieldABCD_wide_sorted_1982%>% dplyr::select(field, exp, plot, year)


##just the veg data##

vegABCD1982<-fieldABCD_wide_sorted_1982[5:181]


# Fill in NA's with zeros in wide data set
vegABCD1982[is.na(vegABCD1982)] <- 0


vegmatrixABCD1982<-as.matrix(vegABCD1982)

logvegABCD1982<-log(1+vegmatrixABCD1982)


## Bray curtis dissimiliarty matrix##
distABCD1982<-vegdist(logvegABCD1982, method="bray")
BCdistABCD1982<-as.matrix(distABCD1982)  ##into a matrix#
BCdistABCDdf1982<-as.data.frame(BCdistABCD1982)  #into a dataframe#

BCdistABCDdf21982<-cbind(plotABCD1982,BCdistABCDdf1982) #merge plot data w/ dissim matrix##

##BCdistABCD2 has the plot info and the BC dissimilarity matrix for all the control plots from E001 ABC and field D ######

###PcOA fields ABCD####

PCOA_ABCD1982 <- cmdscale(BCdistABCD1982, eig=TRUE, add=FALSE) 


## extracting the axis scores for each exp##

PCOAaxes_ABCD1982<- PCOA_ABCD1982$points[,c(1,2)]


##merge the axis scores with the plot info...##

PCOAABCDandplots1982<-cbind(plotABCD1982, PCOAaxes_ABCD1982)

##rename columns##

names(PCOAABCDandplots1982)[5] <- "Axis.1"
names(PCOAABCDandplots1982)[6] <- "Axis.2"


PCOAABCDandplots1982$exp<-as.factor(as.character(PCOAABCDandplots1982$exp))
PCOAABCDandplots1982$field<-as.factor(PCOAABCDandplots1982$field)



#################FIGURES###############################



######### RECREATE FIGURE 1, variation over time  ##########



permmaster3_allnut2<- permmaster3_allnut %>% mutate (year=as.numeric(.id))  %>% mutate (year2=year-1981)

permmaster3_allnut2$sourcevariation<- factor(permmaster3_allnut2$sourcevariation, levels = c("Disturbance", "Fertilization", "Field"))

my_plotlabels = paste0('(', letters, ') ')

abclabels = paste0(my_plotlabels[1:length(levels(permmaster3_allnut2$sourcevariation))]) 

permmaster5<-cbind(permmaster3_allnut2, abclabels)


##make the plot##

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))


variancetimeall_allnut<-permmaster5 %>% ggplot(aes(x=as.numeric(year2), y=R2, group=sourcevariation)) +
  geom_point() +
  facet_wrap(~sourcevariation, strip.position = c("top"))+
  geom_smooth(method="loess", colour="black", size=0.5, fill = "gray1") +
  geom_text(aes(x = 2.0, y = 0.54, label = abclabels, group =     abclabels),size = 4,check_overlap = T)+
  ylab(expression('Explained variation'))+
  xlab("Time since experiment (years)")+
  mytheme+
  theme(strip.text = element_text(face="bold", size=12))+
  theme(strip.background = element_blank())

######### RECREATE FIGURE 2, PCOA ordination by year ########

PCOAandplotsboth_fullN



#first check that there are 18 OR 9 in each ntrt treatment / year##
PCOAandplotscount_fullN<-PCOAandplotsboth_fullN %>% count(year,disk,ntrt)

####2003 for E002 only has 1 field so we remove it###

PCOAandplotsboth2_fullN<-PCOAandplotsboth_fullN %>%filter(year!=2003)

PCOAandplotscount2_fullN<-PCOAandplotsboth2_fullN %>% count(year,disk,ntrt)##good now###

PCOAandplotsboth2_fullN$Axis.1<-as.numeric(PCOAandplotsboth2_fullN$Axis.1)
PCOAandplotsboth2_fullN$Axis.2<-as.numeric(PCOAandplotsboth2_fullN$Axis.2)

PCOAandplotsboth2_fullN$year<-as.numeric(as.character(PCOAandplotsboth2_fullN$year))



PCOAandplotsavg_fullN<-PCOAandplotsboth2_fullN %>% group_by(exp,year,disk,ntrt) %>% dplyr::summarise(meanAxis1=mean(Axis.1, na.rm=TRUE),sdAxis1=sd(Axis.1, na.rm=TRUE),seAxis1=sd(Axis.1, na.rm=TRUE)/sqrt(n()),meanAxis2=mean(Axis.2, na.rm=TRUE),sdAxis2=sd(Axis.2, na.rm=TRUE),seAxis2=sd(Axis.2, na.rm=TRUE)/sqrt(n()), n=n())



####

## plot specifications

PCOAandplotsavg2_fullN<-mutate(PCOAandplotsavg_fullN, ntrt2=ntrt)

PCOAandplotsavg2_fullN$ntrt2<-PCOAandplotsavg2_fullN$ntrt



PCOAandplotsavg2_fullN$ntrt2  <- mapvalues(PCOAandplotsavg2_fullN$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                           to=c("None", "0 N + μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9 N + μ", "17 N + μ", "27.2 N + μ"))

PCOAandplotsavg2_fullN$ntrt2 <- factor(PCOAandplotsavg2_fullN$ntrt2, levels = c("None", "0 N + μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9 N + μ", "17 N + μ", "27.2 N + μ"))




PCOAandplotsavg2_fullN<-PCOAandplotsavg2_fullN  %>% mutate(yearlabel=ifelse(year==1982|year==1992|year==2004,year, NA))%>%mutate(exp_ntrt= as.factor(paste(exp, ntrt, sep = '_')))

PCOAandplotsavg2_fullN$exp_ntrt <- factor(PCOAandplotsavg2_fullN$exp_ntrt , levels = c("1_9", "2_9", "1_1","2_1", "1_2", "2_2","1_3", "2_3", "1_4", "2_4","1_5", "2_5", "1_6", "2_6", "1_7", "2_7", "1_8", "2_8"))

PCOAandplotsavg3_fullN<-PCOAandplotsavg2_fullN %>% mutate(abclabel=exp_ntrt)
levels(PCOAandplotsavg3_fullN$abclabel) <- c( "(a) ", "(b) ", "(c) ", "(d) ", "(e) ", "(f) ", "(g) ", "(h) ", "(i) ", "(j) ", "(k) ", "(l) ", "(m) ", "(n) ", "(o) ", "(p) ", "(q) ", "(r) ")


PCOAandplotsavg3_fullN$disk<-plyr::mapvalues(PCOAandplotsavg3_fullN$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))





mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=7, colour = "black"))+theme(axis.text.y=element_text(size=7, colour = "black"))+theme(axis.title=element_text(size=10))+theme(plot.title=element_text(size=10) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85", size=0.25))+theme(legend.text = element_text(size=8))+theme(legend.key.size = unit(0.45, 'cm'))


##landscape version###

L_PcoAfieldAveragedE001E002_fullN<-ggplot(PCOAandplotsavg3_fullN[order(PCOAandplotsavg3_fullN$year),], aes(x=meanAxis1  , y=meanAxis2 , group=exp_ntrt, color=year, label=year)) +
  geom_point(size=0.35) +
  geom_errorbar(aes(ymin=meanAxis2-seAxis2, ymax=meanAxis2+seAxis2), color="dark gray", size=0.25) +
  geom_errorbarh(aes(xmin=meanAxis1-seAxis1, xmax=meanAxis1+seAxis1), color="dark gray", size=0.25) +
  geom_path(aes(x=meanAxis1 , y=meanAxis2 , group=exp_ntrt, color=year),size=0.25)+
  geom_point(size=0.75) +
  ylab("PCoA 2")+
  xlab("PCoA 1")+
  facet_grid(disk~ntrt2) + 
  mytheme+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D",trans = 'reverse')+
  guides(#reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE),
    #reverse size order (higher diameter on top) 
    size = guide_legend(reverse = TRUE))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(face="bold", size=10))+
  theme(panel.spacing.y=unit(0.5, "lines"))+
  scale_x_continuous(limits = c(-0.4, 0.47), breaks = c(-0.3, 0.0, 0.3))+
  scale_y_continuous(limits = c(-0.48, 0.27), breaks = c(-0.25, 0.0, 0.25))


##portrait version###

P_PcoAfieldAveragedE001E002_fullN<-ggplot(PCOAandplotsavg3_fullN[order(PCOAandplotsavg3_fullN$year),], aes(x=meanAxis1  , y=meanAxis2 , group=exp_ntrt, color=year, label=year)) +
  geom_point(size=0.35) +
  geom_errorbar(aes(ymin=meanAxis2-seAxis2, ymax=meanAxis2+seAxis2), color="dark gray", size=0.25) +
  geom_errorbarh(aes(xmin=meanAxis1-seAxis1, xmax=meanAxis1+seAxis1), color="dark gray", size=0.25) +
  geom_path(aes(x=meanAxis1 , y=meanAxis2 , group=exp_ntrt, color=year),size=0.25)+
  geom_point(size=0.75) +
  ylab("PCoA 2")+
  xlab("PCoA 1")+
  facet_grid(ntrt2~disk) + 
  mytheme+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D",trans = 'reverse')+
  guides(#reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE),
    #reverse size order (higher diameter on top) 
    size = guide_legend(reverse = TRUE))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(face="bold", size=12))+
  theme(panel.spacing.y=unit(0.5, "lines"))+
  scale_x_continuous(limits = c(-0.38, 0.48), breaks = c(-0.3, 0.0, 0.3))+
  scale_y_continuous(limits = c(-0.48, 0.27), breaks = c(-0.25, 0.0, 0.25))+
  geom_text(aes(x = -0.36, y = 0.23, label = abclabel, group =exp_ntrt),size = 3,check_overlap = T)




########## RECREATE FIGURE 4 cd, PCoA by decade #####

PCOAandplotsavg_fullN<-PCOAandplotsboth_fullN %>% group_by(exp,year,disk,ntrt) %>% dplyr::summarise(meanAxis1=mean(Axis.1, na.rm=TRUE),sdAxis1=sd(Axis.1, na.rm=TRUE),seAxis1=sd(Axis.1, na.rm=TRUE)/sqrt(n()),meanAxis2=mean(Axis.2, na.rm=TRUE),sdAxis2=sd(Axis.2, na.rm=TRUE),seAxis2=sd(Axis.2, na.rm=TRUE)/sqrt(n()), n=n())


###get rid of years where we don't have 9 or 18 plots###
####2003 for E002 only has 1 field, so remove it###

PCOAandplotsavg2_fullN<-subset(PCOAandplotsavg_fullN, n>8)


PCOAandplotsavg2_fullN<-PCOAandplotsavg2_fullN  %>% mutate(yearlabel=ifelse(year==1982|year==1988|year==1992|year==2004,year, NA))

######

##by decade####

PcoAplaying_fullN <- PCOAandplotsavg2_fullN %>% filter(year %in% c(1982,2004))

arrow_start_fullN <- PCOAandplotsavg2_fullN %>% filter(year ==1982)
arrow_stop_fullN<- PCOAandplotsavg2_fullN %>% filter(year == 2004)

wide_PcoAplaying_fullN <- pivot_wider(PcoAplaying_fullN, names_from = year, values_from = year)

##just pull out the first year of data##

PCOAandplots1982_fullN<-subset(PCOAandplotsavg2_fullN, year==1982)

##pull out the first, last, and midpoint year##

PCOAandplots198219922004_fullN<-subset(PCOAandplotsavg2_fullN, year==1982|year==1992|year==2004)

PCOAandplots198219922004_fullN<-PCOAandplots198219922004_fullN %>% mutate(abclabel=exp)%>% mutate(ntrt2=ntrt) %>%  mutate(expntrt= paste(exp, ntrt, sep = '_'))



PCOAandplots198219922004_fullN$ntrt2  <- mapvalues(PCOAandplots198219922004_fullN$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                                   to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

PCOAandplots198219922004_fullN$ntrt2 <- factor(PCOAandplots198219922004_fullN$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))




#first arrow and second##
PcoAplayingfirst_fullN <- PCOAandplotsavg2_fullN %>% filter(year %in% c(1982,1992, 2004))

PcoAplayingfirst2_fullN<-PcoAplayingfirst_fullN %>% mutate(abclabel=exp)
levels(PcoAplayingfirst2_fullN$abclabel) <- c( "(c) ", "(d) ")

#make sure its in the right order##
PcoAplayingfirst2_fullN<-PcoAplayingfirst2_fullN %>% arrange(year)%>% arrange(exp)%>% arrange(ntrt)%>% mutate(ntrt2=ntrt) %>%  mutate(expntrt= paste(exp, ntrt, sep = '_'))




PcoAplayingfirst2_fullN$ntrt2  <- mapvalues(PcoAplayingfirst2_fullN$ntrt2, from=c("9", "1", "2", "3", "4","5","6", "7", "8"),
                                            to=c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))

PcoAplayingfirst2_fullN$ntrt2 <- factor(PcoAplayingfirst2_fullN$ntrt2, levels = c("None", "0 N+μ", "1 N + μ", "2 N + μ", "3.4 N + μ","5.4 N + μ","9.5 N + μ", "17 N + μ", "27.2 N + μ"))




wide_PcoAplaying1 <- pivot_wider(PcoAplayingfirst2_fullN, names_from = year, values_from = year)

wide_PcoAplaying1$`1982`<-as.numeric(as.character(wide_PcoAplaying1$`1982`))
wide_PcoAplaying1$`1992`<-as.numeric(as.character(wide_PcoAplaying1$`1992`))
wide_PcoAplaying1$`2004`<-as.numeric(as.character(wide_PcoAplaying1$`2004`))


PCOAandplotslabels_fullN<-subset(PcoAplayingfirst2_fullN, ntrt==8)
PCOAandplotslabels_fullN<-PCOAandplotslabels_fullN%>%  mutate(expntrt= paste(exp, ntrt, sep = '_'))%>%  mutate(yearlabel= year)


######

## plot specifications



PCOAandplotsavg2_fullN$disk<-mapvalues(PCOAandplotsavg2_fullN$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))



ann_textdisk2 <- data.frame(meanAxis1 = c(0.03, 0.03), meanAxis2 =c(0.27,0.27),exp= c("Intact in 1982 (E001)","Disturbed in 1982 (E002)"), ntrt2=c("None", "None"), expntrt=c("1_9" ,"2_9"), abclabel = factor(c("(a) ","(b) " ),levels = c("(c) ", "(d) ")))


####
mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))

PcOAdecade_fullN<-ggplot(data=PcoAplayingfirst2_fullN, aes(x=meanAxis1  , y=meanAxis2 , color=ntrt2)) +
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(data=PCOAandplots198219922004_fullN, aes(x=meanAxis1  , y=meanAxis2 , group=ntrt2, shape=as.factor(year)), size=2)+
  geom_errorbar(data=PCOAandplots198219922004_fullN, aes(ymin=meanAxis2-seAxis2, ymax=meanAxis2+seAxis2,color=ntrt2), size=0.25) +
  geom_errorbarh(data=PCOAandplots198219922004_fullN, aes(xmin=meanAxis1-seAxis1, xmax=meanAxis1+seAxis1,color=ntrt2),  size=0.25) +
  geom_point(data=PCOAandplots198219922004_fullN, aes(x=meanAxis1  , y=meanAxis2 , group=ntrt2,color=ntrt2,shape=as.factor(year)))+
  geom_segment(data=wide_PcoAplaying1, aes(x=`1982`  , y=`1982`, xend =`1992`, yend = `1992`, color=ntrt2), arrow=arrow(length=unit(0.0, "cm"))) +
  geom_path(size=0.75, arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+
  geom_segment(data=wide_PcoAplaying1, aes(x=`1992`  , y=`1992`, xend =`2004`, yend = `2004`, color=ntrt2), arrow=arrow(length=unit(0.0, "cm"))) +
  geom_path(size=0.75, arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+
  
  #ylab("PCoA 2")+
  ylab(expression(atop("", paste("PCoA 2"))))+
  xlab("PCoA 1")+
  facet_grid(~disk) + 
  #xlim(-0.330,0.500)+
  #ylim(-0.500,0.180)+
  mytheme+
  scale_shape_manual(values = c(21, 22, 24))+
  theme(strip.background = element_blank())+
  geom_text(aes(x = -0.31, y = 0.18, label = abclabel, group =     abclabel),size = 4, color= "black", check_overlap = T)+theme(strip.text = element_blank())+scale_colour_manual(values = c("dark gray","#f9cb35", "#f98e09", "#e45a31", "#bc3754", "#8a226a", "#57106e", "#210c4a", "#000004"))+theme(legend.position="none")+scale_y_continuous(limits=c(-0.5, 0.18), labels = scales::number_format(accuracy = 0.05))+scale_x_continuous(limits=c(-0.33, 0.5), labels = scales::number_format(accuracy = 0.05))+geom_text(data=PCOAandplotslabels_fullN, aes(label=yearlabel),position = position_nudge(y = -0.04, x=0.05), color="#000004", size=3.0)



#################SUPPLEMENTAL FIGURES###############################


########## RECREATE supplemental figure PCOA scores #####

###PCOA species scores using biodiversity package###


head(spscores_df)

spscores_df$species<-spscores_df$names

### pull out max and min species for each axis##


spscores_dfsubset<-subset(spscores_df, Dim1>0.1|Dim1< (-0.22)|Dim2>0.1|Dim2<(-0.29))


##merge with functional group information##

specieslookup<-da_full %>% distinct(species, .keep_all = TRUE)

specieslookup2<-specieslookup %>% dplyr::select(species,functional.group,duration,lifeform, pathway, origin)


mergedfuncgroupspscores <- (merge(specieslookup2,spscores_dfsubset, by = 'species', sort=FALSE, all.y=TRUE))


####
mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))


speciesPCOAplot <- ggplot() + 
  geom_segment(data=mergedfuncgroupspscores , 
               aes(x=0, y=0, xend=Dim1*2, yend=Dim2*2, group=functional.group),                  size=0.5, arrow=arrow()) +
  geom_text_repel(data=mergedfuncgroupspscores , 
                  aes(x=Dim1*2, y=Dim2*2, label=names, color=functional.group), nudge_y = -0.04)+mytheme+
  theme(legend.position="none")+
  ylab("PCoA 2")+
  xlab("PCoA 1")


 
########## RECREATE supplemental figure PCOA w/ field D #####

PCOA_ABCD_average1982<-PCOAABCDandplots1982 %>% group_by(exp, field) %>% dplyr::summarise(meanAxis1=mean(Axis.1, na.rm=TRUE),sdAxis1=sd(Axis.1, na.rm=TRUE),seAxis1=sd(Axis.1, na.rm=TRUE)/sqrt(n()),meanAxis2=mean(Axis.2, na.rm=TRUE),sdAxis2=sd(Axis.2, na.rm=TRUE),seAxis2=sd(Axis.2, na.rm=TRUE)/sqrt(n()), n=n())


PCOAABCDandplots1982$exp<-mapvalues(PCOAABCDandplots1982$exp, from=c("0", "1", "2"), to=c("Remnant", "Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


PCOAABCDandplots1982$exp<- factor(PCOAABCDandplots1982$exp, levels = c("Remnant", "Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

PCOAABCDandplots1982$field<-mapvalues(PCOAABCDandplots1982$field, from=c("A", "B", "C", "D"), to=c("A (Last farmed 1968)", "B (Last farmed 1957)", "C (Last farmed in 1934)", "D (Never farmed)"))


###change avg too###

PCOA_ABCD_average1982$exp<-mapvalues(PCOA_ABCD_average1982$exp, from=c("0", "1", "2"), to=c("Remnant", "Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

PCOA_ABCD_average1982$exp<- factor(PCOA_ABCD_average1982$exp, levels = c("Remnant", "Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))


PCOA_ABCD_average1982$field<-mapvalues(PCOA_ABCD_average1982$field, from=c("A", "B", "C", "D"), to=c("A (Last farmed 1968)", "B (Last farmed 1957)", "C (Last farmed in 1934)", "D (Never farmed)"))



mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())


PCOA_ABCD_averageplot1982<-ggplot() +
  geom_point(data=PCOAABCDandplots1982, aes(x=Axis.1, Axis.2 , color=field, shape=as.factor(exp)), size=1.5)+
  geom_point(data=PCOA_ABCD_average1982, aes(x=meanAxis1, y=meanAxis2 , color=field, shape=as.factor(exp)), size=5) +
  geom_errorbar(data=PCOA_ABCD_average1982, aes(x= meanAxis1, ymin=meanAxis2-seAxis2, ymax=meanAxis2+seAxis2), color="gray45", size=0.25) +
  geom_errorbarh(data=PCOA_ABCD_average1982, aes(y = meanAxis2, xmin=meanAxis1-seAxis1, xmax=meanAxis1+seAxis1), color="gray45", size=0.25) +
  ylab("PCoA 2")+
  xlab("PCoA 1")+
  mytheme
