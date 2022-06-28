# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

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


############   CONDUCT PCOA (using all the veg data across years) ##############

#make nutrient and years factors#

exp12subset_1$ntrt<-as.factor(exp12subset_1$ntrt)
exp12subset_1$year<-as.factor(exp12subset_1$year)


#arrange the dataframe so that that columns and rows are in right order#
exp12subset_1_sorted<-exp12subset_1%>%arrange(year)%>%arrange(plot)%>%arrange(field)%>%arrange(exp)


#just the plot data#
plot<-exp12subset_1_sorted%>% dplyr::select('field':'ntrt2')

#just the veg data#
veg<-exp12subset_1_sorted%>% dplyr::select(`mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)`:`mass.above.VIOLA SP.`)

#log transformation#
vegmatrix<-as.matrix(veg)
logveg<-log(1+vegmatrix)

#Bray curtis dissimiliarty matrix#
distE1andE2<-vegdist(logveg, method="bray")
BCdistall<-as.matrix(distE1andE2)  ##into a matrix#
BCdistalldf<-as.data.frame(BCdistall)  #into a dataframe#

BCdistalldf2<-cbind(plot,BCdistalldf) #merge plot data w/ dissim matrix##

##BCdistalldf2 has the plot info and the BC dissimilarity matrix for the entire experiment (E001 and E002, all three fields combined, with a subsetted selection of the ntrt treatments)

#PcOA for both E011 and E002 together#

PCOA_allfieldE12 <- pcoa(distE1andE2)  ## takes a while ##

#PCOA_allfieldE12 <- cmdscale(distE1andE2)  ## Another way of doing the pcoa that gives same results#


#extracting the axis scores for each exp#

PCOAaxesboth <- PCOA_allfieldE12$vectors[,c(1,2)]


#merge the axis scores with the plot data \##

PCOAandplotsboth<-cbind(plot,PCOAaxesboth )

#PCOAandplotsboth has the plot info and the PCOA scores for the entire experiment (E001 and E002, all three fields combined, with a subsetted selection of the ntrt treatments)


############# CONDUCT Permanova: nutrient, disturbance, field effect #########


#for this analysis, getting rid of the years where not every field was sampeled#

exp12subset_2<-exp12subset_1%>%filter(year!=1995)%>%filter(year!=1998)%>%filter(year!=2001)%>%filter(year!=2003)

#log transform the biomass data##

justveg<-as.matrix(exp12subset_2[,12:222])

logveg<-log(1+justveg)

exp12subset_2<-cbind(exp12subset_2[,1:11], logveg)

#drop the levels not sampeled

levels(exp12subset_2$year)
exp12subset_2$year <- as.factor(as.character(exp12subset_2$year))

# Extract one year to get the apply function to work#

da.wide_1990<-subset(exp12subset_2, year=="1990", row.names=NULL)

#Do a PERMANOVA looping over years, looking at effect of exp (disturbed or not), ntrt treatment# and field (A, B, C)###

permallfields<-lapply(X = split.data.frame(x = exp12subset_2, f=exp12subset_2$year),function(t){
  z<-adonis2(vegdist(t[,11:175])~exp+ntrt+field, data= da.wide_1990)
  return(data.frame(z))
}
)


#dataframe of results#
permmaster<-plyr::ldply(permallfields, data.frame)
names<-c("Disturbance", "Fertilization", "Field", "residual", "total")

sourcevariation<-rep(names,length(unique(permmaster$.id)))

permmaster2<-cbind(sourcevariation,permmaster)

permmaster3<-subset(permmaster2,sourcevariation=="Disturbance"|sourcevariation=="Fertilization"|sourcevariation=="Field") 

#This dataframe (permmmaster3) has the results from a permanova for each independent variable (disturbance, field, fertilization) for every year in the dataset w/ everything sampeled. The R2 column is what is plotted in Fig 1 ###


######### RECREATE FIGURE 1, variation over time  ##########

#plot specification#
permmaster4<- permmaster3 %>% mutate (year=as.numeric(.id))  %>% mutate (year2=year-1981)

permmaster4$sourcevariation<- factor(permmaster4$sourcevariation, levels = c("Disturbance", "Fertilization", "Field"))

my_plotlabels = paste0('(', letters, ') ')

abclabels = paste0(my_plotlabels[1:length(levels(permmaster4$sourcevariation))]) 

permmaster5<-cbind(permmaster4, abclabels)

##make the plot##

variancetimeall<-permmaster5 %>% ggplot(aes(x=as.numeric(year2), y=R2, group=sourcevariation)) +
  geom_point() +
  facet_wrap(~sourcevariation, strip.position = c("top"))+
  geom_smooth(method="loess", colour="black", size=0.5, fill = "gray1") +
  geom_text(aes(x = 2.0, y = 0.49, label = abclabels, group =     abclabels),size = 4,check_overlap = T)+
  ylab(expression('Explained variation'))+
  xlab("Time since experiment (years)")+
  mytheme+
  theme(strip.text = element_text(face="bold", size=12))+
  theme(strip.background = element_blank())

######### RECREATE FIGURE 2, PCOA ordination by year ########

#first check that there are 18 in each ntrt treatment / year##
PCOAandplotscount<-PCOAandplotsboth %>% count(year,disk,ntrt)

#2003 for E002 only has 1 field so we remove it##

PCOAandplotsboth2<-PCOAandplotsboth %>%filter(year!=2003)

PCOAandplotscount2<-PCOAandplotsboth2 %>% count(year,disk,ntrt)##good now###

PCOAandplotsboth2$Axis.1<-as.numeric(PCOAandplotsboth2$Axis.1)
PCOAandplotsboth2$Axis.2<-as.numeric(PCOAandplotsboth2$Axis.2)

PCOAandplotsboth2$year<-as.numeric(as.character(PCOAandplotsboth2$year))

PCOAandplotsavg<-PCOAandplotsboth2 %>% group_by(exp,year,disk,ntrt) %>% dplyr::summarise(meanAxis1=mean(Axis.1, na.rm=TRUE),sdAxis1=sd(Axis.1, na.rm=TRUE),seAxis1=sd(Axis.1, na.rm=TRUE)/sqrt(n()),meanAxis2=mean(Axis.2, na.rm=TRUE),sdAxis2=sd(Axis.2, na.rm=TRUE),seAxis2=sd(Axis.2, na.rm=TRUE)/sqrt(n()), n=n())

# plot specifications

PCOAandplotsavg2<-mutate(PCOAandplotsavg, ntrt2=ntrt)

PCOAandplotsavg2<-PCOAandplotsavg2  %>% mutate(yearlabel=ifelse(year==1982|year==1992|year==2004,year, NA))%>%mutate(exp_ntrt= as.factor(paste(exp, ntrt, sep = '_')))

PCOAandplotsavg2$exp_ntrt <- factor(PCOAandplotsavg2$exp_ntrt , levels = c("1_9", "2_9", "1_1","2_1", "1_2", "2_2", "1_4", "2_4", "1_6", "2_6"))

PCOAandplotsavg3<-PCOAandplotsavg2 %>% mutate(abclabel=exp_ntrt)

levels(PCOAandplotsavg3$abclabel) <- c( "(a) ", "(b) ", "(c) ", "(d) ", "(e) ", "(f) ", "(g) ", "(h) ", "(i) ", "(j) ")

PCOAandplotsavg3$ntrt2  <- plyr::mapvalues(PCOAandplotsavg3$ntrt2, from=c("9", "1", "2", "3", "4","5","6"),
                                           to=c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))
PCOAandplotsavg3$ntrt2 <- factor(PCOAandplotsavg3$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

PCOAandplotsavg3$disk<-plyr::mapvalues(PCOAandplotsavg3$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

#to draw points between years#

ann_textdisk3 <- data.frame(meanAxis1 = c(0.00, 0.00), meanAxis2 =c(0.27, 0.27),exp= c("Intact in 1982 (E001)","Disturbed in 1982 (E002)"), year=c(1982, 1982), exp_ntrt=c("1_9", "2_9"),abclabel = factor(c("(a) ","(b) " ),levels = c("(a) ", "(b) ", "(c) ", "(d) ", "(e) ", "(f) ", "(g) ", "(h) ", "(i) ", "(j) ")))

ann_textntrt <- data.frame(meanAxis1 = c(0.45, 0.45, 0.45, 0.45, 0.45), meanAxis2 =c(0.05, 0.05, 0.05, 0.05, 0.05),exp= c("2","2","2","2","2"), year=c(1982, 1982, 1982, 1982, 1982), exp_ntrt=c("2_9", "2_1", "2_2", "2_4", "2_6"),abclabel = factor(c("(b) ","(d) ", "(f) ","(h) ","(j) "),levels = c("(a) ", "(b) ", "(c) ", "(d) ", "(e) ", "(f) ", "(g) ", "(h) ", "(i) ", "(j) ")))

#plot#

PcoAfieldAveragedE001E002<-ggplot(PCOAandplotsavg3[order(PCOAandplotsavg3$year),], aes(x=meanAxis1  , y=meanAxis2 , group=exp_ntrt, color=year, label=year)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=meanAxis2-seAxis2, ymax=meanAxis2+seAxis2), color="gray", size=0.25) +
  geom_errorbarh(aes(xmin=meanAxis1-seAxis1, xmax=meanAxis1+seAxis1), color="gray", size=0.25) +
  geom_path(aes(x=meanAxis1 , y=meanAxis2 , group=exp_ntrt, color=year),size=0.5)+
  geom_point(size=1) +
  ylab("PCoA 2")+
  xlab("PCoA 1")+
  facet_grid(ntrt2~disk) + 
  geom_text(aes(label=yearlabel),position = position_nudge(y = -0.04, x=-0.04), color="black", size=2.5)+
  mytheme+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D",trans = 'reverse')+
  guides(#reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE),
    #reverse size order (higher diameter on top) 
    size = guide_legend(reverse = TRUE))+
  scale_y_continuous(breaks=c(-0.2, 0.0, 0.2, 0.4))+
  geom_text(aes(x = -0.37, y = 0.25, label = abclabel, group =     abclabel),size = 3,check_overlap = T)+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(face="bold", size=10))+
  theme(panel.spacing.y=unit(0.5, "lines"))

########## RECREATE FIGURE 4 cd, PCoA by decade #####

PCOAandplotsavg<-PCOAandplotsboth %>% group_by(exp,year,disk,ntrt) %>% dplyr::summarise(meanAxis1=mean(Axis.1, na.rm=TRUE),sdAxis1=sd(Axis.1, na.rm=TRUE),seAxis1=sd(Axis.1, na.rm=TRUE)/sqrt(n()),meanAxis2=mean(Axis.2, na.rm=TRUE),sdAxis2=sd(Axis.2, na.rm=TRUE),seAxis2=sd(Axis.2, na.rm=TRUE)/sqrt(n()), n=n())

#get rid of years where we don't have all 18 plots#
#2003 for E002 only has 1 field, so remove it#

PCOAandplotsavg2<-subset(PCOAandplotsavg, n==18)

PCOAandplotsavg2<-PCOAandplotsavg2%>%mutate(ntrt2 = as.factor(ntrt))

PCOAandplotsavg2<-PCOAandplotsavg2  %>% mutate(yearlabel=ifelse(year==1982|year==1988|year==1992|year==2004,year, NA))

#by decade#

PcoAplaying <- PCOAandplotsavg2 %>% filter(year %in% c(1982,2004))

arrow_start <- PCOAandplotsavg2 %>% filter(year ==1982)
arrow_stop<- PCOAandplotsavg2 %>% filter(year == 2004)

wide_PcoAplaying <- pivot_wider(PcoAplaying, names_from = year, values_from = year)

#just pull out the first year of data#

PCOAandplots1982<-subset(PCOAandplotsavg2, year==1982)

#pull out the first, last, and midpoint year#

PCOAandplots198219922004<-subset(PCOAandplotsavg2, year==1982|year==1992|year==2004)


#first arrow and second##
PcoAplayingfirst <- PCOAandplotsavg2 %>% filter(year %in% c(1982,1992, 2004))

#make sure its in the right order#
PcoAplayingfirst2<-PcoAplayingfirst2 %>% arrange(year)%>% arrange(exp)%>% arrange(ntrt)%>% mutate(ntrt2=ntrt) %>%  mutate(expntrt= paste(exp, ntrt, sep = '_'))

#wide version#
wide_PcoAplaying1 <- pivot_wider(PcoAplayingfirst2, names_from = year, values_from = year)

wide_PcoAplaying1$`1982`<-as.numeric(as.character(wide_PcoAplaying1$`1982`))
wide_PcoAplaying1$`1992`<-as.numeric(as.character(wide_PcoAplaying1$`1992`))
wide_PcoAplaying1$`2004`<-as.numeric(as.character(wide_PcoAplaying1$`2004`))

#plot specification#

PCOAandplots198219922004<-PCOAandplots198219922004 %>% mutate(abclabel=exp)%>% mutate(ntrt2=ntrt) %>%  mutate(expntrt= paste(exp, ntrt, sep = '_'))

levels(PCOAandplots198219922004$abclabel) <- c( "(a) ", "(b) ")

PCOAandplots198219922004$ntrt2  <- mapvalues(PCOAandplots198219922004$ntrt2, from=c("9", "1", "2", "3", "4","5","6"),
                                             to=c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))
PCOAandplots198219922004$ntrt2 <- factor(PCOAandplots198219922004$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

PcoAplayingfirst2<-PcoAplayingfirst %>% mutate(abclabel=exp)
levels(PcoAplayingfirst2$abclabel) <- c( "(c) ", "(d) ")

PcoAplayingfirst2$ntrt2  <- mapvalues(PcoAplayingfirst2$ntrt2, from=c("9", "1", "2", "3", "4","5","6"),
                                      to=c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))
PcoAplayingfirst$ntrt2 <- factor(PcoAplayingfirst$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

PCOAandplotslabels<-subset(PcoAplayingfirst2, ntrt==4)
PCOAandplotslabels<-PCOAandplotslabels%>%  mutate(expntrt= paste(exp, ntrt, sep = '_'))%>%  mutate(yearlabel= year)

PCOAandplotsavg2$ntrt2  <- mapvalues(PCOAandplotsavg2$ntrt2, from=c("9", "1", "2", "3", "4","5","6"),
                                     to=c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))
PCOAandplotsavg2$ntrt2 <- factor(PCOAandplotsavg2$ntrt2, levels = c("None", "0 N+micro", "1 N + micro", "2 N + micro", "3.4 N + micro","5.4 N + micro","9 N + micro"))

PCOAandplotsavg2$disk<-mapvalues(PCOAandplotsavg2$disk, from=c("0", "1"), to=c("Intact in 1982 (E001)", "Disturbed in 1982 (E002)"))

ann_textdisk2 <- data.frame(meanAxis1 = c(0.03, 0.03), meanAxis2 =c(0.27,0.27),exp= c("Intact in 1982 (E001)","Disturbed in 1982 (E002)"), ntrt2=c("None", "None"), expntrt=c("1_9" ,"2_9"), abclabel = factor(c("(a) ","(b) " ),levels = c("(c) ", "(d) ")))

# make plot# 
PcOAdecade<-ggplot(data=PcoAplayingfirst2, aes(x=meanAxis1  , y=meanAxis2 , color=ntrt2)) +
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(data=PCOAandplots198219922004, aes(x=meanAxis1  , y=meanAxis2 , group=ntrt2, shape=as.factor(year)), size=2)+
  geom_errorbar(data=PCOAandplots198219922004, aes(ymin=meanAxis2-seAxis2, ymax=meanAxis2+seAxis2,color=ntrt2), size=0.25) +
  geom_errorbarh(data=PCOAandplots198219922004, aes(xmin=meanAxis1-seAxis1, xmax=meanAxis1+seAxis1,color=ntrt2),  size=0.25) +
  geom_point(data=PCOAandplots198219922004, aes(x=meanAxis1  , y=meanAxis2 , group=ntrt2,color=ntrt2,shape=as.factor(year)))+
  geom_segment(data=wide_PcoAplaying1, aes(x=`1982`  , y=`1982`, xend =`1992`, yend = `1992`, color=ntrt2), arrow=arrow(length=unit(0.0, "cm"))) +
  geom_path(size=0.75, arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+
  geom_segment(data=wide_PcoAplaying1, aes(x=`1992`  , y=`1992`, xend =`2004`, yend = `2004`, color=ntrt2), arrow=arrow(length=unit(0.0, "cm"))) +
  geom_path(size=0.75, arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+
  geom_text(data=PCOAandplotslabels, aes(label=yearlabel),position = position_nudge(y = -0.018, x=-0.02), color="black", size=3.0)+
  ylab("PCoA 2")+
  xlab("PCoA 1")+
  facet_grid(~disk) + 
  xlim(-0.28,0.375)+
  ylim(-0.35,0.15)+
  scale_colour_manual(values = c("gray45", "#E69F00", "#009E73","#0072B2","#CC79A7"))+
  mytheme+
  theme(legend.position="none")+
  scale_shape_manual(values = c(21, 22, 24))+
  theme(strip.background = element_blank())+
  geom_text(aes(x = -0.25, y = 0.14, label = abclabel, group =     abclabel),size = 4, color= "black", check_overlap = T)+theme(strip.text = element_blank())
