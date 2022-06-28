# MH DeSiervo, LG Shoemaker
# Nutrient supply shifts successional paths but not speed of grassland recovery from disturbance

# Code for functions used in subsequent analyses####

library(ggplot2)

########AIC TABLE FUNCTION#########

aictable<-function(X,m){
  #
  #  This function will create a full model selection table based on AICc.
  #  
  #  Inputs:
  #
  #  X is the AIC table produced by R by using AIC(model1,model2, ...)
  #  m is the sample size
  #
  #
  rnames<-row.names(X)
  AICc<-X$AIC+2*X$df*(X$df+1)/(m-X$df-1)     #small-sample correction
  logL<-X$df-X$AIC/2                         #Log-likelihood
  tab<-data.frame(X[,1],logL,AICc)           #Remove AIC column; add logL and AICc
  colnames(tab)[1]<-c("Params")              #Rename "df" column   
  row.names(tab)<-rnames
  tab<-tab[order(tab$AICc),]                 #Sort by ascending AICc value
  deltaAICc<-tab$AICc-min(tab$AICc)          #Delta AICc
  weight<-exp(-deltaAICc/2)/sum(exp(-deltaAICc/2))  #Weights
  cumwt<-weight                              #Column for cumulative weight
  for(i in 2:dim(X)[1]){                  
    cumwt[i]<-cumwt[i-1]+cumwt[i]              #Accumulate weight from the top
  }
  tab<-data.frame(tab,deltaAICc,weight,cumwt)
  tab<-round(tab,4)
  tab
}

#########THEME FOR GG PLOT###

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))+ theme(panel.grid.major = element_line(colour = "gray 85"))+ theme(panel.grid.minor = element_line(colour = "gray 95"))


