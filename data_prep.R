#####################################################
##      Metapopulation model                       ##
##      northern gannet                            ##
##      Jeglinski et al                            ##
#####################################################


### 1. data preparation


# colony census data
gan1<-read.csv('gannets_master_complete_updated.csv')

# colonies summary
colonies1<-read.csv('Colonies.csv') 

gan1$regional_seas<-as.character(gan1$regional_seas)

gan1<-gan1[-19:-28] # Drops some NA columns
gan1<-gan1[gan1$regional_seas!="Rockall Trough and Bank",]

unique(gan1$regional_seas)


# generate aditional regiosn for remote colonies
unique(gan1$Colony)
gan1$regional_seas<-ifelse(gan1$Colony=="Bjornoya","Bjornoya",
                           ifelse(gan1$Colony=="Runde","Southern Norwegian Sea",
                                  ifelse(gan1$Colony=="Kharlov Kola Peninsula","Kharlov Kola Peninsula",
                                         ifelse(gan1$Colony=="Rouzic","South West Atlantic", 
                                                ifelse(gan1$Colony=="Ortac","South West Atlantic",
                                                       ifelse(gan1$Colony=="Les Etacs","South West Atlantic",
                                                              gan1$regional_seas))))))


unique(gan1$regional_seas) # 15 regional seas

gan1$regional_seas<-as.factor(gan1$regional_seas)
#gan1$regional_seas<-droplevels(gan1$regional_seas)
str(gan1)
names(gan1)

gan1$AOS <-as.integer(gan1$AOS)
gan1$colonization_year<-as.character(gan1$colonization_year)
gan1$colonization_year<-as.numeric(gan1$colonization_year)
gan1$extinction_year<-as.character(gan1$extinction_year)
gan1$extinction_year<-as.numeric(gan1$extinction_year)
gan1<-cbind(gan1,region=as.numeric(as.factor(gan1$regional_seas))) # regional seas become numeric regions for matrix below

cols<-unique(gan1$Colony)

## exclude colonisation attempts


str(colonies1)
nrow(colonies1)                                                              ## 66 colonies & colonisation events
unique(colonies1[which(colonies1$Status=="colonisation attempt"),4])
length(unique(colonies1[which(colonies1$Status=="colonisation attempt"),4])) ## 13 colonisation events to remove
cols1<-subset(colonies1,colonies1$Status!="colonisation attempt")
cols<-unique(cols1$Colony)

##################### prepare data for model

#required
N<-length(cols)
nYrs<-117


par(mfrow=c(3,4))

Pd<-matrix(NA, nrow=N, ncol=nYrs) # Matrix of population sizes (colonies in rows, years in columns)

harvest<-Pd*0 # Matrix with information on hunting (yes/no)


reg<-rep(0,N) # Vector of region to which each colony belongs

cc<-rep(0,N) # Vector of terrestrial ccs based on experts' estimates



# Loop filling data into matrices

for(n in 1:N)
{
  
  Pd[n,]<-gan1[gan1$Colony==cols[n],]$AOS
  reg[n]<-gan1[gan1$Colony==cols[n],]$region[1]               
  harvest[n,]<-gan1[gan1$Colony==cols[n],]$harvest
  harvest[n,is.na(harvest[n,])]=0
  
  cc[n]<-colonies1[as.character(colonies1$Colony)==as.character(cols[n]),]$Terrestrial.CC
  plot(Pd[n,], main=paste(n,". ",cols[n]), col=harvest[n,]+1, ylim=c(0,cc[n]))
}



Pfull<-Pd

nreg<-length(unique(reg)) # Number of regions


# Forecasting
horizon<-20

Pd<-cbind(Pfull,matrix(NA,nrow=N, ncol=horizon))


harvest<-cbind(harvest,matrix(harvest[,nYrs],nrow=N, ncol=horizon))

nYrs<-nYrs+horizon                                                  # 137 years

