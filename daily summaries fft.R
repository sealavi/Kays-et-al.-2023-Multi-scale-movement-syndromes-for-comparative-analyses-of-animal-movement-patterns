install.packages("circmax", repos="http://R-Forge.R-project.org",depdependencies = TRUE)
data=read.csv("C:/Users/salavi/Documents/FFT_Cleaned_Final.csv")
data$timestamp=as.POSIXct(as.character(data$timestamp),format="%Y/%m/%d %H:%M:%S",origin="01-01-1900",tz="UTC")
#data$timestamp<-as.POSIXct(as.character(data$timestamp), format="%Y-%m-%d %H:%M:%S",origin="01-01-1900", tz="UTC")
datasplit=split(data,data$loopID)

datasplit2=split(data,data$individual.local.identifier)


library(circmax)
library(fitdistrplus)

newdat=c()
for(i in 1:length(datasplit)){
  temp=datasplit[[i]]  
  Z=temp$utm.easting+1i*temp$utm.northing
  steplengths=as.numeric(Mod(diff(Z)))
  turnangles=as.numeric(Arg(diff(Z)))
  fit <- try(fitdist(steplengths, "gamma", method="mle"))
  if(class(fit)=="try-error"){
    fit <- fitdist(steplengths, "gamma", method="mme")
  }
  shape=fit$estimate[[1]]
  rate=fit$estimate[[2]]
  
  vm=circfit(turnangles)
  Mu=vm$coefficients[[1]]
  Kappa=vm$coefficients[[2]]
  
  Date=format(mean(temp$timestamp,na.rm=TRUE),format="%Y/%m/%d")
  ID=unique(as.character(temp$individual.local.identifier))
  Species=unique(as.character(temp$individual.taxon.canonical.name))
  outputs=c(ID,Species,Date,shape,rate,Mu,Kappa)
  newdat[i]=list(outputs)
  print(plot(rvonmises(1000, mu=circular(Mu), kappa=Kappa, control.circular=list(units="radians"))))
  
}
newdat=do.call(rbind,newdat)
newdat=data.frame(newdat)
colnames(newdat)=c("ID","Species","Date","shape","rate","Mu","Kappa")
newdat$Date=as.POSIXct(newdat$Date,format="%Y/%m/%d",origin="01-01-1900",tz="UTC")
newdat= newdat[with(newdat, order(ID, Date)), ]

library(circular)
plot(rvonmises(200, mu=circular(Mu), kappa=Kappa, control.circular=list(units="radians")))

write.csv(newdat,file="daily_distributions.csv")



newdat2=c()

for(i in 1:length(datasplit)){
  temp=datasplit[[i]]  
  temp= temp[with(temp, order(individual.local.identifier, timestamp)), ]
  Z=temp$utm.easting+1i*temp$utm.northing
  steplengths=as.numeric(Mod(diff(Z)))
  straight_line_travel_distance=sum(steplengths)
  Date=format(mean(temp$timestamp,na.rm=TRUE),format="%Y/%m/%d")
  StudyDay=unique(temp$day)
  ID=unique(as.character(temp$individual.local.identifier))
  Species=unique(as.character(temp$individual.taxon.canonical.name))
  outputs=c(ID,Species,Date,StudyDay,straight_line_travel_distance)
  newdat2[i]=list(outputs)

}

newdat2=do.call(rbind,newdat2)
newdat2=data.frame(newdat2)
colnames(newdat2)=c("ID","Species","Date","StudyDay","straight_line_travel_distance")
newdat2$Date=as.POSIXct(newdat2$Date,format="%Y/%m/%d",origin="01-01-1900",tz="UTC")
newdat2= newdat2[with(newdat2, order(ID, StudyDay)), ]
write.csv(newdat2,file="daily_straight_line_travel_distance.csv")



newdat3=c()
for(i in 1:length(datasplit2)){
  temp=datasplit2[[i]]  
  Z=temp$utm.easting+1i*temp$utm.northing
  steplengths=as.numeric(Mod(diff(Z)))
  turnangles=as.numeric(Arg(diff(Z)))
  fit <- try(fitdist(steplengths, "gamma", method="mle"))
  if(class(fit)=="try-error"){
    fit <- fitdist(steplengths, "gamma", method="mme")
  }
  shape=fit$estimate[[1]]
  rate=fit$estimate[[2]]
  
  vm=circfit(turnangles)
  Mu=vm$coefficients[[1]]
  Kappa=vm$coefficients[[2]]
  ID=unique(as.character(temp$individual.local.identifier))
  Species=unique(as.character(temp$individual.taxon.canonical.name))
  outputs=c(ID,Species,shape,rate,Mu,Kappa)
  newdat3[i]=list(outputs)

}
newdat3=do.call(rbind,newdat3)
newdat3=data.frame(newdat3)
colnames(newdat3)=c("ID","Species","shape","rate","Mu","Kappa")
newdat3= newdat3[with(newdat3, order(ID, Date)), ]
write.csv(newdat3,file="Study_wide_step_turn_distributions.csv")
