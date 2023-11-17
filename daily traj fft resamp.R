devtools::install_github("JimMcL/trajr")

library(trajr)
library(sp)
library(raster)


data=read.csv("C:/Users/salavi/Documents/FFT_Cleaned_Final.csv")

data$timestamp=as.POSIXct((data$timestamp),format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900",tz="UTC")



data2 <- SpatialPointsDataFrame(coords = data[,c(3,4)], data = data,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
data2 <- spTransform(data2, CRS("+proj=utm +zone=17 +north +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
data2=as.data.frame(data2)


data2$loopID=paste(data2$individual.local.identifier,data2$day,sep="_")
data2$loopID=as.factor(data2$loopID)
data2$Z=data2$location.long.2+1i*data2$location.long.2

datasplit=split(data2,data2$loopID)
datasplit2=split(data2,as.factor(data2$individual.local.identifier))
resolutions=c(36.39182,	28.98044,	32.76763,	37.87446,	42.17405,	12.70379,	23.85797,	40.80566,	31.37489,	29.64003,	20.45771,	30.39861,	36.82671,	32.24596,	20.1786,	34.8374,	46.27092,	40.84662,	20.579,	49.20868,	28.34774,	23.71034,	25.30403,	38.16625,	27.09191,	24.83083,	25.55506,	32.09851,	27.50105,	17.6466,	23.66956,	36.17021,	27.50515,	22.88032,	18.97634,	19.11565,	13.49218,	20.84271,	17.82411,	30.6463,	42.82088,	39.76968,	34.35561,	37.84272,	33.43426,	25.17701,	21.16828,	39.55887,	38.62099)

daily=c()
for(i in 1:length(datasplit)){
  track=datasplit[[i]]
  coords=datasplit[[i]][,c(35,36,2)]
  colnames(coords)=c("x","y","time")
  day=unique(track$day)
  date=as.character(as.Date(min(track$timestamp)))
  ID=unique(track$individual.local.identifier)
  Species=unique(track$individual.taxon.canonical.name)
  
  track2=TrajFromCoords(coords)
  track2=TrajRediscretize(track2, 1)
  straightness=TrajStraightness(track2)
  sinuosity=TrajSinuosity2(track2)
  x=c(ID,Species,day,date,straightness,sinuosity)
  daily[i]=list(x)
  
}
daily=do.call(rbind,daily)
daily=data.frame(daily)

daily2=c()
for(i in 1:length(datasplit)){
  track=datasplit[[i]]
  coords=datasplit[[i]][,c(35,36,2)]
  colnames(coords)=c("x","y","time")
  day=unique(track$day)
  date=as.character(as.Date(min(track$timestamp)))
  ID=unique(track$individual.local.identifier)
  Species=unique(track$individual.taxon.canonical.name)
  
  track2=TrajFromCoords(coords)
  track2=try(TrajRediscretize(track2, 21.45333))
  if(class(track2)=="try-error"){next}
  straightness=TrajStraightness(track2)
  sinuosity=TrajSinuosity2(track2)
  x=c(ID,Species,day,date,straightness,sinuosity)
  daily2[i]=list(x)
  
}
daily2=do.call(rbind,daily2)
daily2=data.frame(daily2)
colnames(daily)=c("ID","Species","Study_Day","Date","Straightness","Sinuosity")
colnames(daily2)=c("ID","Species","Study_Day","Date","Straightness","Sinuosity")

write.csv(daily,file="FFT Straightness and Sinuosity.csv")



daily3=c()
count=1
for(i in 1:length(datasplit2)){
  temp=datasplit2[[i]]
  datasplit3=split(temp,temp$loopID)
  
  for(j in 1:length(datasplit3)){
    track=datasplit3[[j]]
    coords=datasplit3[[j]][,c(35,36,2)]
    colnames(coords)=c("x","y","time")
    day=unique(track$day)
    date=as.character(as.Date(min(track$timestamp)))
    ID=unique(track$individual.local.identifier)
    Species=unique(track$individual.taxon.canonical.name)
    
    track2=TrajFromCoords(coords)
    track2=try(TrajRediscretize(track2, resolutions[i]))
    if(class(track2)=="try-error"){next}
    straightness=TrajStraightness(track2)
    sinuosity=TrajSinuosity2(track2)
    x=c(ID,Species,day,date,straightness,sinuosity)
    daily3[count]=list(x)
    count=count+1
  }
  
  
}
daily3=do.call(rbind,daily3)
daily3=data.frame(daily3)
colnames(daily3)=c("ID","Species","Study_Day","Date","Straightness","Sinuosity")


newdata <- merge(daily, daily2, by.x = c("ID","Species","Study_Day","Date"), by.y = c("ID","Species","Study_Day","Date"), all.x = TRUE)
newdata <- merge(newdata, daily3, by.x = c("ID","Species","Study_Day","Date"), by.y = c("ID","Species","Study_Day","Date"), all.x = TRUE)
colnames(newdata)=c("ID","Species","Study_Day","Date","Straightness.1", "Sinuosity.1", "Straightness.MSL", "Sinuosity.MSL",    "Straightness.TAU",   "Sinuosity.TAU")


plot(track2, lwd = 2)
points(track2, draw.start.pt = FALSE, pch = 16, col = "black", cex = 1.2)

# Resample to step length 1
resampled1 <- TrajRediscretize(track2, 1)
resampled2 <- TrajRediscretize(track2, 21)
resampled3 <- TrajRediscretize(track2, resolutions[i])

# Plot rediscretized trajectory in red
lines(resampled3, col = "#FF0000A0", lwd = 1)
points(resampled3, type = 'p', col = "#FF0000A0", pch = 16)


daily4=c()
for(i in 1:length(datasplit)){
  track=datasplit[[i]]
  coords=datasplit[[i]][,c(35,36,2)]
  colnames(coords)=c("x","y","time")
  day=unique(track$day)
  date=as.character(as.Date(min(track$timestamp)))
  ID=unique(track$individual.local.identifier)
  Species=unique(track$individual.taxon.canonical.name)
  
  track2=TrajFromCoords(coords)
  track2=try(TrajRediscretize(track2, 29.84784))
  if(class(track2)=="try-error"){next}
  straightness=TrajStraightness(track2)
  sinuosity=TrajSinuosity2(track2)
  x=c(ID,Species,day,date,straightness,sinuosity)
  daily4[i]=list(x)
  
}
daily4=do.call(rbind,daily4)
daily4=data.frame(daily4)
colnames(daily4)=c("ID","Species","Study_Day","Date","Straightness","Sinuosity")


newdata <- merge(newdata, daily4, by.x = c("ID","Species","Study_Day","Date"), by.y = c("ID","Species","Study_Day","Date"), all.x = TRUE)
newdata2=na.omit(newdata)

colnames(newdata)=c("ID","Species","Study_Day","Date","Straightness.1", "Sinuosity.1", "Straightness.MSL", "Sinuosity.MSL",    "Straightness.TAU",   "Sinuosity.TAU","Straightness.MT",   "Sinuosity.MT")
newdata2=na.omit(newdata)

write.csv(newdata2,file="FFT_sinuosity_straightness_different_resamping.csv")

for(i in 5:ncol(newdata2)){
  newdata2[,i]=as.numeric(newdata2[,i])
}
newdata2=newdata2[-which(newdata2$Sinuosity.MSL>1),]
newdata2=newdata2[-which(newdata2$Sinuosity.TAU>3),]
newdata2=newdata2[which(newdata2$Sinuosity.MT>5),]

library(ggplot2)
ggplot(newdata2)+geom_boxplot(aes(x=Species,y=Sinuosity.1,fill=Species))+theme_classic()
ggplot(newdata2)+geom_boxplot(aes(x=Species,y=Sinuosity.MSL,fill=Species))+theme_classic()
ggplot(newdata2)+geom_boxplot(aes(x=Species,y=Sinuosity.TAU,fill=Species))+theme_classic()
ggplot(newdata2)+geom_boxplot(aes(x=Species,y=Sinuosity.MT,fill=Species))+theme_classic()
