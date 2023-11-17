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
colnames(daily)=c("ID","Species","Study_Day","Date","Straightness","Sinuosity")
write.csv(daily,file="FFT Straightness and Sinuosity.csv")
