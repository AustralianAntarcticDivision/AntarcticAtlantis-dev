library(raadtools)
files  <- icefiles()
range(files$date)
## use calendar months 1998-2018
library(dplyr)
files <- files %>%
  dplyr::filter(between(as.Date(date), as.Date("1998-01-01"), as.Date("2018-12-31"))) %>%
  mutate(month = format(date, "%B"))

range(files$date)
umonths <- unique(files$month)
## loop over months
monthlist <- setNames(vector("list", length(umonths)),
                      umonths)
for (i in seq_along(umonths)) {
  f0 <- files %>% dplyr::filter(month == umonths[i])
  monthlist[[i]] <- mean(readice(f0$date))
}
saveRDS(monthlist, "monthly_seaice_mean_1998-2018.rds")
library(rbgm)
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- bgmfile(bfile)
boxes <- boxSpatial(bgm)
## turn the monthly ice grids into a brick

ice <- brick(monthlist)
names(ice) <- umonths

ice_means <- raster::extract(ice, boxes, fun = mean, na.rm = TRUE)

saveRDS(ice_means, "table_ice_means_28.rds")



######################################
############### Mao eddit #######################
######################################

setwd("/Users/mmori/Atlantis_East_Ant/Mao_seaice")

monthlist<-readRDS("monthly_seaice_mean_1998-2018.rds", refhook = NULL)

load("/Users/mmori/Atlantis_East_Ant/creating inputfile/BGM_box.Rdata")
ice <- brick(monthlist)
season<- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
names(ice) <- season
ice_all <- raster::extract(ice,box, weights = TRUE)


tablei<-list()
for (m in 1:12){
month<-matrix(0,2,28)
for (p in 1:28){
sample<-ice_all[[p]][,m]
badna<-which(is.na(sample))
bad<-length(badna)
good<-which(sample[!is.na(sample)] >=50)
month[1,p]<-round(length(good)/(length(sample)-bad)*100,1)
month[2,p]<-round(length(badna)/(length(sample))*100,1)
}
tablei[[m]]<-month
}

tableice<-readRDS("table_ice_means_28.rds", refhook = NULL)

######## Plot box #############
library(fields)

pdf("/Users/mmori/Atlantis_East_Ant/Mao_seaice/mean_ice.pdf",width=14, height=14)
plot(crop(ice, spTransform(box, projection(ice))), addfun = function() plot(spTransform(box, projection(ice)), add = TRUE))
dev.off()

colp<-c(c(colorRampPalette(c("white","blue","cyan","yellow","orange","red","red4"))(90)),c(colorRampPalette(c("black"))(11)))

pdf("/Users/mmori/Atlantis_East_Ant/Mao_seaice/ice_50_1.pdf",width=14, height=14)
par(mfrow=c(3,4))

for (m in 1:12){
par(mar=c(5.1,2.1,4.1,5.1))
plot(box,col=colp[round(tablei[[m]][1,])+1],main=season[m])#total number 36429
image.plot(box,legend.only=TRUE, zlim=c(0,100), col=colp,legend.cex=1,legend.line = 0.1,legend.width=1.0,legend.mar=c(6.1,2.1,4.1,5.1),legend.shrink=0.6)
}
par(mfrow=c(1,1))

dev.off()


pdf("/Users/mmori/Atlantis_East_Ant/Mao_seaice/ice_15nan.pdf",width=14, height=14)
par(mfrow=c(3,4))

for (m in 1:12){
  par(mar=c(5.1,2.1,4.1,5.1))
  plot(box,col=colp[round(tablei[[m]][2,])+1],main=season[m])#total number 36429
  image.plot(box,legend.only=TRUE, zlim=c(0,100), col=colp,legend.cex=1,legend.line = 0.1,legend.width=1.0,legend.mar=c(6.1,2.1,4.1,5.1),legend.shrink=0.6)
}
par(mfrow=c(1,1))

dev.off()
