#### monthly sea ice
setwd("/Users/mmori/Atlantis_East_Ant/Sea_ice/Mao_seaice")

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


####

library(ncdf4)
setwd("/Users/mmori/Atlantis_East_Ant/Sea_ice")
# Define netcdf dimensions: name, units and value
dim1 <- ncdim_def( # create a time dimension
  name = 't',
  units = 'seconds since 2000-01-01 00:00:00 +10',
  unlim = TRUE,
  vals = as.double(seq(0,86400*365,86400))
)
#as.double(seq(0,86400*365,86400)
dim2 <- ncdim_def( # create a box dimension
  name = 'b',
  units = 'num',
  vals = as.double(1:28)
)

var.dim <- list(dim2, dim1)

ver1<-ncvar_def(name = as.character("Ice_Class1"),
                units = '(none)', dim = var.dim, missval = 0, prec="double")

ver2<-ncvar_def(name = as.character("t(t)"),units = 'seconds since 2000-01-01 00:00:00 +10', dim = dim1)
ver3<-ncvar_def(name = as.character("total_depth"),
                units = 'm', dim = var.dim, missval =0, prec="double")

filename <- "/Users/mmori/Atlantis_East_Ant/Sea_ice/sea_ice_test.nc"
outnc <- nc_create(filename,  list(ver1, ver3), force_v4 = TRUE)




# add global attributes
ncatt_put(nc = outnc, varid = 0, attname = 'title', attval = "trivial")
ncatt_put(nc = outnc, varid = 0, attname = 'geometry', attval = "Antarctica_28.bgm")
ncatt_put(nc = outnc, varid = 0, attname = 'parameters', attval = "")
ncatt_put(nc = outnc, varid = 0, attname = 'history', attval = "Fri Jan  25 13:33:40 2019: ncks -d t,0,0 sea_ice_test.nc")
ncatt_put(nc = outnc, varid = 0, attname = 'NCO', attval = "4.0.8")


ice_50<-ice_D<-matrix(0, nrow = 28, ncol = length(seq(0,86400*365,86400)))
monthd<-c(31,28,31,30,31,30,31,31,30,31,30,32)
st<-c(1,32,60,91,121,151,181,212,243,273,304,334)
ed<-c(31,59,90,120,150,180,211,242,272,303,333,365)
for (m in 1:12){
  
who<-  which(tablei[[m]][1,]>=90)
for (i in st[m]:ed[m]){
  if (length(who)>0){
  ice_50[who,i]<-1
  ice_D[who,i]<-1
  }
}
}  
  
  


ncvar_put(outnc,varid ="Ice_Class1",ice_50)
ncvar_put(outnc,varid ="total_depth",ice_D)

nc_close(outnc)