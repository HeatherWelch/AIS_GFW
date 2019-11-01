## data exploration 10_29_19
# summarizing points in kernal density polygons from arcgis

library(tidyverse)
library(tiff)
library(raster)
library(gbm)
library(dismo)
library(ncdf4)
library(tidyverse)
library(MASS)
library(raster)
library(ggmap)
library(cowplot)
library(d3heatmap)
library(maps)
library(fields)
library(mapdata)
library(rasterVis)
library(corrplot)
library(factoextra)
library(gplots)
library(sf)
library(transformr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(maptools)
library(rmapshaper)

dev_eval=function(model_object){
  null <- model_object$self.statistics$null.deviance
  res <- model_object$self.statistics$resid.deviance
  dev=((null - res)/null)*100 
  return(dev)
}

extracto=function(pnts,raster,colname){
  pnts=pnts %>% as.data.frame()
  b=pnts%>% as.data.frame() %>% dplyr::select(x,y) 
  a=extract(raster,b)
  pnts$new=a
  position=grep("new",names(pnts))
  colnames(pnts)[position]=colname
  return(pnts)
}

# all kernal density polygons ####
# a=read.csv("~/Desktop/Win7Desktop/GFW_IUU/all_density.csv")
a=read.csv("~/Desktop/Win7Desktop/GFW_IUU/all_density_fd.csv")
b=read.csv("~/Desktop/Win7Desktop/GFW_IUU/filtered_gaps_10_24_29.csv") %>% filter(frac_day_coverage>.2)
# write.csv(b,"~/Desktop/Win7Desktop/GFW_IUU/filtered_gaps_10_24_29_frac_day.csv")
a$location=NA
a <- within(a, location[ID==1] <- 'West_Atlantic')
a <- within(a, location[ID==2] <- 'East_Atlantic')
a <- within(a, location[ID==3] <- 'West_Pacific')
# a <- within(a, location[ID==4] <- 'Indian')
a <- within(a, location[ID==5] <- 'East_Pacific')
a <- within(a, location[ID==6] <- 'Patagonian_Shelf')

c=a %>% .[complete.cases(.),] %>% group_by(flag,location) %>% summarise(n=n())
ggplot(c,aes(x=reorder(flag,-n),y=n))+geom_bar(stat="identity")+facet_wrap(~as.factor(location),scales="free_x")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()

c=a %>% .[complete.cases(.),] %>% group_by(vessel_cla,location) %>% summarise(n=n())
ggplot(c,aes(x=reorder(vessel_cla,-n),y=n))+geom_bar(stat="identity")+facet_wrap(~as.factor(location),scales="free_x")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()

c=a %>% .[complete.cases(.),] %>% mutate(vessel_flag=paste0(vessel_cla,"_",flag)) %>% group_by(vessel_flag,location) %>% summarise(n=n())%>%group_by(location) %>% 
  mutate(rank = rank(desc(n))) %>% arrange(rank) %>% filter(rank<6)%>% as.data.frame() 
ggplot(c,aes(x=vessel_flag,y=n))+geom_bar(aes(fill=vessel_flag),stat="identity")+facet_wrap(~as.factor(location),scales="free_x")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()

c=a %>% .[complete.cases(.),] %>% group_by(location) %>% summarise(n=mean(frac_day_c)) %>% arrange(n)


# building some models ####
# turn raster into points
# extract bathy, chla, sd_bathy, sd_chla
# covariates: bathy, chla, sd_bathy, sd_chla, distance from shore
# possible covariates: flag, gear, distance from port, gap_hrs

# turn raster into points
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
b=a %>% filter(frac_day_coverage>.2)
c=b %>% .[complete.cases(.),] %>%  dplyr::select(on_lat,on_lon)
c=c[,2:1]
x=raster()
res(x)=1
crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
ras=rasterize(c,x,fun=sum)
ras_pnt=rasterToPoints(ras) %>% as.data.frame()

dist_shore=b %>% .[complete.cases(.),]%>% dplyr::select(c(on_distance_from_shore_m))
ras_dist_shore=rasterize(c,x,vals=dist_shore,fun=mean)
master_pnt=rasterToPoints(ras_dist_shore) %>% as.data.frame() %>% rename(dist_shore=layer) %>% mutate(sum_gaps=ras_pnt$layer)

# spatial_points=b %>% .[complete.cases(.),] %>% dplyr::select(c(on_lon,on_lat,on_distance_from_shore_m))
# coordinates(spatial_points)=~on_lon+on_lat
# crs(spatial_points)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
# ras_dist_shore=raster::rasterize(spatial_points,raster=x,field="on_distance_from_shore_m",fun=mean)
# ras_pnt_dist_shore=rasterToPoints(ras_dist_shore)

# creating covariate rasters for extraction ####
# bathy=raster("~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/ETOPO1_Bed_g_geotiff.tif")
# modis=raster("/Users/heatherwelch/Dropbox/JPSS/global/global_averages_mask/modis.grd")
# 
# bathy_res=resample(bathy,ras)
# modis_res=resample(modis,ras)
# bathy_res_sd=focal(bathy_res,w=matrix(1,nrow=3,ncol=7),fun=sd,na.rm=T,pad=T)
# modis_res_sd=focal(modis_res,w=matrix(1,nrow=3,ncol=7),fun=sd,na.rm=T,pad=T)

# writeRaster(bathy_res_sd,"/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1_sd.grd",overwrite=T)
# writeRaster(modis_res_sd,"/Volumes/SeaGate/IUU_GRW/global_shp/modis_1_sd.grd",overwrite=T)
# writeRaster(bathy_res,"/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1.grd",overwrite=T)
# writeRaster(modis_res,"/Volumes/SeaGate/IUU_GRW/global_shp/modis_1.grd",overwrite=T)

# #temperature
# nc=nc_open("/Volumes/SeaGate/IUU_GRW/global_shp/otemp.anal1deg.nc")
# temp_ras=raster("/Volumes/SeaGate/IUU_GRW/global_shp/otemp.anal1deg.nc") %>% rotate()
# values(x)=1
# temp_ras=resample(temp_ras,x)
# writeRaster(temp_ras,"/Volumes/SeaGate/IUU_GRW/global_shp/sst.grd",overwrite=T)

# # distance to shore
# data(wrld_simpl)
# #Rasterize and set land pixels to NA ---> use one of these to crease eez buffer [writen out raster, commented out code]
# r2 <- rasterize(wrld_simpl, x, 1)
# #r3 <- mask(is.na(r2), r2, maskvalue=1, updatevalue=NA) don't use this

# Calculate distance to nearest non-NA pixel
# d <- distance(r2)
# writeRaster(d,"/Volumes/SeaGate/IUU_GRW/global_shp/dist_shore.grd",overwrite=T)

# # lat and long rasters
# lat <- lon <- x
# xy <- coordinates(x)
# lon[] <- xy[, 1]
# lat[] <- xy[, 2]
# writeRaster(lon,"/Volumes/SeaGate/IUU_GRW/global_shp/lon.grd",overwrite=T)
# writeRaster(lat,"/Volumes/SeaGate/IUU_GRW/global_shp/lat.grd",overwrite=T)

## trying a buffer
# buf=buffer(r2,width=370400)
# writeRaster(buf,"/Volumes/SeaGate/IUU_GRW/global_shp/land_w_200nm_buffer.grd",overwrite=T)

# running extraction ####
bathy_res_sd=raster("/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1_sd.grd")
modis_res_sd=raster("/Volumes/SeaGate/IUU_GRW/global_shp/modis_1_sd.grd")
bathy_res=raster("/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1.grd")
modis_res=raster("/Volumes/SeaGate/IUU_GRW/global_shp/modis_1.grd")
temp=raster("/Volumes/SeaGate/IUU_GRW/global_shp/sst.grd")
dist_shore=raster("/Volumes/SeaGate/IUU_GRW/global_shp/dist_shore.grd")
lat=raster("/Volumes/SeaGate/IUU_GRW/global_shp/lat.grd")
lon=raster("/Volumes/SeaGate/IUU_GRW/global_shp/lon.grd")

c=extracto(pnts=master_pnt,raster=bathy_res,colname="bathy")
c=extracto(pnts=c,raster=bathy_res_sd,colname="bathy_sd")
c=extracto(pnts=c,raster=modis_res,colname="modis")
c=extracto(pnts=c,raster=modis_res_sd,colname="modis_sd")
c=extracto(pnts=c,raster=temp,colname="temp")
c=extracto(pnts=c,raster=lat,colname="lat")
c=extracto(pnts=c,raster=lon,colname="lon")
c=extracto(pnts=c,raster=dist_shore,colname="dist_shore")
# b=b %>% .[complete.cases(.),]

master=c %>% mutate(random=1:nrow(c))

# model with covariates: bathy, chla, sd_bathy, sd_chla, distance from shore. TC=1 ####
# gap_brt_poiss = gbm.step(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random"),gbm.y="sum_gaps",family="poisson",learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 1, bag.fraction = 0.6)

dev_eval(gap_brt_poiss_fixed) ## 53.19

gap_brt = gbm.fixed(master,gbm.x=c("bathy","bathy_sd","modis","modis_sd"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 1, bag.fraction = 0.6)
dev_eval(gap_brt) # 48.19
summary(gap_brt)
gbm.plot(gap_brt)

## need to mask land
land_mask=bathy_res
land_mask[values(land_mask)>0]=NA
co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd) %>% mask(.,land_mask)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd")
gap_brt_ras=predict(co_stack,gap_brt,n.trees=2000,type="response")

plot(log(gap_brt_ras))

# model testing lat_lon, tc=3 ####
gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random","x","y"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) # 48.19
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)
# gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")

# ## add temp? --->>> seems to really f up predictions, removing for now
# nc=nc_open("/Volumes/SeaGate/IUU_GRW/global_shp/otemp.anal1deg.nc") 
# temp_ras=raster("/Volumes/SeaGate/IUU_GRW/global_shp/otemp.anal1deg.nc") %>% rotate()
# values(x)=1
# temp_ras=resample(temp_ras,x)
# c=extracto(pnts=c,raster=temp_ras,colname="temp")
# master=c %>% mutate(random=1:nrow(c))
# 
# gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random","x","y","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 1, bag.fraction = 0.6)
# dev_eval(gap_brt_poiss_fixed) # 54.46
# summary(gap_brt_poiss_fixed)
# gbm.plot(gap_brt_poiss_fixed)
# 
# gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("bathy","bathy_sd","modis","modis_sd","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 1, bag.fraction = 0.6)
# dev_eval(gap_brt_poiss_fixed) # 49.40
# summary(gap_brt_poiss_fixed)
# gbm.plot(gap_brt_poiss_fixed)
# 
# co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd,temp_ras) %>% mask(.,land_mask)
# names(co_stack)=c("bathy","bathy_sd","modis","modis_sd","temp")
# 
# gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
# plot(log(gap_brt_ras))
# plot(temp_ras)

# y, performed worse than random, refitting w.o.
gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("bathy","bathy_sd","modis","modis_sd","x","dist_shore"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) # 48.19
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

# need to create global distance to land and latitude rasters
data(wrld_simpl)
# Create a raster template for rasterizing the polys. 
# (set the desired grid resolution with res)
# r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, res=1)

# Rasterize and set land pixels to NA ---> use one of these to crease eez buffer [writen out raster, commented out code] ####
# r2 <- rasterize(wrld_simpl, x, 1)
# r3 <- mask(is.na(r2), r2, maskvalue=1, updatevalue=NA)

# Calculate distance to nearest non-NA pixel
# d <- distance(r2)
# writeRaster(d,"/Volumes/SeaGate/IUU_GRW/global_shp/dist_shore.grd",overwrite=T)
# dist_shore=d

## lat and long rasters
# lat <- lon <- x
# xy <- coordinates(x)
# lon[] <- xy[, 1]
# lat[] <- xy[, 2]
# writeRaster(lon,"/Volumes/SeaGate/IUU_GRW/global_shp/lon.grd",overwrite=T)
# writeRaster(lat,"/Volumes/SeaGate/IUU_GRW/global_shp/lat.grd",overwrite=T)
#####

## okay now predict model
co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd,lon,dist_shore) %>% mask(.,land_mask)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd","x","dist_shore")
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
plot(log(gap_brt_ras))

## trying a buffer
# buf=buffer(r2,width=370400)
# writeRaster(buf,"/Volumes/SeaGate/IUU_GRW/global_shp/land_w_200nm_buffer.grd",overwrite=T)

gap_brt_ras_mask=mask(gap_brt_ras,buf,inverse=T)
plot(log(gap_brt_ras_mask))

rasterVis::gplot(gap_brt_ras_mask)+geom_tile(aes(fill=log(value)))+
  # scale_fill_continuous(low="blue", high="darkred", 
  #                       guide="colorbar",na.value="white")+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# BRT, all vars including temp, TC=3 ####
gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random","x","y","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

# refitting w.o. y (worese than random)
gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","x","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd,lon,dist_shore,temp) %>% mask(.,buf,inverse=T)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd","x","dist_shore","temp")
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
plot(log(gap_brt_ras))

rasterVis::gplot(gap_brt_ras_mask)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

rasterVis::gplot(ras)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# masking modeled to real
model_mask=mask(gap_brt_ras_mask,ras)

modeled=rasterVis::gplot(model_mask)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),na.value="white",limits=c(2,17.3))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

real=rasterVis::gplot(ras)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),na.value="white",limits=c(2,17.3))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# BRT, all vars including temp, TC=3, this time including zeros ####
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
b=a %>% filter(frac_day_coverage>.2)
c=b %>% .[complete.cases(.),] %>%  dplyr::select(on_lat,on_lon)
c=c[,2:1]
x=raster()
res(x)=1
crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
ras=rasterize(c,x,fun=sum)
ras[is.na(ras[])]=0
ras=ras %>% mask(.,buf,inverse=T)
ras_pnt=rasterToPoints(ras) %>% as.data.frame()%>% rename(sum_gaps=layer)

# running extraction #
bathy_res_sd=raster("/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1_sd.grd")
modis_res_sd=raster("/Volumes/SeaGate/IUU_GRW/global_shp/modis_1_sd.grd")
bathy_res=raster("/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1.grd")
modis_res=raster("/Volumes/SeaGate/IUU_GRW/global_shp/modis_1.grd")
temp=raster("/Volumes/SeaGate/IUU_GRW/global_shp/sst.grd")
dist_shore=raster("/Volumes/SeaGate/IUU_GRW/global_shp/dist_shore.grd")
lat=raster("/Volumes/SeaGate/IUU_GRW/global_shp/lat.grd")
lon=raster("/Volumes/SeaGate/IUU_GRW/global_shp/lon.grd")

c=extracto(pnts=ras_pnt,raster=bathy_res,colname="bathy")
c=extracto(pnts=c,raster=bathy_res_sd,colname="bathy_sd")
c=extracto(pnts=c,raster=modis_res,colname="modis")
c=extracto(pnts=c,raster=modis_res_sd,colname="modis_sd")
c=extracto(pnts=c,raster=temp,colname="temp")
c=extracto(pnts=c,raster=lat,colname="lat")
c=extracto(pnts=c,raster=lon,colname="lon")
c=extracto(pnts=c,raster=dist_shore,colname="dist_shore")

master=c %>% mutate(random=1:nrow(c))

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random","x","y","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

# refitting w.o. y, x, temp, dist_shore (worese than random)
gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("bathy","bathy_sd","modis","modis_sd"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd) %>% mask(.,buf,inverse=T)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd")
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
plot(log(gap_brt_ras))

rasterVis::gplot(gap_brt_ras)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),limits=c(7.175879,18.3))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# masking modeled to real
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
b=a %>% filter(frac_day_coverage>.2)
c=b %>% .[complete.cases(.),] %>%  dplyr::select(on_lat,on_lon)
c=c[,2:1]
x=raster()
res(x)=1
crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
ras=rasterize(c,x,fun=sum)
model_mask=mask(gap_brt_ras,ras)

modeled=rasterVis::gplot(model_mask)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),na.value="white",limits=c(2,18.0))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

real=rasterVis::gplot(ras)+geom_tile(aes(fill=log(value)))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),na.value="white",limits=c(2,18.0))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 



# BRT, binomial


