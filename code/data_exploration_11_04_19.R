## attempting to model clusters ahead of GFW meeting 11.06.19

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
  a=raster::extract(raster,b)
  pnts$new=a
  position=grep("new",names(pnts))
  colnames(pnts)[position]=colname
  return(pnts)
}

# read in rasters ####
buf=raster("/Volumes/SeaGate/IUU_GRW/global_shp/land_w_200nm_buffer.grd")
bathy_res_sd=raster("/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1_sd.grd")
modis_res_sd=raster("/Volumes/SeaGate/IUU_GRW/global_shp/modis_1_sd.grd")
bathy_res=raster("/Volumes/SeaGate/IUU_GRW/global_shp/bathy_1.grd")
modis_res=raster("/Volumes/SeaGate/IUU_GRW/global_shp/modis_1.grd")
temp=raster("/Volumes/SeaGate/IUU_GRW/global_shp/sst.grd")
dist_shore=raster("/Volumes/SeaGate/IUU_GRW/global_shp/dist_shore.grd")
lat=raster("/Volumes/SeaGate/IUU_GRW/global_shp/lat.grd")
lon=raster("/Volumes/SeaGate/IUU_GRW/global_shp/lon.grd")

rasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras/"

# cluster 3 ####

toMatch=c("ALB_squid_jigger","VUT_squid_jigger","KOR_Squid_jigger","RUS_fishing","FLK_trawlers","CHN_trawlers","CHN_fishing","KOR_trawlers","CHN_tuna_purse_seines","TWN_squid_jigger","KOR_set_longlines","COL_tuna_purse_seines","CHN_squid_jigger","ECU_tuna_purse_seines","CRI_drifting_longlines","ESP_trawlers")
ras=lapply(toMatch,function(x)paste0(rasterdir,x,".grd")) %>% unlist() %>% stack() %>% sum(na.rm=T)
ras[values(ras)==0]<-NA
ras=ras %>% mask(.,buf,inverse=T)
ras_pnt=rasterToPoints(ras) %>% as.data.frame()%>% rename(sum_gaps=layer)

c=extracto(pnts=ras_pnt,raster=bathy_res,colname="bathy")
c=extracto(pnts=c,raster=bathy_res_sd,colname="bathy_sd")
c=extracto(pnts=c,raster=modis_res,colname="modis")
c=extracto(pnts=c,raster=modis_res_sd,colname="modis_sd")
c=extracto(pnts=c,raster=temp,colname="temp")
c=extracto(pnts=c,raster=lat,colname="lat")
c=extracto(pnts=c,raster=lon,colname="lon")
c=extracto(pnts=c,raster=dist_shore,colname="dist_shore")

master=c %>% mutate(random=sample(1:nrow(.),nrow(.),replace=F))

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","x","y","temp","random"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","x","y","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)

co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd,lon,dist_shore,temp,lat) %>% mask(.,buf,inverse=T)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd","x","dist_shore","temp","y")
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
plot(gap_brt_ras)

rasterVis::gplot(gap_brt_ras)+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# cluster 4 ####

toMatch=c("ECU_drifting_longlines","JPN_drifting_longlines","ESP_trawlers","PRT_set_longlines","PRT_trawlers","EST_trawlers","ESP_other_purse_seines","ESP_set_gillnets","PRT_drifting_longlines","ESP_drifting_longlines","BLZ_tuna_purse_seines","CUW_tuna_purse_seines")
ras=lapply(toMatch,function(x)paste0(rasterdir,x,".grd")) %>% unlist() %>% stack() %>% sum(na.rm=T)
ras[values(ras)==0]<-NA
ras=ras %>% mask(.,buf,inverse=T)
ras_pnt=rasterToPoints(ras) %>% as.data.frame()%>% rename(sum_gaps=layer)

c=extracto(pnts=ras_pnt,raster=bathy_res,colname="bathy")
c=extracto(pnts=c,raster=bathy_res_sd,colname="bathy_sd")
c=extracto(pnts=c,raster=modis_res,colname="modis")
c=extracto(pnts=c,raster=modis_res_sd,colname="modis_sd")
c=extracto(pnts=c,raster=temp,colname="temp")
c=extracto(pnts=c,raster=lat,colname="lat")
c=extracto(pnts=c,raster=lon,colname="lon")
c=extracto(pnts=c,raster=dist_shore,colname="dist_shore")

master=c %>% mutate(random=sample(1:nrow(.),nrow(.),replace=F))

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","random","x","y","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","x","y","temp"),gbm.y="sum_gaps",family="poisson",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
plot(gap_brt_ras)

rasterVis::gplot(gap_brt_ras)+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# point-based models - zeros come from 1 degree bins wi.thout points ####
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
ras_pnt_abs=ras_pnt %>% filter(sum_gaps==0) %>% mutate(presAbs=0) %>% select(-c(sum_gaps)) %>% .[sample(nrow(.),10000),] ## take a random sample, was doing too good with all points
pres=a %>% select(on_lat,on_lon) %>% mutate(presAbs=1) %>% rename(x=on_lon) %>% rename(y=on_lat)

master=rbind(pres,ras_pnt_abs) %>% mutate(random=sample(1:nrow(.),nrow(.),replace=F))

c=extracto(pnts=master,raster=bathy_res,colname="bathy")
c=extracto(pnts=c,raster=bathy_res_sd,colname="bathy_sd")
c=extracto(pnts=c,raster=modis_res,colname="modis")
c=extracto(pnts=c,raster=modis_res_sd,colname="modis_sd")
c=extracto(pnts=c,raster=temp,colname="temp")
c=extracto(pnts=c,raster=lat,colname="lat")
c=extracto(pnts=c,raster=lon,colname="lon")
c=extracto(pnts=c,raster=dist_shore,colname="dist_shore")

master=c
# write.csv(master,"/Users/heatherwelch/Desktop/Win7Desktop/GFW_IUU/binary.csv")

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","x","y","temp","random"),gbm.y="presAbs",family="bernoulli",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd,lon,dist_shore,temp,lat) %>% mask(.,buf,inverse=T)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd","x","dist_shore","temp","y")
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=2000,type="response")
plot(gap_brt_ras)

rasterVis::gplot(gap_brt_ras)+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# refitting without lat and long

gap_brt_poiss_fixed = gbm.fixed(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","temp"),gbm.y="presAbs",family="bernoulli",n.trees=2000,learning.rate = 0.01, tree.complexity = 3, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)

gap_brt_poiss_fixed = gbm.step(master,gbm.x=c("dist_shore","bathy","bathy_sd","modis","modis_sd","temp"),gbm.y="presAbs",family="bernoulli",learning.rate = 0.1, tree.complexity = 5, bag.fraction = 0.6)
dev_eval(gap_brt_poiss_fixed) 
summary(gap_brt_poiss_fixed)
gbm.plot(gap_brt_poiss_fixed)
write_rds(gap_brt_poiss_fixed,"/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/brt_bernoulli_lr.1,tc5,bf.6.rds")


co_stack=stack(bathy_res,bathy_res_sd,modis_res,modis_res_sd,dist_shore,temp) %>% mask(.,buf,inverse=T)
names(co_stack)=c("bathy","bathy_sd","modis","modis_sd","dist_shore","temp")
gap_brt_ras=predict(co_stack,gap_brt_poiss_fixed,n.trees=gap_brt_poiss_fixed$gbm.call$best.trees,type="response")
plot(gap_brt_ras)

rasterVis::gplot(gap_brt_ras)+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# attempting some validation

kfolds_eval <- function(dataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=10,ncol=4))
  colnames(Evaluations_kfold) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family="bernoulli", tree.complexity=tc,
                                 learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$PresAbs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    e <- evaluate(p=pres, a=abs)
    Evaluations_kfold[counter,1] <- k
    Evaluations_kfold[counter,2] <- dev
    Evaluations_kfold[counter,3] <- e@auc
    Evaluations_kfold[counter,4] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_kfold)}

test=kfolds_eval(master,gbm.x =c("dist_shore","bathy","bathy_sd","modis","modis_sd","temp"),gbm.y="presAbs",lr=0.1,tc=5 )
