library(corrplot)
library(ggbiplot)
library(cluster)
source("load_functions.r")
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")

map.US <- map_data(map="state")
map.world = map_data(map="world")
shift <- 180+20
map.world2=map_data('world', wrap=c(0,360), ylim=c(-55,75))

### gap start points####
map=ggplot()+geom_map(data=map.US,map=map.US,aes(map_id=region,x=long,y=lat,fill="usa"),color="black")#+coord_cartesian()
map=map+geom_point(data=a,aes(x=lon,y=lat,fill="points"),alpha=.5,size=.7)
map


### gap end points? ####
map=ggplot()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat,fill="world"),color="black")
map=map+geom_point(data=a,aes(x=next_lon,y=next_lat,fill="points"),alpha=.5,size=.7)
map


### segments ####
c=a
c$on_lon.new <- ifelse(c$lon < -25, c$lon + 360, c$lon)
c$off_lon.new <- ifelse(c$next_lon < -25, c$next_lon + 360, c$next_lon)

map=ggplot()+geom_polygon(data = map.world2, aes(x=long, y = lat, group = group),color="white")+
  geom_segment(data=c,aes(x=on_lon.new,y=lat, xend=off_lon.new, yend=next_lat,color=flag))
map


### exploring the data ####
summary(a)
b=a %>% group_by(vessel_class) %>% summarise(n=n())
ggplot(a,aes(vessel_class))+geom_bar()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(a,aes(flag))+geom_bar()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))


### focues only of fishing boats ####
b=a %>% filter(on_fishing_list_best=="True") %>% .[complete.cases(.),] %>% filter(flag!="")
c=b %>% select(-c(ssvid,type,lat,lon,receiver_type,timestamp,next_timestamp,seg_id,next_seg_id,on_fishing_list_best))

## cluster ####
#https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
gower_dist <- daisy(c,
                    metric = "gower",
                    type = list(logratio = 3))

### gam ####
hist(log10(c$gap_km))
c$gap_km_log=log10(c$gap_km)
d=c %>% select(-c(next_lat,next_lon))

models2=mgcv::gam(gap_km~s(implied_speed_km_per_hour)+s(gap_hours)+s(distance_from_shore_m)+s(next_lat)+s(next_lon)+next_receiver_type+factor(vessel_class)+factor(flag),family="poisson",data=c)
models=mgcv::gam(gap_km~s(implied_speed_km_per_hour)+s(gap_hours)+s(distance_from_shore_m)+s(next_lat)+s(next_lon),family="poisson",data=c)

brt=gbm.fixed(c,gbm.x=c("implied_speed_km_per_hour","gap_hours","next_lat","distance_from_shore_m","next_lon","next_receiver_type","vessel_class","flag"),
             gbm.y="gap_km",family="gaussian",learning.rate = 0.00001, tree.complexity = 3, bag.fraction = 0.6,n.trees=1000)
              
brt2=gbm.fixed(c,gbm.x=c("implied_speed_km_per_hour","gap_hours","next_lat","distance_from_shore_m","next_receiver_type","vessel_class","flag"),
              gbm.y="gap_km",family="gaussian",learning.rate = 0.00001, tree.complexity = 3, bag.fraction = 0.6,n.trees=1000)
gbm.plot(brt2)
gbm.plot.fits(brt2)
