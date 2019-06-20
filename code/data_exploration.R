### script exploring california data

library(corrplot)
library(ggbiplot)
library(cluster)
source("load_functions.r")
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/AIS_gaps_California.csv")
names(a)

map.US <- map_data(map="state")
map.world = map_data(map="world")
shift <- 180+20
# test=map.world %>% mutate(long.new=long+shift) %>% mutate(long.new,ifelse(long.new > 180, long.new-360, long.new)) %>% 
#   mutate(to.split=sum(diff(long.new) > 300, na.rm=TRUE) > 0, by=group) %>% 
#   mutate(gr.split=ifelse(to.split & long.new < 0, paste0(group, ".", 1), group))

map.world2=map_data('world', wrap=c(0,360), ylim=c(-55,75))

b=read.csv("/Volumes/SeaGate/IUU_GRW/data/AIS_gaps_California.csv")

coordinates(a)=~on_lon+on_lat

### gap start points####
map=ggplot()+geom_map(data=map.US,map=map.US,aes(map_id=region,x=long,y=lat,fill="usa"),color="black")#+coord_cartesian()
map=map+geom_point(data=b,aes(x=on_lon,y=on_lat,fill="points"),alpha=.5,size=.7)
map


### gap end points ####
map=ggplot()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat,fill="world"),color="black")#+coord_cartesian()
map=map+geom_point(data=b,aes(x=off_lon,y=off_lat,fill="points"),alpha=.5,size=.7)
map

## segments ####
# c=b %>% mutate(on_lon.new=on_lon+shift) %>% mutate(on_lon.new=ifelse(on_lon.new > 180, on_lon.new-360, on_lon.new)) %>% 
#   mutate(off_lon.new=off_lon+shift) %>% mutate(off_lon.new=ifelse(off_lon.new > 180, off_lon.new-360, off_lon.new))

c=b
c$on_lon.new <- ifelse(c$on_lon < -25, c$on_lon + 360, c$on_lon)
c$off_lon.new <- ifelse(c$off_lon < -25, c$off_lon + 360, c$off_lon)

# map=ggplot() + 
#   geom_polygon(data=test, 
#                aes(x=long.new, y=lat, group=gr.split), 
#                colour="black", fill="gray80", size = 0.25) +
#   coord_equal()
# map=map+geom_segment(data=c,aes(x=on_lon.new,y=on_lat, xend=off_lon.new, yend=off_lat))
# map


map=ggplot()+geom_polygon(data = map.world2, aes(x=long, y = lat, group = group),color="white")+
  geom_segment(data=c,aes(x=on_lon.new,y=on_lat, xend=off_lon.new, yend=off_lat,color=flag))
  

# map=ggplot()+geom_polygon(data = map.world2,map=map.world2, aes(map_id=region,x=long, y = lat,fill="world"),color="white")+
#   geom_segment(data=c,aes(x=on_lon.new,y=on_lat, xend=off_lon.new, yend=off_lat,color=flag))
map


### histograms ###
days_low=b %>% select(gap_days) %>% filter(gap_days<10)
days_high=b %>% select(gap_days) %>% filter(gap_days>=10&gap_days<700)
hist(days_low$gap_days)
hist(days_high$gap_days)

km_low=b %>% select(gap_km) %>% filter(gap_km<18)
km_high=b %>% select(gap_km) %>% filter(gap_km>=18&gap_km<1000)
hist(km_low$gap_km)
hist(km_high$gap_km)

### scatter ###
plot(b$gap_days,b$gap_km)

### corrplot ###
e=b %>% select(-c(flag,vessel_class,rank,on_start,off_start,off_type,on_type)) %>% .[complete.cases(.),]
m=cor(e)
corrplot(m,type="upper",order = "hclust",method="square")

### pca ###
e=b %>% select(-c(flag,vessel_class,rank,on_start,off_start,off_type,on_type,off_lat_bin,on_lat_bin,off_class_a_pos_per_day,on_class_a_pos_per_day)) %>% .[complete.cases(.),]
m=prcomp(e,center=T,scale.=T)
summary(m)
flag=b %>% select(-c(vessel_class,rank,on_start,off_start,off_type,on_type,off_lat_bin,on_lat_bin,off_class_a_pos_per_day,on_class_a_pos_per_day)) %>% .[complete.cases(.),] %>% select(flag) 
flag=flag$flag%>% as.factor()

vessel_class=b %>% select(-c(flag,rank,on_start,off_start,off_type,on_type,off_lat_bin,on_lat_bin,off_class_a_pos_per_day,on_class_a_pos_per_day)) %>% .[complete.cases(.),] %>% select(vessel_class) 
vessel_class=vessel_class$vessel_class%>% as.factor()
ggbiplot(m,groups=flag,varname.size=5,circle=T,ellipse=T)
ggbiplot(m,groups=vessel_class,varname.size=5,circle=T,ellipse=T)
# ggplot2::autoplot(prcomp(e,center=T,scale.=T))

### cluster ####
# https://towardsdatascience.com/hierarchical-clustering-on-categorical-data-in-r-a27e578f2995
e=b %>% select(-c(rank,on_start,off_start,off_type,on_type,off_lat_bin,on_lat_bin,off_class_a_pos_per_day,on_class_a_pos_per_day,mmsi,off_distance_from_shore,on_distance_from_shore,off_distance_from_port,on_distance_from_port,off_lon_bin,on_lon_bin)) %>% .[complete.cases(.),]
ee=e[1:20000,]
d=daisy(ee,metric = c("gower"))

aggl.clust.c <- hclust(d, method = "complete")
plot(aggl.clust.c,
     main = "Agglomerative, complete linkages")

### removing vessels at port ----> ask if there's a best way to remove ports ####
e=b %>% select(off_avg_distance_port_last_8_hrs) %>% filter(off_avg_distance_port_last_8_hrs<.5)
hist(e$off_avg_distance_port_last_8_hrs) ### filter off_ave >.5

e=b %>% select(on_avg_distance_port_next_8_hrs) %>% filter(on_avg_distance_port_next_8_hrs<.5)
hist(e$on_avg_distance_port_next_8_hrs) ### filter off_ave >.5

d=b %>% filter(off_avg_distance_port_last_8_hrs>.5,on_avg_distance_port_next_8_hrs>.5,off_distance_from_port >.5,on_distance_from_port>.5)%>% .[complete.cases(.),]

### also removing vessels that went nowhere in zerodays
plot(d$gap_days,d$gap_km)
d=d %>% filter(gap_days>0&&gap_km>1)
d=d %>% filter(gap_days>2&&gap_km<200) ### this logic is messed up



e=d %>% select(-c(rank,on_start,off_start,off_type,on_type,off_lat_bin,on_lat_bin,off_class_a_pos_per_day,on_class_a_pos_per_day,mmsi,off_distance_from_shore,on_distance_from_shore,off_distance_from_port,on_distance_from_port,off_lon_bin,on_lon_bin)) %>% .[complete.cases(.),]

d=daisy(e,metric = c("gower"))

aggl.clust.c <- hclust(d, method = "complete")
plot(aggl.clust.c,
     main = "Agglomerative, complete linkages")

