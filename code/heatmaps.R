# heatmap. Flag x time, values = number of gaps
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")
master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),] %>% filter(flag!="")
master=master %>% dplyr::select(flag,timestamp)
library(lubridate)
test=master %>% mutate(month=month(timestamp),year=year(timestamp),ymd=ymd_hms(timestamp)) %>% mutate(ymd=date(ymd)) %>% group_by(flag,ymd) %>% .[complete.cases(.),]%>% summarise(n=n()) 
ggplot(test,aes(ymd,flag))+geom_tile(aes(fill=n))


master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),]
test=master %>% mutate(lat_r=round(lat,1),lon_r=round(lon,1)) %>% group_by(lat_r,lon_r)%>% summarise(n=n()) 
ggplot(test,aes(lon_r,lat_r))+geom_tile(aes(fill=n))


master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),]
test=master %>% mutate(lat_r=round(lat,0),ymd=ymd_hms(timestamp))%>% mutate(ymd=month(ymd))  %>% group_by(lat_r,ymd)%>% summarise(n=log(n())) 
ggplot(test,aes(ymd,lat_r))+geom_tile(aes(fill=n))+scale_fill_gradientn(colours = terrain.colors(10))


master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),]
test=master %>% mutate(gap_hr_r=log(round(gap_hours,1)),gap_km_r=log(round(gap_km,1)))%>% group_by(gap_hr_r,gap_km_r)%>% summarise(n=n()) %>% filter(gap_hr_r>=0)
ggplot(test,aes(gap_hr_r,gap_km_r))+geom_tile(aes(fill=n))+scale_fill_gradientn(colours = terrain.colors(10))


## clustering heatmap ####
master=a %>% filter(on_fishing_list_best=="True") %>% filter(type=="AIS.1"|type=="AIS.2"|type=="AIS.3") %>% filter(distance_from_shore_m>4000&next_distance_from_shore_m>4000) %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[complete.cases(.),]
rownames(master)=master$seg_id
test=master %>% dplyr::select(gap_hours, gap_km,lat,lon) %>% as.matrix() %>% scale()
test=master %>% dplyr::select(gap_hours, gap_km) %>% as.matrix() %>% scale()
heatmap(test,Colv = F,scale='row')

master=a %>% filter(on_fishing_list_best=="True") %>% filter(type=="AIS.1"|type=="AIS.2"|type=="AIS.3") %>% filter(distance_from_shore_m>4000&next_distance_from_shore_m>4000) %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[complete.cases(.),]
test=master %>% dplyr::select(flag,vessel_class,gap_hours,gap_km,lat,lon) %>% group_by(flag,vessel_class) %>% summarise(hours=sum(gap_hours),km=sum(gap_km),lat=mean(lat),lon=mean(lon))
test2=test %>% mutate(flag_class=paste0(flag,"_",vessel_class)) %>% dplyr::select(-c(flag,vessel_class)) %>% as.data.frame()
rownames(test2)=test2$flag_class
test3=test2[,2:ncol(test2)-1] %>% dplyr::select(-c(flag)) %>% as.matrix() %>% scale()
heatmap(test3,Colv = F,scale='row')

# install.packages("d3heatmap")
library(d3heatmap)
d3heatmap(test3, colors = "RdYlBu",
          k_row = 5, # Number of groups in rows
          k_col = 4 # Number of groups in columns
)


master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[complete.cases(.),]
test=master %>% dplyr::select(gap_hours,gap_km,lat,lon,distance_from_shore_m,next_lat,next_lon,next_distance_from_shore_m,vessel_class,flag)%>% group_by(flag,vessel_class) %>% 
  summarise(hours=mean(gap_hours),km=mean(gap_km),lat=mean(lat),lon=mean(lon),distance_from_shore_m=mean(distance_from_shore_m),next_lat=mean(next_lat),next_lon=mean(next_lon),next_distance_from_shore_m=mean(next_distance_from_shore_m))
test2=test %>% mutate(flag_class=paste0(flag,"_",vessel_class))%>% as.data.frame() %>% dplyr::select(-c(flag,vessel_class)) 
rownames(test2)=test2$flag_class
test3=test2[4:nrow(test2)-1,] %>% dplyr::select(-c(flag_class)) %>% as.matrix() %>% scale() %>% .[2:nrow(.),]
heatmap(test3,Colv = F,scale='none')

d3heatmap(test3, colors = "RdYlBu",
          k_row = 5, # Number of groups in rows
          k_col = 5 # Number of groups in columns
)


master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[complete.cases(.),]
test=master %>% dplyr::select(gap_hours,gap_km,distance_from_shore_m,next_distance_from_shore_m,vessel_class,flag)%>% group_by(flag,vessel_class) %>% 
  summarise(hours=mean(gap_hours),km=mean(gap_km),distance_from_shore_m=mean(distance_from_shore_m),next_distance_from_shore_m=mean(next_distance_from_shore_m))
test2=test %>% mutate(flag_class=paste0(flag,"_",vessel_class))%>% as.data.frame() %>% dplyr::select(-c(flag,vessel_class)) 
rownames(test2)=test2$flag_class
test3=test2[4:nrow(test2)-1,] %>% dplyr::select(-c(flag_class)) %>% as.matrix() %>% scale() %>% .[2:nrow(.),]
heatmap(test3,Colv = F,scale='none')

d3heatmap(test3, colors = "RdYlBu",
          k_row = 5, # Number of groups in rows
          k_col = 3 # Number of groups in columns
)



master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[complete.cases(.),]
test=master %>% dplyr::select(gap_hours,gap_km,lat,lon,distance_from_shore_m,next_lat,next_lon,next_distance_from_shore_m,vessel_class,flag,timestamp) %>% mutate(m=month(timestamp))%>% group_by(flag,vessel_class,m) %>% 
  summarise(hours=mean(gap_hours),km=mean(gap_km),lat=mean(lat),lon=mean(lon),distance_from_shore_m=mean(distance_from_shore_m),next_lat=mean(next_lat),next_lon=mean(next_lon),next_distance_from_shore_m=mean(next_distance_from_shore_m))
test2=test %>% mutate(flag_class=paste0(flag,"_",vessel_class,"_",m))%>% as.data.frame() %>% dplyr::select(-c(flag,vessel_class,m)) 
rownames(test2)=test2$flag_class
test3=test2[23:nrow(test2)-1,] %>% dplyr::select(-c(flag_class)) %>% as.matrix() %>% scale()
heatmap(test3,Colv = F,scale='none')

d3heatmap(test3, colors = "RdYlBu",
          k_row =8, # Number of groups in rows
          k_col = 5 # Number of groups in columns
)
