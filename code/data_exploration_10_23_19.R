### new data following GFW workshop. Email from Tyler: 
# Here's a link to download an updated raw dataset of all AIS gap events longer than six hours since the beginning of 2018.
# https://drive.google.com/a/globalfishingwatch.org/file/d/1jT1zLqZIF2UtpctBvDMHimVjlRb6koW-/view?usp=sharing 
# 
# The fields should mostly be self-explanatory, but the following are important to note:
# off/on_type: The AIS device type (Class A/B). The off_type and on_type fields should almost always match, but not always. For now, probably best to filter out any events where these don't match.
# off/on_receiver_type: What receiver type (satellite/terrestrial) received the AIS message.
# frac_day_coverage: On average, what fraction of 5 minute windows in a day are AIS messages received from the quarter degree grid cell (off_lat/lon_bin) where the off event occurred? This is our first approximation of the reception quality in the area and a rough proxy for the likelihood of seeing a gap in a cell (e.g. higher coverage, less likely a gap should occur).
# vessel_class: geartype
# flag: ISO3 country code for the vessel
# The data file is rather large (~2 GB uncompressed). My suggestion is to start by subsetting to gaps where the off event was broadcast by a Class A device and received by a satellite receiver in an area with relatively high coverage.


### data were too big, ran this to filter ####
# a=read.csv("/Users/EcoCast/Downloads/raw_gaps_v20191021.csv")
# b=a %>% filter(on_type==off_type) %>% filter(on_type=="A")%>% filter(off_type=="A")  %>% filter(on_receiver_type=="satellite") %>% filter(off_receiver_type=="satellite") %>% filter(frac_day_coverage>.2)
# c=a %>% filter(on_type==off_type) %>% filter(on_type=="A")%>% filter(off_type=="A")  %>% filter(on_receiver_type=="satellite") %>% filter(off_receiver_type=="satellite") #%>% filter(frac_day_coverage>.2)
# write.csv(c,"filtered_gaps.csv")
#####

library(tidyverse)
library(MASS)
library(raster)
library(ggmap)

a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps.csv")
b=a %>% filter(frac_day_coverage>.2)

#basic histograms ####
hist(b$gap_hours)
summary(b$gap_hours)

hist(b$on_distance_from_port_m)
summary(b$on_distance_from_port_m)

hist(b$frac_day_coverage)
summary(b$frac_day_coverage)

hist(b$on_distance_from_shore_m)
summary(b$on_distance_from_shore_m)

# flag and class hists ####
c=b %>% .[complete.cases(.),] %>% group_by(flag) %>% summarise(n=n())
ggplot(c,aes(x=reorder(flag,-n),y=n))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

c=b %>% .[complete.cases(.),] %>% group_by(vessel_class) %>% summarise(n=n())
ggplot(c,aes(x=reorder(vessel_class,-n),y=n))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

# kernal density maps ####
map.world = map_data(map="world")

ggplot(b,aes(x=on_lon,y=on_lat))+ geom_point(color="blue") + scale_colour_gradient(low="blue",high="red")+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))
ggplot(b,aes(x=off_lon,y=off_lat))+ geom_point(color="blue") + scale_colour_gradient(low="blue",high="red")+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))

map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=b)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+
  scale_fill_gradient(low = "blue", high = "red")+theme_bw()
map

# facet wrapping
map=ggplot()+stat_density_2d(aes(x=off_lon,y=off_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=b)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+facet_wrap(~vessel_class)+
  scale_fill_gradient(low = "blue", high = "red")+theme_bw()
map

c=b %>% filter(flag=="CHN"|flag=="ESP"|flag=="RUS"|flag=="TWN"|flag=="USA"|flag=="ARG"|flag=="NOR"|flag=="KOR"|flag=="PRT")
map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=c)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+facet_wrap(~flag)+
  scale_fill_gradient(low = "blue", high = "red")+theme_bw()
map
    