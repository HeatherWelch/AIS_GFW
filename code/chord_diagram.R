### chord diagram
# vignette from https://github.com/mattflor/chorddiag ####
devtools::install_github("mattflor/chorddiag", build_vignettes = TRUE)

m <- matrix(c(11975,  5871, 8916, 2868,
              1951, 10048, 2060, 6171,
              8010, 16145, 8090, 8045,
              1013,   990,  940, 6907),
            byrow = TRUE,
            nrow = 4, ncol = 4)
haircolors <- c("black", "blonde", "brown", "red")
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)
m


library(chorddiag)
groupColors <- c("#000000", "#FFDD89", "#957244", "#F26223")
chorddiag(m, groupColors = groupColors, groupnamePadding = 40)
# doing it for AIS Data 1. get EEZ file, 2. extract EEZ names for start and end points, 3. plot chord ####
library(tidyverse)
library(sf)
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")

# clean this up
master=a %>% filter(on_fishing_list_best=="True") %>% filter(type=="AIS.1"|type=="AIS.2"|type=="AIS.3") %>% filter(distance_from_shore_m>4000&next_distance_from_shore_m>4000) %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[0:500,]
coordinates(master)=~lon+lat
eez=st_read("/Volumes/SeaGate/IUU_GRW/data/World_EEZ_v10_20180221/eez_v10.shp")%>% st_set_crs(.,"+proj=longlat +datum=WGS84 +no_defs") %>% select(Sovereign1)
eezstart=st_as_sf(master) %>% st_set_crs(.,"+proj=longlat +datum=WGS84 +no_defs")
start_int=st_intersection(eezstart,eez)
start=start_int %>% as.data.frame() 

master=a %>% filter(on_fishing_list_best=="True") %>% filter(type=="AIS.1"|type=="AIS.2"|type=="AIS.3") %>% filter(distance_from_shore_m>4000&next_distance_from_shore_m>4000) %>% mutate(seg_id=substr(seg_id,11,20)) %>% .[0:500,]
coordinates(master)=~next_lon+next_lat
eezend=st_as_sf(master) %>% st_set_crs(.,"+proj=longlat +datum=WGS84 +no_defs")
end_int=st_intersection(eezend,eez)
end=end_int %>% as.data.frame()

master=left_join(start,end, by="seg_id") %>% select(Sovereign1.x,Sovereign1.y) %>% .[complete.cases(.),] %>% group_by(Sovereign1.x,Sovereign1.y) %>% summarise(n=n())
test=master %>% spread(Sovereign1.y,n) %>% .[,c(1,2,4,5)]
test[is.na(test)] <- 0
# t=as.matrix(test)
# 
# rownames(t)=colnames(t)[2:4]
# t=t[,2:4]
test=as.data.frame(test)
rownames(test)=colnames(test)
a=as.matrix(test)

# haircolors=colnames(test)[2:4]
# dimnames(t) <- list(have = haircolors,
#                     prefer = haircolors)

library(chorddiag)
# groupColors <- c("#000000", "#FFDD89", "#957244")
# chorddiag(a, groupColors = groupColors, groupnamePadding = 40)

library(circlize)
chordDiagram(a)
