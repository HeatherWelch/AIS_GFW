### sandkey diagram

library(dplyr)
library(networkD3)
library(tidyr)

# vignette from https://towardsdatascience.com/using-networkd3-in-r-to-create-simple-and-clear-sankey-diagrams-48f8ba8a4ace ####
refresults=read.csv("/Volumes/SeaGate/IUU_GRW/data/EU-referendum-result-data.csv")
results <- refresults %>% 
  dplyr::group_by(Region) %>% 
  dplyr::summarise(Remain = sum(Remain), Leave = sum(Leave))
results <- tidyr::gather(results, result, vote, -Region)
regions <- unique(as.character(results$Region))
nodes <- data.frame(node = c(0:13), 
                    name = c(regions, "Leave", "Remain"))

results <- merge(results, nodes, by.x = "Region", by.y = "name")
results <- merge(results, nodes, by.x = "result", by.y = "name")
links <- results[ , c("node.x", "node.y", "vote")]
colnames(links) <- c("source", "target", "value")

networkD3::sankeyNetwork(Links = links, Nodes = nodes, 
                         Source = 'source', 
                         Target = 'target', 
                         Value = 'value', 
                         NodeID = 'name',
                         units = 'votes')

# trying it out myself ####
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
## the above is code from the chord diagram, could probably streamline
nodes=unique(master[,1]) %>% as.data.frame() %>% mutate(node=seq(from=0,to=2)) %>% rename(name=Sovereign1.x)
nodes_end=unique(master[,1]) %>% as.data.frame() %>% mutate(node=seq(from=3,to=5)) %>% rename(name=Sovereign1.x)
nodes=rbind(nodes,nodes_end)
links= left_join(master,nodes,by=c("Sovereign1.x"="name")) %>% rename(source=node) %>% left_join(.,nodes_end,by=c("Sovereign1.y"="name"))%>% rename(target=node) %>% .[c(1:4,6:nrow(.)),]

networkD3::sankeyNetwork(Links = links, Nodes = nodes, 
                         Source = 'source', 
                         Target = 'target', 
                         Value = 'n', 
                         NodeID = 'name',
                         units = 'gaps')
