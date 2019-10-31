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

# I don't have a good answer yet for what a good coverage threshold is and it probably depends on how long of gaps you're looking at (e.g. a 6 hour gap is more feasible than a 24 hour gap). However, the heatmap looks reasonable based on what we know about gaps (e.g. Argentina, Japan/Russia EEZ boundary) and anecdotes about where vessels turn of their AIS (e.g. leaving ports in Alaska). I suggest restricting satellite gaps to offshore (maybe start with 370400 meters = 200 nautical miles) - vessel density becomes an issue closer to shore for satellite reception.

### data were too big, ran this to filter ####
# a=read.csv("/Users/EcoCast/Downloads/raw_gaps_v20191021.csv")
# b=a %>% filter(on_type==off_type) %>% filter(on_type=="A")%>% filter(off_type=="A")  %>% filter(on_receiver_type=="satellite") %>% filter(off_receiver_type=="satellite") %>% filter(frac_day_coverage>.2)
# c=a %>% filter(on_type==off_type) %>% filter(on_type=="A")%>% filter(off_type=="A")  %>% filter(on_receiver_type=="satellite") %>% filter(off_receiver_type=="satellite") %>% filter(on_distance_from_shore_m>370400)
# write.csv(c,"filtered_gaps_10_24_29.csv")
#####

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

outfolder="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/"#;dir.create(outfolder)
rasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras/"#;dir.create(rasterdir)
pngdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_png/"#;dir.create(pngdir)

a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
b=a %>% filter(frac_day_coverage>.2)

#basic histograms ####
hist1=ggplot()+geom_density(data=b,aes(x=gap_hours),fill="blue",alpha=.2)+theme_bw()
hist2=ggplot()+geom_density(data=b,aes(x=on_distance_from_port_m),fill="blue",alpha=.2)+theme_bw()
hist3=ggplot()+geom_density(data=b,aes(x=frac_day_coverage),fill="blue",alpha=.2)+theme_bw()
hist4=ggplot()+geom_density(data=b,aes(x=on_distance_from_shore_m),fill="blue",alpha=.2)+theme_bw()
plot_grid(hist1,hist2,hist3,hist4,nrow = 1,ncol = 4)

df <- data.frame(Min=as.numeric(),
                 a1st_Qu=as.numeric(), 
                 Median=as.numeric(), 
                 Mean=as.numeric(), 
                 a3rd_Qu=as.numeric(), 
                 Max=as.numeric(), 
                 stringsAsFactors=FALSE) 
df[1,]=summary(b$gap_hours) %>% as.character() %>% as.numeric() %>% round(.,digits=3)
df[2,]=summary(b$on_distance_from_port_m) %>% as.character()%>% as.numeric() %>% round(.,digits=3)
df[3,]=summary(b$frac_day_coverage) %>% as.character()%>% as.numeric() %>% round(.,digits=3)
df[4,]=summary(b$on_distance_from_shore_m) %>% as.character()%>% as.numeric() %>% round(.,digits=3)

df$metric=c("gap_hours","on_distance_from_port_m","frac_day_coverage","on_distance_from_shore_m")
df


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

map=ggplot()+stat_density_2d(aes(x=off_lon,y=off_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=b)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+
  scale_fill_gradient(low = "blue", high = "red")+theme_bw()
map

# facet wrapping
c=b %>% filter(vessel_class!="dredge_fishing"|vessel_class!="other_fishing")
map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=c)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+facet_wrap(~vessel_class)+
  scale_fill_gradient(low = "blue", high = "red")+theme_bw()

png(paste0(outfolder,"heat_facet_vessel_class.png"),width=36,height=16,units='cm',res=400)
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
map
dev.off()

c=b %>% filter(flag=="CHN"|flag=="TWN"|flag=="ESP"|flag=="KOR"|flag=="PRT"|flag=="ECU"|flag=="FRA"|flag=="VUT"|flag=="JPN")
map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=c)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+facet_wrap(~flag)+
  +   scale_fill_gradient(low = "blue", high = "red")+theme_bw()

png(paste0(outfolder,"heat_facet_flag.png"),width=36,height=16,units='cm',res=400)
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
map
dev.off()

# digging into the 10 most gappy countries and gears ####
c=b %>% .[complete.cases(.),] %>% group_by(flag,vessel_class) %>% summarise(n=n()) %>% filter(flag=="CHN"|flag=="TWN"|flag=="ESP"|flag=="KOR"|flag=="PRT"|flag=="ECU"|flag=="FRA"|flag=="VUT"|flag=="JPN"|flag=="BLZ")
ggplot(c,aes(x=reorder(flag,-n),y=n))+geom_bar(aes(fill=vessel_class),stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

c=b %>% .[complete.cases(.),] %>% group_by(flag,vessel_class) %>% summarise(n=n()) %>% filter(vessel_class=="squid_jigger"|vessel_class=="drifting_longlines"|vessel_class=="trawlers"|vessel_class=="tuna_purse_seines"|vessel_class=="fishing")
ggplot(c,aes(x=reorder(vessel_class,-n),y=n))+geom_bar(aes(fill=flag),stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

# heatmaps to see which flag/vessels have gaps in similar areas ####
test=b %>% .[complete.cases(.),] %>%  dplyr::select(flag,vessel_class,on_lat,on_lon) %>% group_by(flag,vessel_class) %>% summarise(lat=mean(on_lat),lon=mean(on_lon))
test2=test %>% mutate(flag_class=paste0(flag,"_",vessel_class)) %>% dplyr::select(-c(flag,vessel_class)) %>% as.data.frame()
rownames(test2)=test2$flag_class
test3=test2[,2:ncol(test2)-1] %>% dplyr::select(-c(flag)) %>% as.matrix() %>% scale()
heatmap(test3,Colv = F,scale='row')

d3heatmap(test3, colors = "RdYlBu",
          k_row = 5, # Number of groups in rows
          k_col = 2 # Number of groups in columns
)

# cluster analysis of gear type and flag by location (lat, lon) ####
c=b %>% .[complete.cases(.),] %>%  dplyr::select(flag,vessel_class,on_lat,on_lon) %>% mutate(flag_class=paste0(flag,"_",vessel_class)) %>% dplyr::select(-c(flag,vessel_class)) 
d=c %>% group_by(flag_class) %>% summarise(n=n())
ggplot(d,aes(x=reorder(flag_class,-n),y=n))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1,size=7))
e=d %>% filter(n>20)
ggplot(e,aes(x=reorder(flag_class,-n),y=n))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1,size=7))

# creating rasters 1 dg ####
# blank raster
x=raster()
res(x)=1
crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

f=c %>% filter(flag_class=="CHN_squid_jigger") %>% dplyr::select(-flag_class)
f=f[,2:1]
ggplot(f,aes(x=on_lon,y=on_lat))+ geom_point(color="blue") + scale_colour_gradient(low="blue",high="red")+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))
# coordinates(f)=~on_lon+on_lat
g=rasterize(f,x)
plot(g)

rasterVis::gplot(g)+geom_tile(aes(fill=value))+
scale_fill_continuous(low="blue", high="darkred", 
                         guide="colorbar",na.value="white")+
coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 

# now batching
c=b %>%  dplyr::select(flag,vessel_class,on_lat,on_lon)%>% .[complete.cases(.),]  %>% mutate(flag_class=paste0(flag,"_",vessel_class)) %>% dplyr::select(-c(flag,vessel_class)) 
d=c %>% group_by(flag_class) %>% summarise(n=n())%>% filter(n>20)
FC_list=d$flag_class

for(i in 1:length(FC_list)){
  # for(i in 1:5){
  print(FC_list[i])
  a1=c %>% filter(flag_class==FC_list[i]) %>% dplyr::select(-flag_class)
  a1=a1[,2:1]
  g=rasterize(a1,x)
  writeRaster(g,paste0(rasterdir,FC_list[i]),overwrite=T)
  
  ras_plot=rasterVis::gplot(g)+geom_tile(aes(fill=value))+
    scale_fill_continuous(low="blue", high="darkred", 
                          guide="colorbar",na.value="white")+
    coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +ggtitle(paste0("1 degree gaps for ",FC_list[i]," (num gaps= ",d[i,2],")"))+theme(text = element_text(size=7),plot.title = element_text(size = 7))+ guides(shape = guide_legend(override.aes = list(size = 0.3)))
  
  png(paste0(pngdir,FC_list[i],".png"),width=16,height=8,units='cm',res=400)
  par(ps=10)
  par(mar=c(1,1,1,1))
  par(cex=1)
  print(ras_plot)
  dev.off()
  
}

# raster correlation

rasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras/"
a=list.files(rasterdir,pattern = ".grd",full.names = T)%>% stack()
names=list.files(rasterdir,pattern = ".grd") %>% gsub(".grd","",.)
names(a)=names

b=as.data.frame(a)
c=cor(b,method = 'pearson',use="pairwise.complete.obs")
c[is.na(c)] <- 0
corrplot(c, method="square",type = "upper",order="hclust",tl.col="black",tl.cex = .6)

# creating rasters 4 dg ####
# blank raster
x=raster()
res(x)=4
crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

# f=c %>% filter(flag_class=="CHN_squid_jigger") %>% dplyr::select(-flag_class)
# f=f[,2:1]
# ggplot(f,aes(x=on_lon,y=on_lat))+ geom_point(color="blue") + scale_colour_gradient(low="blue",high="red")+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))
# # coordinates(f)=~on_lon+on_lat
# g=rasterize(f,x)
# plot(g)
# 
# rasterVis::gplot(g)+geom_tile(aes(fill=value))+
#   scale_fill_continuous(low="blue", high="darkred", 
#                         guide="colorbar",na.value="white")+
#   coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   scale_y_continuous(expand = c(0,0)) + 
#   scale_x_continuous(expand = c(0,0)) 

# now batching
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
b=a %>% filter(frac_day_coverage>.2)

c=b %>%  dplyr::select(flag,vessel_class,on_lat,on_lon)%>% .[complete.cases(.),]  %>% mutate(flag_class=paste0(flag,"_",vessel_class)) %>% dplyr::select(-c(flag,vessel_class)) 
d=c %>% group_by(flag_class) %>% summarise(n=n())%>% filter(n>100)
FC_list=d$flag_class

newrasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras_4/";dir.create(newrasterdir)
newpngdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_png_4/";dir.create(newpngdir)


for(i in 1:length(FC_list)){
  # for(i in 1:5){
  print(FC_list[i])
  a1=c %>% filter(flag_class==FC_list[i]) %>% dplyr::select(-flag_class)
  a1=a1[,2:1]
  g=rasterize(a1,x)
  writeRaster(g,paste0(newrasterdir,FC_list[i]),overwrite=T)
  
  ras_plot=rasterVis::gplot(g)+geom_tile(aes(fill=value))+
    scale_fill_continuous(low="blue", high="darkred", 
                          guide="colorbar")+
    coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +ggtitle(paste0("1 degree gaps for ",FC_list[i]," (num gaps= ",d[i,2],")"))+theme(text = element_text(size=7),plot.title = element_text(size = 7))+ guides(shape = guide_legend(override.aes = list(size = 0.3)))
  
  png(paste0(newpngdir,FC_list[i],".png"),width=16,height=8,units='cm',res=400)
  par(ps=10)
  par(mar=c(1,1,1,1))
  par(cex=1)
  print(ras_plot)
  dev.off()
  
}

# raster correlation

rasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras_4/"
a=list.files(rasterdir,pattern = ".grd",full.names = T)%>% stack()
names=list.files(rasterdir,pattern = ".grd") %>% gsub(".grd","",.)
names(a)=names

b=as.data.frame(a)
c=cor(b,method = 'pearson',use="pairwise.complete.obs")
c[is.na(c)] <- 0
corrplot(c, method="color",type = "upper",order="hclust",tl.col="black",tl.cex = .6,hclust.method="ward.D")





# unsupervised classification 4 dg ####
# starting with 4 deg and NAs, might switch NAs to zero as this is likely swamping correlations
# d.n. work, switching to zeros
# https://rspatial.org/raster/rs/4-unsupclassification.html
# http://remote-sensing.eu/unsupervised-classification-with-r/
# http://www.sthda.com/english/wiki/print.php?id=234#clustering-validation

rasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras_4/"
a=list.files(rasterdir,pattern = ".grd",full.names = T)%>% stack()
names=list.files(rasterdir,pattern = ".grd") %>% gsub(".grd","",.)
names(a)=names
v=getValues(a)
v[is.na(v)] <- 0
i<- which(!is.na(v))
v <- na.omit(v)
E <- kmeans(v, 10, iter.max = 100, nstart = 10)
kmeans_raster <- raster(a)
kmeans_raster[i] <- E$cluster
plot(kmeans_raster)
kmeans_raster[values(kmeans_raster)==6]<-NA

ras_plot=rasterVis::gplot(kmeans_raster)+geom_tile(aes(fill=factor(value)))+
  # scale_fill_manual("",values=c("1"="#fef65b","2"="#ff0000", "3"="#daa520","4"="#0000ff","5"="#0000ff","7"="#00ff00","8"="#cbbeb5","9"="#c3ff5b", "10"="#ff7373"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +ggtitle(paste0("10 clusters"))+theme(text = element_text(size=7),plot.title = element_text(size = 7))+ guides(shape = guide_legend(override.aes = list(size = 0.3)))

png(paste0(outfolder,"clusters_10.png"),width=16,height=8,units='cm',res=400)
par(ps=10)
par(mar=c(1,1,1,1))
par(cex=1)
ras_plot
dev.off()

# cluster validation
wss <- 0

# For 1 to 15 cluster centers
for (i in 1:15) {
  km.out <- kmeans(v, centers = i, nstart = 20)
  # Save total within sum of squares to wss variable
  wss[i] <- km.out$tot.withinss
}

# Plot total within sum of squares vs. number of clusters
plot(1:15, wss, type = "b", 
     xlab = "Number of Clusters", 
     ylab = "Within groups sum of squares")

# 6 clusters loooks good
v=getValues(a)
v[is.na(v)] <- 0
i<- which(!is.na(v))
v <- na.omit(v)
E <- kmeans(v, 6, iter.max = 100, nstart = 10)
kmeans_raster <- raster(a)
kmeans_raster[i] <- E$cluster
plot(kmeans_raster)
writeRaster(kmeans_raster,"/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/kmeans_raster_6.grd")
# kmeans_raster[values(kmeans_raster)==6]<-NA

ras_plot=rasterVis::gplot(kmeans_raster)+geom_tile(aes(fill=factor(value)))+
  # scale_fill_manual("",values=c("1"="#fef65b","2"="#ff0000", "3"="#daa520","4"="#0000ff","5"="#0000ff","7"="#00ff00","8"="#cbbeb5","9"="#c3ff5b", "10"="#ff7373"))+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +ggtitle(paste0("6 clusters"))+theme(text = element_text(size=7),plot.title = element_text(size = 7))+ guides(shape = guide_legend(override.aes = list(size = 0.3)))

png(paste0(outfolder,"clusters_6.png"),width=16,height=8,units='cm',res=400)
par(ps=10)
par(mar=c(1,1,1,1))
par(cex=1)
ras_plot
dev.off()

# visualizing clusters
Emeans=E$center 
heatmap(Emeans,scale="column")
Emeans=E$center %>% scale()
mapp=d3heatmap(Emeans, colors = "RdYlBu",
          k_row = 5, # Number of groups in rows
          k_col = 2 # Number of groups in columns
)

rows=mapp$x$matrix$rows
columns=mapp$x$matrix$cols
data=mapp$x$matrix$data %>% as.vector()
# dataT=t(data)

test=matrix(data,nrow=6,ncol=39,byrow = T)
# test=as.data.frame(as.vector(data),nrow=6,ncol=39)
dimnames(test)=list(rows,columns)

# option 2. assign each gear_flag to clusters where it has value >.6 (determined by lookign at what makes sense in heat map)
test2=test
test2[test2>.6]=1
test2[test2<.6]=0

# option 2. assign each gear_flag to the cluster where it has the highest value
# test[which(apply(test, 2, function(x) x == max(x,na.rm=TRUE)))] <- 10
# test[test<10]=0
class(test2)="numeric"
d3heatmap(test2, colors = "RdYlBu",
               k_row = 5, # Number of groups in rows
               k_col = 2 # Number of groups in columns
)
par(mar = c(10, 2, 2, 2))
heatmap.2(test2,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
colsums=apply(test2,1,sum)
names(test2)[which(test2[1,] == 1, arr.ind=T)[, "col"]]
# lapply(apply(test2, 1, function(x)which(x=1)), names)
test2T=t(test2) %>% as.data.frame() %>% mutate(vessel_class=rownames(.))
colnames(test2T)=c("c3","c4","c6","c1","c5","c2","vessel_class")
rasterdir="/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/gear_flag_ras_4/"
kmeans_raster=raster("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/kmeans_raster_6.grd")

for(i in 1:6){
  if(i!=2){
    col=paste0("c",i)
    colpos=grep(col,names(test2T))
    c1=test2T %>% filter(.[,colpos]>0)
    toMatch=c1$vessel_class
    print(col)
    print(toMatch)
    ras=lapply(toMatch,function(x)paste0(rasterdir,x,".grd")) %>% unlist() %>% stack() %>% sum(na.rm=T)
    ras[values(ras)==0]<-NA 
    clus_ras=kmeans_raster
    clus_ras[values(clus_ras)!=i]<-NA 
    clus_ras2=clus_ras %>% rasterToPolygons()
    cr=fortify(clus_ras2)

ras_plot=rasterVis::gplot(ras)+geom_tile(aes(fill=value))+
  scale_fill_continuous(low="blue", high="darkred", 
                        guide="colorbar",na.value="white")+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +ggtitle(paste0(col))+theme(text = element_text(size=7),plot.title = element_text(size = 7))+ guides(shape = guide_legend(override.aes = list(size = 0.3)))+
  geom_polygon(data=cr,aes(x=long,y=lat,group=id),color="yellow",fill="yellow",alpha=.2)


png(paste0(outfolder,"cluster_",col,".png"),width=16,height=8,units='cm',res=400)
par(ps=10)
par(mar=c(1,1,1,1))
par(cex=1)
print(ras_plot)
dev.off()
  }
}

# wow that was exhausting. going back to bigger picture ####
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
b=a %>% filter(frac_day_coverage>.2)
c=b %>% .[complete.cases(.),] %>%  dplyr::select(on_lat,on_lon)
c=c[,2:1]
x=raster()
res(x)=1
crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
ras=rasterize(c,x)
writeRaster(ras,"/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/plots_10_24_10/total_gaps.tif",format="GTiff")
# writeRaster(ras,"/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/raw_gaps_v20191021.csv/total_gaps.grd")

# fname<-system.file("/Volumes/SeaGate/IUU_GRW/global_shp/WDPA_Oct2019-shapefile/WDPA_Oct2019-shapefile-polygons.shp",package = "sf")
# nc <- st_read("/Volumes/SeaGate/IUU_GRW/global_shp/WDPA_Oct2019-shapefile/WDPA_Oct2019-shapefile-polygons.shp")
# 
# nc1=nc %>% filter(IUCN_CAT=="Ia"|IUCN_CAT=="Ib"|IUCN_CAT=="II"|IUCN_CAT=="III"|IUCN_CAT=="IV") %>% filter(STATUS=="Designated") %>% filter(MARINE==2)
# st_write(nc1,"/Volumes/SeaGate/IUU_GRW/global_shp/filtered_mpas.shp")
# nc1=st_read(("/Volumes/SeaGate/IUU_GRW/global_shp/filtered_mpas.shp"))
nc1=readShapeSpatial(("/Volumes/SeaGate/IUU_GRW/global_shp/filtered_mpas.shp"))
nc2=fortify(nc1)
nc3=nc2 %>% filter(id==0|id==1)

# nc3=nc1[1:3,]
# simplepolys <- rmapshaper::ms_simplify(input = as(nc1, 'Spatial')) %>%
#   st_as_sf()


# ncF=fortify(nc1)
# broken attempts at sf
# nc1=st_read(("/Volumes/SeaGate/IUU_GRW/global_shp/filtered_mpas.shp"))
# nc3=nc1[1:3,]
# ggplot()+geom_sf()+geom_sf(data=nc3,fill="yellow")


ras_plot=rasterVis::gplot(ras)+geom_tile(aes(fill=value))+
  scale_fill_distiller(palette="Spectral",na.value="white")+
  coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +ggtitle(paste0("Total gaps"))+theme(text = element_text(size=7),plot.title = element_text(size = 7))+ guides(shape = guide_legend(override.aes = list(size = 0.3)))+
 geom_polygon(data=nc2,aes(x=long,y=lat,group=id),color="yellow",fill="yellow",alpha=.2)


png(paste0(outfolder,"mpa_total_gaps.png"),width=16,height=8,units='cm',res=400)
par(ps=10)
par(mar=c(1,1,1,1))
par(cex=1)
ras_plot
dev.off()

a=raster("/Users/heatherwelch/Dropbox/JPSS/global/global_averages_mask/modis.grd")
writeRaster(a,"/Users/heatherwelch/Desktop/Win7Desktop/GFW_IUU/modis.tif",format="GTiff")