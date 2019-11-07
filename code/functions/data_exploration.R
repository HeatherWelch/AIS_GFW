## basic analysis of a GAP csv

# create four histograms: gap_hours, on_dist_from_port_m, frac_day_caverage, on_distance_from_shore_m
generic_histograms=function(data,outdir,saveName){
  hist1=ggplot()+geom_density(data=b,aes(x=gap_hours),fill="blue",alpha=.2)+theme_bw()
  hist2=ggplot()+geom_density(data=b,aes(x=on_distance_from_port_m),fill="blue",alpha=.2)+theme_bw()
  hist3=ggplot()+geom_density(data=b,aes(x=frac_day_coverage),fill="blue",alpha=.2)+theme_bw()
  hist4=ggplot()+geom_density(data=b,aes(x=on_distance_from_shore_m),fill="blue",alpha=.2)+theme_bw()
  cool=plot_grid(hist1,hist2,hist3,hist4,nrow = 1,ncol = 4)
  # return(cool)
  png(paste0(outdir,saveName,".png"),width=36,height=16,units='cm',res=400)
  par(ps=10)
  par(mar=c(4,4,1,1))
  par(cex=1)
  print({cool})
  dev.off()
}

# flag and class hists
flag_vesselClass_hist=function(data,outdir,saveName){
c=data %>% .[complete.cases(.),] %>% group_by(flag) %>% summarise(n=n())
flag=ggplot(c,aes(x=reorder(flag,-n),y=n))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
c=data %>% .[complete.cases(.),] %>% group_by(vessel_class) %>% summarise(n=n())
vessel=ggplot(c,aes(x=reorder(vessel_class,-n),y=n))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

cool=plot_grid(flag,vessel,nrow = 2,ncol = 1)
png(paste0(outdir,saveName,".png"),width=36,height=16,units='cm',res=400)
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({cool})
dev.off()
}

# global kernal density map
global_KD=function(data,outdir,saveName){
map.world = map_data(map="world")
map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=data)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+
  scale_fill_gradient(low = "blue", high = "red",guide = F)+theme_bw()+theme(legend.position = "none")

png(paste0(outdir,saveName,".png"),width=36,height=22,units='cm',res=400)
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({map})
dev.off()
}

# facet wrapped vessel class kernal density
vessel_facet_KD=function(data,outdir,saveName){
map.world = map_data(map="world")
map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=data)
map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+facet_wrap(~vessel_class)+
  scale_fill_gradient(low = "blue", high = "red")+theme_bw()+theme(legend.position = "none")

png(paste0(outdir,saveName,".png"),width=36,height=22,units='cm',res=400)
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({map})
dev.off()
}

# facet wrapped top 9 most gappy flags kernal density
flag_facet_KD=function(data,outdir,saveName){
  map.world = map_data(map="world")
  d=data %>% .[complete.cases(.),] %>% group_by(flag) %>% summarise(n=n()) %>% arrange(desc(n)) %>% .[1:9,1]
  toMatch=d$flag
  c=data %>% filter(flag %in% toMatch)
  map=ggplot()+stat_density_2d(aes(x=on_lon,y=on_lat,fill=stat(level),alpha=stat(level)),geom="polygon",bins=10,n=100,contour=T,data=c)
  map=map+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+facet_wrap(~flag)+
    scale_fill_gradient(low = "blue", high = "red")+theme_bw()+theme(legend.position = "none")
  
  png(paste0(outdir,saveName,".png"),width=36,height=22,units='cm',res=400)
  par(ps=10)
  par(mar=c(4,4,1,1))
  par(cex=1)
  print({map})
  dev.off()
}

# plot as abundance raster, use log=T to see more variation in abundance
abundance_raster=function(data,outdir,saveName,desired.resolution,log=T){
  map.world = map_data(map="world")
  c=data %>% .[complete.cases(.),] %>%  dplyr::select(on_lat,on_lon)
  c=c[,2:1]
  x=raster()
  res(x)=desired.resolution
  crs(x)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
  ras=rasterize(c,x,fun=sum)
  
  if(log==T) {
  map=rasterVis::gplot(ras)+geom_tile(aes(fill=log(value)))+
    scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),na.value="white")+
    coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) 
  }
  
  if(log==F) {
    map=rasterVis::gplot(ras)+geom_tile(aes(fill=value))+
      scale_fill_gradientn(colours = c("purple","blue","cyan","green","yellow","red"),na.value="white")+
      coord_equal()+geom_map(data=map.world,map=map.world,aes(map_id=region,x=long,y=lat))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      scale_y_continuous(expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0)) 
    
  }
  
  png(paste0(outdir,saveName,".png"),width=36,height=22,units='cm',res=400)
  par(ps=10)
  par(mar=c(4,4,1,1))
  par(cex=1)
  print({map})
  dev.off()
  
}

## examples of use: 
data=read.csv("/Volumes/SeaGate/IUU_GRW/data/raw_gaps_2018-10_2019/filtered_gaps_10_24_29.csv")
outdir="/Volumes/SeaGate/IUU_GRW/plots_testing_fcns/"
saveName="historgram"

generic_histograms(data,outdir = outdir,saveName = saveName)

saveName="historgram_flag_vessel"
flag_vesselClass_hist(data,outdir = outdir,saveName = saveName)

saveName="global_kd"
global_KD(data,outdir = outdir,saveName = saveName)

saveName="vessel_facet_kd"
vessel_facet_KD(data,outdir = outdir,saveName = saveName)

saveName="flag_facet_kd"
flag_facet_KD(data,outdir = outdir,saveName = saveName)

saveName="abundance_raster_4deg"
abundance_raster(data,outdir = outdir,saveName = saveName,desired.resolution=4)

saveName="abundance_raster_1deg"
abundance_raster(data,outdir = outdir,saveName = saveName,desired.resolution=1,log=F)
