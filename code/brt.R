#### brts 2018 data
# 1 fishing only
# 2 Class A only
# Remove port

library(corrplot)
library(ggbiplot)
library(cluster)
source("load_functions.r")
library(gbm)
library(dismo)

a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")

# clean this up
master=a %>% filter(on_fishing_list_best=="True") %>% filter(type=="AIS.1"|type=="AIS.2"|type=="AIS.3") %>% filter(distance_from_shore_m>4000&next_distance_from_shore_m>4000) %>% mutate(seg_id=substr(seg_id,11,20))

### kernal density for ansences ####
# Convert xy points to SpatialPoints class:
xy <- master[, c('lat', 'lon')]
xy <- xy[xy$lat>=30 & xy$lat<=48,] 
coordinates(xy)=~lon+lat
#xy.sp <- sp::SpatialPoints(xy, proj4string=CRS("+proj=longlat +ellps=WGS84"))
# xy.sp <- sp::SpatialPoints(xy)

# Calculate 95% UD kernel for home range boundary
UD<-adehabitatHR::kernelUD(xy, kern="bivnorm") # Creates bivariate normal utilization distribution
image(UD) #visualize
UD95<-adehabitatHR::getverticeshr(UD, 95.9) # Pulls out 95% kernel boundary, or whichever percentage you enter

# Create random points
rpts<-sp::spsample(UD95, 1842, type="random") #can also type="regular", "hexogonal", "clustered", etc.
#here N = 1000, can set to any number

# Assign random dates within time range of data
dates <- sample(seq(min(as.Date(master$seg_id)), max(as.Date(master$seg_id)), by="day"), 1842, replace = T)
rpts$dt <- dates
absence=rpts %>% as.data.frame() %>% mutate(presabs=0) %>% rename(lon=x,lat=y)
presence=master %>% select(c(seg_id,lat,lon)) %>% mutate(dt=as.Date(seg_id)) %>% mutate(presabs=1) %>% select(-c(seg_id))
df=rbind(presence,absence)

### extract ####
source("/Users/heatherwelch/Dropbox/Eco-ROMS/heather_working/Eco-ROMS-private/Extracto_Scripts/Extracto_ROMS.R",chdir = TRUE)
getvarROMS <- function(nc,varname,inpts,desired.resolution,FUN,name){
  if ((desired.resolution*10) %% 2 ==0) stop("Desired Resolution must be an odd number (e.g. 0.3 not 0.2)")
  #inpts$dt <- as.POSIXct(inpts$dt, '%Y-%m-%d',tz='UTC')
  nc.data <- nc_open(nc, write=FALSE)
  lat <- ncvar_get(nc.data,'lat'); lat <- lat[1,]
  lon <- ncvar_get(nc.data,'lon'); lon <- lon[,1]
  nrows <- length(lat); ncols <- length(lon)
  yr <- ncvar_get(nc.data,'year'); mth <- ncvar_get(nc.data,'month'); day <- ncvar_get(nc.data,'day')
  tim <- as.POSIXct(paste(yr,mth,day,sep='-'),tz='UTC') %>% as.Date()
  desired.resolution = desired.resolution/2
  for (i in 1:nrow(inpts)){
    print(paste(varname,inpts$dt[i],sep=' '))
    if (inpts$dt[i] %in% tim){
      xdate <- which(inpts$dt[i]==tim)
      c <- which.min(abs(lon-inpts$lon[i]))
      c_low <- which.min(abs(lon-(inpts$lon[i]-desired.resolution)))
      c_up <- which.min(abs(lon-(inpts$lon[i]+desired.resolution)))
      r <- which.min(abs(lat-inpts$lat[i]))
      r_low <- which.min(abs(lat-(inpts$lat[i]-desired.resolution)))
      r_up <- which.min(abs(lat-(inpts$lat[i]+desired.resolution)))
      numcols=abs(c_up-c_low); numrows=abs(r_up-r_low)
      
      if (desired.resolution*2!=0.1){
        data.var  <-  ncvar_get(nc.data,varname,start=c(c_low,r_low,xdate),
                                count=c(numcols,numrows,1),verbose=FALSE)
        inpts[i,paste(varname,'_',name,sep='')] <- FUN(data.var[!is.nan(data.var)])
      } else{
        data.var.point  <-  ncvar_get(nc.data,varname,start=c(c,r,xdate),
                                      count=c(1,1,1),verbose=FALSE)
        inpts[i,paste(varname,'_',name,sep='')] <- data.var.point
       }
    } #else{
    #   inpts[i,paste(varname,'_',name,sep='')] <- NA ### commenting this out so that already extracted values are not overwritted with NAs
    # }
  }
  nc_close(nc.data)
  return(inpts)
}
input_file=df
# input_file$dt2 <- as.Date(input_file$dt, format="%Y-%m-%d") %>% format(.,"%m/%d/%y")
input_file$dt <- as.Date(input_file$dt, format="%Y-%m-%d") #%>% format(.,"%m/%d/%y")) 
# inpts=input_file

# input_file$dt2=as.POSIXct(input_file$dt, '%Y-%m-%d',tz='UTC')

input_file <- input_file[input_file$lat>=30 & input_file$lat<=48,] #remove points outside of ROMS boundary
ROMS_files_newNRT <- list.files("~/Dropbox/Eco-ROMS/ROMS & Bathym Data/WCNRT",pattern=".nc", full.names=T)
xtracto_NRT=function(input_file_newRT,netcdf_list){
  #Mean 0.1
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[1], 'BV', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[2], 'curl', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  #input_file_newRT <- getvarROMS(ROMS_files_newNRT[3], 'd26', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[4], 'ild', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[5], 'ssh', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[7], 'sst', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[8], 'su', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[9], 'sustr', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[10], 'sv', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[11], 'svstr', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  
  #SD0.3
  # input_file_newRT <- getvarROMS(ROMS_files_newNRT[7], 'sst', input_file_newRT, desired.resolution = 0.3, sd, 'sd_0.3')
  # input_file_newRT <- getvarROMS(ROMS_files_newNRT[5], 'ssh', input_file_newRT, desired.resolution = 0.3, sd, 'sd_0.3')
  
  #Now get bathymetry (z).
  #Mean 0.1 
  rastIn <- '~/Dropbox/Eco-ROMS/ROMS & Bathym Data/BATHYMETRY ETOPO1/z_.1.grd'
  input_file_newRT <- getZ(input_file_newRT,rastIn,name="z_0.1")
  
  #Now get bathymetry sd, known as zsd
  #SD0.3
  rastIn <- '~/Dropbox/Eco-ROMS/ROMS & Bathym Data/BATHYMETRY ETOPO1/zsd_.3.grd'
  input_file_newRT <- getZsd(input_file_newRT,rastIn,name="zsd_0.3")
  
  #Create EKE from surface geostrophic velocity fields
  # eke = (u2 + v2) / 2
  input_file_newRT$EKE_0.1 <- (input_file_newRT$su_mean_0.1^2 + input_file_newRT$sv_mean_0.1^2)/2
  
  #Create log EKE
  input_file_newRT$logEKE_0.1 <- log(input_file_newRT$EKE_0.1)
  
  #Get lunar illumination
  library(lunar)
  col_index <- length(input_file_newRT)+1
  counter=1
  for (i in 1:nrow(input_file_newRT)){
    ill <- lunar.illumination(input_file_newRT$dt[i])
    input_file_newRT[counter,col_index] <- ill
    counter=counter+1
  }
  colnames(input_file_newRT)[col_index] <- "lunar"
  
  #Get curl & w at additional resolutions
  # input_file_newRT <- getvarROMS(ROMS_files_newNRT[2], 'curl', input_file_newRT, desired.resolution = 0.5, mean, 'mean_0.5')
  
  return(input_file_newRT)
}
AIS_ROMS_NRT <- xtracto_NRT(input_file, ROMS_files_newNRT)

ROMS_files_newNRT <- list.files("~/Dropbox/Eco-ROMS/ROMS & Bathym Data/WCNRT_2017-2018",pattern="731.nc", full.names=T)
xtracto_NRT2=function(input_file_newRT,netcdf_list){
  #Mean 0.1
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[1], 'BV', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[2], 'curl', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  #input_file_newRT <- getvarROMS(ROMS_files_newNRT[3], 'd26', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[3], 'ild', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[4], 'ssh', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[5], 'sst', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[6], 'su', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[7], 'sustr', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[8], 'sv', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  input_file_newRT <- getvarROMS(ROMS_files_newNRT[9], 'svstr', input_file_newRT, desired.resolution = 0.1, mean, 'mean_0.1')
  
  #SD0.3
  # input_file_newRT <- getvarROMS(ROMS_files_newNRT[5], 'sst', input_file_newRT, desired.resolution = 0.3, sd, 'sd_0.3')
  # input_file_newRT <- getvarROMS(ROMS_files_newNRT[4], 'ssh', input_file_newRT, desired.resolution = 0.3, sd, 'sd_0.3')
  
  #Now get bathymetry (z).
  #Mean 0.1 
  rastIn <- '~/Dropbox/Eco-ROMS/ROMS & Bathym Data/BATHYMETRY ETOPO1/z_.1.grd'
  input_file_newRT <- getZ(input_file_newRT,rastIn,name="z_0.1")
  
  #Now get bathymetry sd, known as zsd
  #SD0.3
  rastIn <- '~/Dropbox/Eco-ROMS/ROMS & Bathym Data/BATHYMETRY ETOPO1/zsd_.3.grd'
  input_file_newRT <- getZsd(input_file_newRT,rastIn,name="zsd_0.3")
  
  #Create EKE from surface geostrophic velocity fields
  # eke = (u2 + v2) / 2
  input_file_newRT$EKE_0.1 <- (input_file_newRT$su_mean_0.1^2 + input_file_newRT$sv_mean_0.1^2)/2
  
  #Create log EKE
  input_file_newRT$logEKE_0.1 <- log(input_file_newRT$EKE_0.1)
  
  #Get lunar illumination
  library(lunar)
  col_index <- length(input_file_newRT)+1
  counter=1
  for (i in 1:nrow(input_file_newRT)){
    ill <- lunar.illumination(input_file_newRT$dt[i])
    input_file_newRT[counter,col_index] <- ill
    counter=counter+1
  }
  colnames(input_file_newRT)[col_index] <- "lunar"
  
  #Get curl & w at additional resolutions
  # input_file_newRT <- getvarROMS(ROMS_files_newNRT[2], 'curl', input_file_newRT, desired.resolution = 0.5, mean, 'mean_0.5')
  
  return(input_file_newRT)
}
AIS_ROMS_NRT2 <- xtracto_NRT2(AIS_ROMS_NRT, ROMS_files_newNRT)
# master=AIS_ROMS_NRT2[complete.cases(AIS_ROMS_NRT2),]

ROMS_files_newNRT <- list.files("~/Dropbox/Eco-ROMS/ROMS & Bathym Data/WCNRT_2017-2018",pattern="231.nc", full.names=T)
# a=AIS_ROMS_NRT2[c(1:66, 68:122, 125:278, 280:734, 736:757, 759:861, 863:989, 991:1060,1063:nrow(AIS_ROMS_NRT2)),]
# # a=a %>% filter(dt!=2018-09-21)
AIS_ROMS_NRT3 <- xtracto_NRT2(AIS_ROMS_NRT2, ROMS_files_newNRT)
master=AIS_ROMS_NRT3[complete.cases(AIS_ROMS_NRT3),]

colchange=function(x){
  c=colnames(x)
  a=gsub("presabs","PresAbs",c)%>%lapply(.,function(x)gsub("_mean_0.1","",x))%>%unlist()%>%lapply(.,function(x)gsub("_0.3","",x))%>%unlist()%>%lapply(.,function(x)gsub("_0.1","",x))%>%unlist()%>%gsub("zsd","z_sd",.)%>%gsub("EKE","unlogged_EKE",.)%>%gsub("logunlogged_EKE","EKE",.)%>%gsub("BV_frequency","bv",.)
  colnames(x)=a
  return(x)
}
a=colchange(master)

### fit ####
gbm.x <- c("curl","ild", "ssh", "sst", "su", "sustr", "sv", "svstr", 
           "z", "z_sd", "EKE", "lunar","BV")

brtt=gbm.fixed(data=a,gbm.x = gbm.x,gbm.y = "PresAbs",tree.complexity = 3, bag.fraction = 0.6, n.trees=1000)

gbm.x <- c("curl","ild", "ssh", "sst", "su", "sustr", "sv", "svstr", 
            "z_sd", "EKE", "lunar","BV")

brtt2=gbm.fixed(data=a,gbm.x = gbm.x,gbm.y = "PresAbs",tree.complexity = 3, bag.fraction = 0.6, n.trees=1000)

gbm.x <- c("curl","ild", "ssh", "sst", "su", "sustr", "sv", "svstr", 
           "z_sd", "EKE", "lunar","BV","lat")

brtt3=gbm.fixed(data=a,gbm.x = gbm.x,gbm.y = "PresAbs",tree.complexity = 3, bag.fraction = 0.6, n.trees=1000)


## predict ####
source("/Users/heatherwelch/Dropbox/Eco-ROMS/heather_working/Eco-ROMS-private/EcoROMs/create_ROMS_daily_stack.R",chdir = TRUE)
source("/Users/heatherwelch/Dropbox/Eco-ROMS/heather_working/Eco-ROMS-private/EcoROMs/create_ROMS_raster.R")
droppath="/Users/heatherwelch/Dropbox/"
predDir=paste0(droppath,"Eco-ROMS/ROMS & Bathym Data/daily_prediction_layers/")
staticDir=paste0(droppath,"Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/")
studyarea=readShapeSpatial(paste0(staticDir,"sa_square_coast3.shp"))
template= raster("/Users/heatherwelch/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")
get_date="2016-08-01"
daily_stack=create_ROMS_daily_stack(get_date = get_date,predDir = predDir,staticDir=staticDir,droppath=droppath,template=template)
stack=as.data.frame(daily_stack,stringsAsFactors=F) %>% rename(BV=bv)

### brt 1 2016 ####

brtt_pred=predict.gbm(brtt,newdata=stack,n.trees=1000,type='response')
brtt_pred <- setValues(template,brtt_pred)%>%mask(.,studyarea)
plot(brtt_pred)
plot(xy,add=T)

xy <- master[, c('lat', 'lon')]

layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
#layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
par(mar=c(4,4,1,1)) # I usually set my margins before each plot
pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
ncolors <- 100
breaks <- seq(0,1,,ncolors+1)
image(brtt_pred, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,48))
maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
points(xy,cex=4,pch = ".",col="black")
box()

par(mar=c(4,4,0,1)) # I usually set my margins before each plot
levs <- breaks[-1] - diff(breaks)/2
image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
mtext(paste0("Probability of gap presence "), side=1, line=2.5)

box()

### brt 2 2016 ####

brtt_pred=predict.gbm(brtt2,newdata=stack,n.trees=1000,type='response')
brtt_pred <- setValues(template,brtt_pred)%>%mask(.,studyarea)
plot(brtt_pred)
# plot(xy,add=T)
# 
# xy <- master[, c('lat', 'lon')] %>% mutate(layer=1)

layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
#layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
par(mar=c(4,4,1,1)) # I usually set my margins before each plot
pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
ncolors <- 100
breaks <- seq(0,1,,ncolors+1)
image(brtt_pred, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,48))
maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
points(xy,cex=4,pch = ".",col="black")
box()

par(mar=c(4,4,0,1)) # I usually set my margins before each plot
levs <- breaks[-1] - diff(breaks)/2
image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
mtext(paste0("Probability of gap presence "), side=1, line=2.5)

box()


### brt 3 2016 ####

brtt_pred=predict.gbm(brtt3,newdata=stack,n.trees=1000,type='response')
brtt_pred <- setValues(template,brtt_pred)%>%mask(.,studyarea)
plot(brtt_pred)
# plot(xy,add=T)
# 
# xy <- master[, c('lat', 'lon')] %>% mutate(layer=1)

layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
#layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
par(mar=c(4,4,1,1)) # I usually set my margins before each plot
pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
ncolors <- 100
breaks <- seq(0,1,,ncolors+1)
image(brtt_pred, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,48))
maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
points(xy,cex=4,pch = ".",col="black")
box()

par(mar=c(4,4,0,1)) # I usually set my margins before each plot
levs <- breaks[-1] - diff(breaks)/2
image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
mtext(paste0("Probability of gap presence "), side=1, line=2.5)

box()




get_date="2018-08-01"
daily_stack=create_ROMS_daily_stack(get_date = get_date,predDir = predDir,staticDir=staticDir,droppath=droppath,template=template)
stack=as.data.frame(daily_stack,stringsAsFactors=F) %>% rename(BV=bv)

### brt 1 2018 ####

brtt_pred=predict.gbm(brtt,newdata=stack,n.trees=1000,type='response')
brtt_pred <- setValues(template,brtt_pred)%>%mask(.,studyarea)
plot(brtt_pred)
# plot(xy,add=T)
# 
# xy <- master[, c('lat', 'lon')]

layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
#layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
par(mar=c(4,4,1,1)) # I usually set my margins before each plot
pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
ncolors <- 100
breaks <- seq(0,1,,ncolors+1)
image(brtt_pred, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,48))
maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
points(xy,cex=4,pch = ".",col="black")
box()

par(mar=c(4,4,0,1)) # I usually set my margins before each plot
levs <- breaks[-1] - diff(breaks)/2
image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
mtext(paste0("Probability of gap presence "), side=1, line=2.5)

box()

### brt 2 2018 ####

brtt_pred=predict.gbm(brtt2,newdata=stack,n.trees=1000,type='response')
brtt_pred <- setValues(template,brtt_pred)%>%mask(.,studyarea)
plot(brtt_pred)
# plot(xy,add=T)
# 
# xy <- master[, c('lat', 'lon')] %>% mutate(layer=1)

layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
#layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
par(mar=c(4,4,1,1)) # I usually set my margins before each plot
pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
ncolors <- 100
breaks <- seq(0,1,,ncolors+1)
image(brtt_pred, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,48))
maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
points(xy,cex=4,pch = ".",col="black")
box()

par(mar=c(4,4,0,1)) # I usually set my margins before each plot
levs <- breaks[-1] - diff(breaks)/2
image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
mtext(paste0("Probability of gap presence "), side=1, line=2.5)

box()


### brt 3 2018 ####

brtt_pred=predict.gbm(brtt3,newdata=stack,n.trees=1000,type='response')
brtt_pred <- setValues(template,brtt_pred)%>%mask(.,studyarea)
plot(brtt_pred)
# plot(xy,add=T)
# 
# xy <- master[, c('lat', 'lon')] %>% mutate(layer=1)

layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
#layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
par(mar=c(4,4,1,1)) # I usually set my margins before each plot
pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
ncolors <- 100
breaks <- seq(0,1,,ncolors+1)
image(brtt_pred, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,48))
maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
points(xy,cex=4,pch = ".",col="black")
box()

par(mar=c(4,4,0,1)) # I usually set my margins before each plot
levs <- breaks[-1] - diff(breaks)/2
image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
mtext(paste0("Probability of gap presence "), side=1, line=2.5)

box()



