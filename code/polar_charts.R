# roseplot

library(ggplot2)
a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")
master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),] %>% filter(flag!="") %>% mutate(month=month(timestamp,label=T))
ggplot(data=master,aes(x=month,y=log(gap_hours),fill=vessel_class))+
  geom_bar(stat="identity",width=1,size=0.1)+
  xlab("")+ylab("")#+coord_polar()+theme_dark()

install.packages("highcharter")
library(magrittr)
library(highcharter)

data("weather")

x <- c("Min", "Mean", "Max")
y <- sprintf("{point.%s}", c("min_temperaturec", "mean_temperaturec", "max_temperaturec"))
tltip <- tooltip_table(x, y)

hchart(weather, type = "columnrange",
       hcaes(x = date, low = min_temperaturec, high = max_temperaturec,
             color = mean_temperaturec)) %>%
  hc_chart(polar = TRUE) %>%
  hc_yAxis( max = 30, min = -10, labels = list(format = "{value} C"),
            showFirstLabel = FALSE) %>%
  hc_xAxis(
    title = list(text = ""), gridLineWidth = 0.5,
    labels = list(format = "{value: %b}")) #%>%
  hc_tooltip(useHTML = TRUE, pointFormat = tltip,
             headerFormat = as.character(tags$small("{point.x:%d %B, %Y}")))
  
  ####
  a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")
  master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),] %>% filter(flag!="") %>% mutate(ymd=ymd_hms(timestamp)) %>% mutate(ymd=date(ymd)) %>% 
    group_by(ymd) %>% summarise(min=min(gap_km),max=max(gap_km),mean=mean(gap_km))

  hchart(master, type = "columnrange",
         hcaes(x = ymd, low = log(min), high = log(max),
               color = log(mean))) %>%
    hc_chart(polar = TRUE) %>%
    hc_yAxis( max = log(10182102), min = 0, labels = list(format = "{value} hours"),
              showFirstLabel = FALSE) %>%
    hc_xAxis(
      title = list(text = ""), gridLineWidth = 0.5,
      labels = list(format = "{value: %b}"))%>%
  hc_tooltip(useHTML = TRUE, pointFormat = tltip,
             headerFormat = as.character(tags$small("{point.x:%d %B, %Y}")))
  
  #####
  a=read.csv("/Volumes/SeaGate/IUU_GRW/data/gaps_2018_v20190618.csv")
  master=a %>% filter(on_fishing_list_best=="True") %>% mutate(seg_id=substr(seg_id,11,20))%>% .[complete.cases(.),] %>% filter(flag!="") %>% mutate(ymd=ymd_hms(timestamp)) %>% mutate(ymd=date(ymd)) %>% 
    group_by(ymd) %>% summarise(min=min(gap_hours),max=max(gap_hours),mean=mean(gap_hours))
  
  x <- c("Min", "Mean", "Max")
  y <- sprintf("{point.%s}", c("min", "mean", "max"))
  tltip <- tooltip_table(x, y)
  
  hchart(master, type = "columnrange",
         hcaes(x = ymd, low = log(min), high = log(max),
               color = log(mean))) %>%
    hc_chart(polar = TRUE) %>%
    hc_yAxis( max = log(6872.64), min = 0, labels = list(format = "{value} hours"),
              showFirstLabel = FALSE) %>%
    hc_xAxis(
      title = list(text = ""), gridLineWidth = 0.5,
      labels = list(format = "{value: %b}"))%>%
    hc_tooltip(useHTML = TRUE, pointFormat = tltip,
               headerFormat = as.character(tags$small("{point.x:%d %B, %Y}")))
  ######
  hchart(master, type = "column",
         hcaes(x = month, low = min, high = max,
               color = mean)) %>%
    hc_chart(polar = TRUE) %>%
    hc_yAxis( max = 10182102, min = 0, labels = list(format = "{value} hours"),
              showFirstLabel = FALSE) %>%
    hc_xAxis(
      title = list(text = ""), gridLineWidth = 0.5,
      labels = list(format = "{value}"))%>%
    hc_tooltip(useHTML = TRUE, pointFormat = tltip,
               headerFormat = as.character(tags$small("{point.x:%d %B, %Y}")))
  