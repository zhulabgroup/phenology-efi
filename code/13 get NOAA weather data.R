source("code/utilities/downloadNOAAFiles.R")
source("code/utilities/noaa_gefs_read.R")
coord_df<-read_csv(paste0(path,"coord.csv"))
site_list <- coord_df$siteID %>% unlist()

if(update) {
  cycle="00"
  base_dir <- paste0(getwd(),"/noaa/NOAAGEFS_1hr/")
  dir.create(base_dir, recursive = T)
  
  all_date_list<-seq(as.Date("2020-09-25"),today,by=1)
  if (length(list.files(base_dir))==0) {
    download_date_list<-all_date_list
  } else {
    downloaded_date_list<-list.files(paste0(base_dir,"BART"))
    download_date_list<-as.Date(setdiff(as.character(all_date_list),downloaded_date_list))
  }

  for(k in 1:length(site_list)){
    for(j in 1:length(download_date_list)){
      if (!dir.exists(paste0(base_dir,site_list[k],"/",download_date_list[j],"/",cycle))|
          length(list.files(paste0(base_dir,site_list[k],"/",download_date_list[j],"/",cycle)))<31) {
        try(download_noaa_files_s3(siteID = coord_df$siteID[k],
                                   date = download_date_list[j],
                                   cycle = cycle,
                                   local_directory =getwd() ))
      }
    }
  }
  
  # last
  path_last<-stringr::str_replace(path, as.character(today), as.character(today-1))
  if (!file.exists(paste0(path_last,"NOAA_weather.csv"))) {
    noaa_weather_df<-data.frame()
    time_last<-as.Date("2020-09-24")
  } else {
    noaa_weather_df<-read_csv(paste0(path_last,"NOAA_weather.csv"))
    time_last<-max(noaa_weather_df$forecast_date)
  }
  
  weather_list<-vector(mode="list")
  new_date_list<-seq(time_last+1,today,by=1)
  for (i in 1:length(new_date_list)){
    try(weather_list[[i]]<-noaa_gefs_read(base_dir, new_date_list[i], cycle, site_list) %>% 
          dplyr::select(siteID, time, 
                        temp=air_temperature,
                        humi=relative_humidity,
                        prcp=precipitation_flux) %>%
          mutate(temp=temp-273.15) %>% #how to transform humidity?
          mutate(prcp=prcp*60*60/30) %>% 
          mutate(date=as.Date(time)) %>% 
          dplyr::select(-time) %>% 
          group_by(siteID, date) %>% 
          summarize(tmax=max(temp),
                    tmin=min(temp),
                    tmean=mean(temp),
                    humi=mean(humi),
                    prcp=sum(prcp)) %>% 
          mutate(forecast_date=new_date_list[i]) %>% 
          ungroup()
    )
    print(new_date_list[i])
  }

  noaa_weather_df<-rbind(noaa_weather_df,data.table::rbindlist(weather_list) )%>% 
    group_by(siteID, date) %>% 
    arrange(desc(forecast_date)) %>% 
    slice(1) %>% 
    ungroup() 
  write_csv(noaa_weather_df,paste0(path,"NOAA_weather.csv"))


  neon_weather_df<-read_csv(paste0(path,"NEON_weather.csv")) %>% 
    mutate(group="NEON")
  noaa_weather_df<-read_csv(paste0(path,"NOAA_weather.csv")) %>% 
    mutate(group="NOAA") %>% 
    dplyr::select(-forecast_date)
  
  weather_df<-bind_rows(neon_weather_df,noaa_weather_df) %>% 
    gather(key="var", value="value", -siteID, -date, -group) %>%
    tidyr::complete(siteID, date, var,group)
  # group_by(siteID, date, var) %>% 
  # arrange(group) %>% 
  # slice(1) %>% 
  # ungroup()

  p1<-ggplot(weather_df %>% filter(var=="tmax"))+
    geom_line(aes(x=date,y=value, col=group, group=group), alpha=0.5)+
    theme_classic()+
    facet_wrap(~siteID,nrow=2)
  p2<-ggplot(weather_df %>% filter(var=="tmin"))+
    geom_line(aes(x=date,y=value, col=group, group=group), alpha=0.5)+
    theme_classic()+
    facet_wrap(~siteID,nrow=2)
  p3<-ggplot(weather_df %>% filter(var=="tmean"))+
    geom_line(aes(x=date,y=value, col=group, group=group), alpha=0.5)+
    theme_classic()+
    facet_wrap(~siteID,nrow=2)+
    ylab("temperature (degree Celsius)")
  p4<-ggplot(weather_df %>% filter(var=="prcp"))+
    geom_line(aes(x=date,y=value, col=group, group=group), alpha=0.5)+
    theme_classic()+
    facet_wrap(~siteID,nrow=2)+
    ylab("precipitation (mm)")
  
  cairo_pdf(paste0(path,"weather.pdf"))
  print(gridExtra::grid.arrange(p3,p4,ncol=1))
  dev.off()
  
  write_csv(weather_df,paste0(path,"weather.csv"))
}

weather_df<-read_csv(paste0(path, "weather.csv")) %>% 
  spread(key="group", value="value") %>% 
  mutate(value=case_when(
    is.na(NEON) ~ NOAA,
    !is.na(NEON) ~ NEON
  )) %>% 
  dplyr::select(-NEON, -NOAA)