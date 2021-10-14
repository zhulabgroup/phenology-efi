# devtools::install_github("NEONScience/NEON-utilities/neonUtilities") # need to install Bioconductor

# https://www.neonscience.org/resources/learning-hub/tutorials/neondatastackr
library(neonUtilities)
library(tidyverse)

base_path<-getwd()
if (update) {
  dir.create(paste0(base_path,"/NEON weather"))
  coord_df<-read_csv(paste0(path,"coord.csv"))
  site_list <- coord_df$siteID %>% unlist()
  
  # triple-aspirated air temperature
  trip.temp <- loadByProduct(dpID="DP1.00003.001", 
                             site=site_list,
                             package="basic",
                             check.size = F)
  # names(trip.temp)
  # write_csv(trip.temp$TAAT_30min,paste0(base_path,"NEON weather/NEON_temp_raw.csv"))
  
  temp_df<-trip.temp$TAAT_30min %>% 
    dplyr::select(siteID, startDateTime,tempTripleMean) %>% 
    mutate(date=as.Date(startDateTime)) %>% 
    dplyr::select(-startDateTime) %>% 
    group_by(siteID, date) %>% 
    summarize(tmin=min(tempTripleMean),
              tmax=max(tempTripleMean),
              tmean=mean(tempTripleMean)) %>% 
    drop_na()
  write_csv(temp_df,paste0(base_path,"/NEON weather/NEON_temp.csv"))
  
  ###
  rel.hum <- loadByProduct(dpID="DP1.00098.001", 
                           site=site_list,
                           package="basic",
                           check.size = F)
  # names(rel.hum)
  # write_csv(rel.hum$RH_30min,paste0(base_path,"NEON weather/NEON_humi_raw.csv"))
  
  humi_df<-rel.hum$RH_30min %>% 
    dplyr::select(siteID, startDateTime,RHMean) %>% 
    mutate(date=as.Date(startDateTime)) %>% 
    dplyr::select(-startDateTime) %>% 
    group_by(siteID, date) %>% 
    summarize(humi=mean(RHMean)) %>% 
    mutate(humi=humi/100) %>% 
    drop_na()
  write_csv(humi_df,paste0(base_path,"/NEON weather/NEON_humi.csv"))
  
  
  
  ### precipitation
  prcp <- loadByProduct(dpID="DP1.00006.001", 
                        site=site_list,
                        package="basic",
                        check.size = F)
  # names(prcp)
  # write_csv(prcp$PRIPRE_30min,paste0(base_path,"NEON weather/NEON_prcp_primary_raw.csv"))
  # write_csv(prcp$SECPRE_30min,paste0(base_path,"NEON weather/NEON_prcp_primary_raw.csv"))
  
  prcp_df1<-prcp$PRIPRE_30min %>% 
    dplyr::select(siteID, startDateTime,priPrecipBulk) %>% 
    mutate(date=as.Date(startDateTime)) %>% 
    dplyr::select(-startDateTime) %>% 
    group_by(siteID, date) %>% 
    summarize(prcp=sum(priPrecipBulk)) %>% 
    drop_na() %>% 
    mutate(group="primary")
  
  prcp_df2<-prcp$SECPRE_30min %>% 
    # filter(!(secPrecipBulk==0&secPrecipExpUncert==0)) %>% 
    dplyr::select(siteID, startDateTime,secPrecipBulk) %>% 
    mutate(date=as.Date(startDateTime)) %>% 
    dplyr::select(-startDateTime) %>% 
    group_by(siteID, date) %>% 
    summarize(prcp=sum(secPrecipBulk)) %>% 
    drop_na() %>% 
    mutate(group="secondary")
  
  prcp_df<-rbind(prcp_df1,prcp_df2) %>% 
    group_by(siteID, date) %>% 
    arrange(desc(group)) %>% 
    slice(1) %>% 
    ungroup()
  
  write_csv(prcp_df,paste0(base_path,"/NEON weather/NEON_prcp.csv"))
  
  neon_weather_df<-temp_df %>% 
    full_join(humi_df, by=c("siteID", "date")) %>% 
    full_join(prcp_df, by=c("siteID", "date")) 
  
  write_csv(neon_weather_df,paste0(path,"NEON_weather.csv"))
}

if (!update) {
  if(!file.exists(paste0(path,"NEON_weather.csv"))) {
    path_last<-stringr::str_replace(path, as.character(today), as.character(today-1))
    neon_weather_df<-read_csv(paste0(path_last,"NEON_weather.csv"))
    write_csv(neon_weather_df,paste0(path,"NEON_weather.csv"))
  } else {
    neon_weather_df<-read_csv(paste0(path,"NEON_weather.csv"))
  }
}
