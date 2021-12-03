# https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-usage

if (update) {
  ###temp
  # what did we have last time?
  path_last<-stringr::str_replace(path, as.character(today), as.character(today-1))
  if (!file.exists(paste0(path_last,"NEON_weather.csv"))) {
    neon_weather_df<-data.frame()
    time_last<-as.Date("1900-01-01")
  } else {
    neon_weather_df<-read_csv(paste0(path_last,"NEON_weather.csv"))
    time_last<-max(neon_weather_df$date[which(!is.na(neon_weather_df$tmean))])
  }
  release_last<-lubridate::floor_date(time_last, unit="month")
  
  # what release could we possibly have now?
  release_this<-lubridate::floor_date(today-15, unit="month")-months(1)
  
  # Request data availability info
  err<-T
  while (err) {
    tryCatch( {
      req.temp <- GET("http://data.neonscience.org/api/v0/products/DP1.00003.001")
      err<-F
    },
    error = function(e) {
      err<-T
      message(conditionMessage(e))
    })
  }
  
  # make this JSON readable
  avail.temp <- jsonlite::fromJSON(content(req.temp, as="text"), simplifyDataFrame=T, flatten=T)
  # get data availability list for the product
  temp.urls <- unlist(avail.temp$data$siteCodes$availableDataUrls)
  
  temp_df_all_sites<-vector(mode="list")
    for (i in 1:length(site_list))  {
               urls_all_time<-temp.urls[grep(site_list[i], temp.urls)]
               releases<-as.Date(paste0(stringr::str_sub(urls_all_time,-7,-1),"-01"))
               urls_new<-urls_all_time[which (releases>release_last&releases<=release_this)]
               if(length(urls_new)>0) {
                 temp_df_one_site<-vector(mode="list",length = length(urls_new))
                 for (j in 1:length(urls_new)) {
                   err<-T
                   while (err) {
                     tryCatch( {
                       tmp <- GET(urls_new[j])
                       err<-F
                     },
                     error = function(e) {
                       err<-T
                       message(conditionMessage(e))
                     })
                   }
                   tmp.files <- jsonlite::fromJSON(content(tmp, as="text"))
                   # tmp.files$data$files$name
                   if (length(grep("TAAT_30min", tmp.files$data$files$name))>0) {
                     temp_df <- read.delim(tmp.files$data$files$url
                                           [intersect(grep("TAAT_30min", 
                                                           tmp.files$data$files$name),
                                                      grep("basic", 
                                                           tmp.files$data$files$name))], 
                                           sep=",") %>%
                       dplyr::select(startDateTime,tempTripleMean) %>% 
                       mutate(date=as.Date(startDateTime)) %>%
                       dplyr::select(-startDateTime) %>% 
                       mutate(siteID=site_list[i])  %>% 
                       group_by(siteID, date) %>% 
                       summarize(tmin=min(tempTripleMean),
                                 tmax=max(tempTripleMean),
                                 tmean=mean(tempTripleMean)) %>% 
                       drop_na()
                   } else {
                     temp_df<-data.frame()
                   }
                   temp_df_one_site[[j]]<-temp_df
                   
                   print(paste0(i,",",j))
                 }
                 temp_df_one_site<-bind_rows(temp_df_one_site)
               } else {
                 temp_df_one_site<-NA
               }
               temp_df_all_sites[[i]]<-temp_df_one_site
             }
  temp_df_all_sites<-temp_df_all_sites[which(!is.na(temp_df_all_sites))]
  temp_df_all_sites<-bind_rows(temp_df_all_sites)
  if(nrow(temp_df_all_sites)==0) {
    temp_df_all_sites<-temp_df_all_sites %>% 
      mutate(siteID=character(0),
             date=as.Date(character(0)),
             tmin=double(0),
             tmax=double(0),
             tmean=double(0)
      )
  }
  print("temp")
  
  ###humi
  # what did we have last time?
  path_last<-stringr::str_replace(path, as.character(today), as.character(today-1))
  if (!file.exists(paste0(path_last,"NEON_weather.csv"))) {
    neon_weather_df<-data.frame()
    time_last<-as.Date("1900-01-01")
  } else {
    neon_weather_df<-read_csv(paste0(path_last,"NEON_weather.csv"))
    time_last<-max(neon_weather_df$date[which(!is.na(neon_weather_df$humi))])
  }
  release_last<-lubridate::floor_date(time_last, unit="month")
  
  # what release could we possibly have now?
  release_this<-lubridate::floor_date(today-15, unit="month")
  
  # Request data availability info
  err<-T
  while (err) {
    tryCatch( {
      req.humi <- GET("http://data.neonscience.org/api/v0/products/DP1.00098.001")
      err<-F
    },
    error = function(e) {
      err<-T
      message(conditionMessage(e))
    })
  }
  
  # make this JSON readable
  avail.humi <- jsonlite::fromJSON(content(req.humi, as="text"), simplifyDataFrame=T, flatten=T)
  # get data availability list for the product
  humi.urls <- unlist(avail.humi$data$siteCodes$availableDataUrls)
  
  humi_df_all_sites<-vector(mode="list")
    for (i in 1:length(site_list)) {
               urls_all_time<-humi.urls[grep(site_list[i], humi.urls)]
               releases<-as.Date(paste0(stringr::str_sub(urls_all_time,-7,-1),"-01"))
               urls_new<-urls_all_time[which (releases>release_last&releases<=release_this)]
               if(length(urls_new)>0) {
                 humi_df_one_site<-vector(mode="list")
                 for (j in 1:length(urls_new)) {
                   err<-T
                   while (err) {
                     tryCatch( {
                       tmp <- GET(urls_new[j])
                       err<-F
                     },
                     error = function(e) {
                       err<-T
                       message(conditionMessage(e))
                     })
                   }
                   tmp.files <- jsonlite::fromJSON(content(tmp, as="text"))
                   # tmp.files$data$files$name
                   if (length(grep("000.060.030.RH_30min", tmp.files$data$files$name))>0) {
                     humi_df <- read.delim(tmp.files$data$files$url
                                           [intersect(grep("000.060.030.RH_30min", 
                                                           tmp.files$data$files$name),
                                                      grep("basic", 
                                                           tmp.files$data$files$name))], 
                                           sep=",")%>%
                       dplyr::select(startDateTime,RHMean) %>% 
                       mutate(date=as.Date(startDateTime)) %>%
                       dplyr::select(-startDateTime) %>% 
                       mutate(siteID=site_list[i])  %>% 
                       group_by(siteID, date) %>% 
                       summarize(humi=mean(RHMean)) %>% 
                       drop_na()
                   } else {
                     humi_df<-data.frame()
                   }
                   
                   humi_df_one_site[[j]]<-humi_df
                   
                   print(paste0(i,",",j))
                 }
                 humi_df_one_site<-bind_rows(humi_df_one_site)
               } else {
                 humi_df_one_site<-NA
               }
               humi_df_all_sites[[i]]<-humi_df_one_site
             }
  humi_df_all_sites<-humi_df_all_sites[which(!is.na(humi_df_all_sites))]
  humi_df_all_sites<-bind_rows(humi_df_all_sites)
  if(nrow(humi_df_all_sites)==0) {
    humi_df_all_sites<-humi_df_all_sites %>% 
      mutate(siteID=character(0),
             date=as.Date(character(0)),
             humi=double(0)
      )
  }
  print("humi")
  
  ### prcp
  # what did we have last time?
  path_last<-stringr::str_replace(path, as.character(today), as.character(today-1))
  if (!file.exists(paste0(path_last,"NEON_weather.csv"))) {
    neon_weather_df<-data.frame()
    time_last<-as.Date("1900-01-01")
  } else {
    neon_weather_df<-read_csv(paste0(path_last,"NEON_weather.csv"))
    time_last<-max(neon_weather_df$date[which(!is.na(neon_weather_df$prcp))])
  }
  release_last<-lubridate::floor_date(time_last, unit="month")
  
  # what release could we possibly have now?
  release_this<-lubridate::floor_date(today-15, unit="month")
  
  # Request data availability info
  err<-T
  while (err) {
    tryCatch( {
      req.prcp <- GET("http://data.neonscience.org/api/v0/products/DP1.00006.001")
      err<-F
    },
    error = function(e) {
      err<-T
      message(conditionMessage(e))
    })
  }
  
  # make this JSON readable
  avail.prcp <- jsonlite::fromJSON(content(req.prcp, as="text"), simplifyDataFrame=T, flatten=T)
  # get data availability list for the product
  prcp.urls <- unlist(avail.prcp$data$siteCodes$availableDataUrls)
  prcp_df_all_sites<-vector(mode="list")
    for (i in 1:length(site_list)) {
               urls_all_time<-prcp.urls[grep(site_list[i], prcp.urls)]
               releases<-as.Date(paste0(stringr::str_sub(urls_all_time,-7,-1),"-01"))
               urls_new<-urls_all_time[which (releases>release_last&releases<=release_this)]
               if(length(urls_new)>0) {
                 prcp_df_one_site<-vector(mode="list")
                 for (j in 1:length(urls_new)) {
                   err<-T
                   while (err) {
                     tryCatch( {
                       tmp <- GET(urls_new[j])
                       err<-F
                     },
                     error = function(e) {
                       err<-T
                       message(conditionMessage(e))
                     })
                   }
                   tmp.files <- jsonlite::fromJSON(content(tmp, as="text"))
                   # tmp.files$data$files$name
                   if (length(grep("PRIPRE_30min", tmp.files$data$files$name))>0) {
                     prcp_df1 <- read.delim(tmp.files$data$files$url
                                            [intersect(grep("PRIPRE_30min", 
                                                            tmp.files$data$files$name),
                                                       grep("basic", 
                                                            tmp.files$data$files$name))], 
                                            sep=",") %>%
                       dplyr::select(startDateTime,priPrecipBulk) %>% 
                       mutate(date=as.Date(startDateTime)) %>%
                       dplyr::select(-startDateTime) %>% 
                       mutate(siteID=site_list[i])  %>% 
                       group_by(siteID, date) %>% 
                       summarize(prcp=sum(priPrecipBulk)) %>% 
                       drop_na()%>% 
                       mutate(group="primary")
                   } else {
                     prcp_df1<-data.frame()
                   }
                   if (length(grep("SECPRE_30min", tmp.files$data$files$name))>0) {
                     prcp_df2 <- read.delim(tmp.files$data$files$url
                                            [intersect(grep("SECPRE_30min", 
                                                            tmp.files$data$files$name),
                                                       grep("basic", 
                                                            tmp.files$data$files$name))], 
                                            sep=",")%>%
                       dplyr::select(startDateTime,secPrecipBulk) %>% 
                       mutate(date=as.Date(startDateTime)) %>%
                       dplyr::select(-startDateTime) %>% 
                       mutate(siteID=site_list[i])  %>% 
                       group_by(siteID, date) %>% 
                       summarize(prcp=sum(secPrecipBulk)) %>% 
                       drop_na()%>% 
                       mutate(group="secondary")
                   }else {
                     prcp_df2<-data.frame()
                   }
                   prcp_df<-bind_rows(prcp_df1,prcp_df2)
                   if (nrow(prcp_df)>0) {
                     prcp_df<-prcp_df%>% 
                       group_by(siteID, date) %>% 
                       arrange(desc(group)) %>% 
                       slice(1) %>% 
                       ungroup()
                   }
                   prcp_df_one_site[[j]]<-prcp_df
                   
                   print(paste0(i,",",j))
                 }
                 prcp_df_one_site<-bind_rows(prcp_df_one_site)
               } else {
                 prcp_df_one_site<-NA
               }
               prcp_df_all_sites[[i]]<-prcp_df_one_site
             }
  prcp_df_all_sites<-prcp_df_all_sites[which(!is.na(prcp_df_all_sites))]
  prcp_df_all_sites<-bind_rows(prcp_df_all_sites)
  if(nrow(prcp_df_all_sites)==0) {
    prcp_df_all_sites<-prcp_df_all_sites %>% 
      mutate(siteID=character(0),
             date=as.Date(character(0)),
             prcp=double(0),
             group=character(0)
             )
  }
  print("prcp")
  
  df_all_sites<-temp_df_all_sites %>% 
    full_join(humi_df_all_sites, by=c("siteID", "date")) %>% 
    full_join(prcp_df_all_sites, by=c("siteID", "date"))
  
  if (nrow(neon_weather_df)>0) {
    df_all_sites<-df_all_sites %>% 
      left_join(neon_weather_df, by=c("siteID", "date")) %>% 
      mutate(tmin=case_when(is.na(tmin.y)~tmin.y, TRUE ~tmin.x),
             tmax=case_when(is.na(tmax.y)~tmax.y, TRUE ~tmax.x),
             tmean=case_when(is.na(tmean.y)~tmean.y, TRUE ~tmean.x),
             humi=case_when(is.na(humi.y)~humi.y, TRUE ~humi.x),
             prcp=case_when(is.na(prcp.y)~prcp.y, TRUE ~prcp.x)) %>% # to make sure I keep the latest version
      dplyr::select(siteID, date, tmin, tmax, tmean, humi, prcp)
  }
  neon_weather_df_new<-bind_rows(neon_weather_df,df_all_sites) %>% 
    arrange(desc(row_number())) %>% # to make sure I keep the latest version
    distinct(siteID, date, .keep_all = T) %>% 
    arrange(siteID, date)
  
  write_csv(neon_weather_df_new,paste0(path,"NEON_weather.csv"))
}

neon_weather_df<-read_csv(paste0(path,"NEON_weather.csv"))
