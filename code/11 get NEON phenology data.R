library(tidyverse)

source("code/utilities/downloadPhenoCam.R")
source("code/utilities/calculatePhenoCamUncertainty.R")
if (update) {
    ##Selected Sites for Challenge
  siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
               "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
               "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")
  
  site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
  
  allData <- data.frame(matrix(nrow = 0, ncol = 5))
  
  message(paste0("Downloading and generating phenology targets ", Sys.time()))
  
  for(i in 1:length(siteIDs)){
    siteName <- siteIDs[i]
    message(siteName)
    if(siteName != "NEON.D11.CLBJ.DP1.00033"){
      URL_gcc90 <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_1day.csv",sep="") ##URL for daily summary statistics
      URL_individual <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_roistats.csv",sep="") ##URL for individual image metrics
    }else{
      URL_gcc90 <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_2000_1day.csv",sep="") ##URL for daily summary statistics
      URL_individual <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_2000_roistats.csv",sep="") ##URL for individual image metrics
    }
    phenoData <- download.phenocam(URL = URL_gcc90)
    dates <- unique(phenoData$date)
    phenoData_individual <- download.phenocam(URL=URL_individual,skipNum = 17)
    gcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates, var="gcc_90") ##Calculates standard deviations on daily gcc90 values
    rcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates, var="rcc_90") ##Calculates standard deviations on daily rcc90 values
    
    subPhenoData <- phenoData %>% 
      mutate(siteID = stringr::str_sub(siteName, 10, 13), 
             time = date) %>% 
      dplyr::select(time, siteID, gcc_90, rcc_90) %>% 
      cbind(gcc_sd, rcc_sd)
    
    allData <- rbind(allData,subPhenoData)
    
  }
  
  full_time <- seq(min(allData$time),max(allData$time), by = "1 day")
  
  full_time <- tibble(time = rep(full_time, 8),
                      siteID = c(rep("HARV", length(full_time)),
                                 rep("BART", length(full_time)),
                                 rep("SCBI", length(full_time)),
                                 rep("STEI", length(full_time)),
                                 rep("UKFS", length(full_time)),
                                 rep("GRSM", length(full_time)),
                                 rep("DELA", length(full_time)),
                                 rep("CLBJ", length(full_time))))
  
  
  allData <- left_join(full_time, allData, by = c("time", "siteID"))
  write_csv(allData, paste0(path, "phenology-targets.csv"))
}

####

ts_all<-read_csv(paste0(path,"phenology-targets.csv")) %>% 
  tidyr::complete(siteID,time) %>% 
  filter(time<today)

if (update) {
  p<-ggplot(ts_all)+
    geom_line(aes(x=time, y=rcc_90), col="orange")+
    geom_ribbon(aes(x=time, ymin=rcc_90-1.96*rcc_sd, ymax=rcc_90+1.96*rcc_sd), fill="orange", alpha=0.5)+
    geom_line(aes(x=time, y=gcc_90), col="darkgreen")+
    geom_ribbon(aes(x=time, ymin=gcc_90-1.96*gcc_sd, ymax=gcc_90+1.96*gcc_sd), fill="darkgreen", alpha=0.5)+
    theme_classic()+
    facet_wrap(~siteID)+
    guides(col=F, fill=F)
  cairo_pdf(paste0(path, "phenology-targets.pdf"))
  print(p)
  dev.off()
}

if(update) {
  rois_coord <- get_rois() %>% 
    distinct(site, lat, lon)
  
  site_list<-ts_all%>% distinct(siteID) %>% unlist()
  coord_list<-vector(mode="list")
  for (i in 1:length(site_list)) {
    coord_list[[i]]<-rois_coord %>% 
      filter(stringr::str_detect(rois_coord$site,site_list[i])) %>% 
      sample_n(1) %>% 
      mutate(siteID=site_list[i])
  }
  coord_df<-rbindlist(coord_list) %>% 
    dplyr::select(-site) %>% 
    mutate(id=row_number())
  # ggplot(as.data.frame(coord_df))+
  #   geom_point(aes(x=lon, y=lat), alpha=1)+
  #   theme_minimal()
  
  write_csv(as.data.frame(coord_df), paste0(path, "coord.csv"))
} else {
  coord_df<-read_csv(paste0(path, "coord.csv"))
}

if (update) {
  distMat<-GeoDistanceInMetresMatrix(coord_df)/1000 
  distMat<-distMat/1000
  write_csv(as.data.frame(distMat),paste0(path, "distMat.csv"))
} else {
  distMat<-as.matrix(read_csv(paste0(path, "distMat.csv")))
}

