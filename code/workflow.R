today_list<-seq(as.Date("2021-04-15"),as.Date("2021-04-15"),by=1)
for (d in 1:length(today_list)) {
  today<-today_list[d]
  
  forecast_project_id<-"UCSC_P_EDM" 
  forecast_model_id   <- "v1.4"
  forecast_iteration_id<- uuid::UUIDgenerate() 
  
  source("code/01 EFI utils.R")
  source("code/02 EFI settings.R")
  
  update<-T
  source("code/11 get NEON phenology data.R")
  update<-T
  source("code/12 get NEON weather data api.R")
  update<-T
  source("code/13 get NOAA weather data.R")
  
  source("code/21 preprocess data.R")
  source("code/22 prepare embeddings.R")
  source("code/23 train GP model.R")
  source("code/24 fill previous gaps.R")
  source("code/25 forecast.R")
  source("code/26 output for submission.R")
}



