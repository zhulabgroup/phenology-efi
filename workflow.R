today_list<-seq(as.Date("2021-02-12"),as.Date("2021-02-18"),by=1)
for (d in 1:length(today_list)) {
  today<-today_list[d]
  
  forecast_project_id<-"UCSC_P_EDM" 
  forecast_model_id   <- "v0.1"
  forecast_iteration_id<- uuid::UUIDgenerate() 
  
  source("EFI utils.R")
  source("EFI settings.R")
  
  update<-T
  source("get NEON phenology data.R")
  update<-F
  source("get NEON weather data.R")
  update<-T
  source("get NOAA weather data.R")
  
  source("preprocess data.R")
  source("prepare embeddings.R")
  source("train GP model.R")
  source("fill previous gaps.R")
  source("forecast.R")
  source("output for submission.R")
}



